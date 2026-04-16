import sys
import numpy as np
import rmsd
from enum import Enum, auto

from typing import TextIO, Optional

from pymerk.ensemble import Ensemble
from pymerk.driver import BaseDriver
from pymerk.molecule import Molecule


class SelectDriver(Enum):
    """Enumeration for selecting which driver computes a correction term.
    
    Used to route solvation and thermal corrections to main or auxiliary driver.
    
    Attributes:
        NONE: No correction is included.
        MAIN: Correction computed by main driver.
        AUX: Correction computed by auxiliary driver.
    """
    NONE = auto()
    MAIN = auto()
    AUX = auto()


class BaseFilter:
    """Abstract base class for ensemble filtering operations.
    
    Provides shared infrastructure for computing energies, reporting results,
    and filtering ensembles based on various criteria. Supports dual-driver evaluation
    where solvation and thermal corrections can come from different sources.
    
    Attributes:
        main_driver: Primary `BaseDriver` for computing electronic energies.
        aux_driver: Optional secondary `BaseDriver` for auxiliary corrections.
        label: Energy label for reporting (e.g., 'E', 'G'). Defaults to 'E'.
    """

    def __init__(self, driver: BaseDriver, aux_driver: Optional[BaseDriver] = None, label: str = 'E'):
        """Initialize a BaseFilter.
        
        Args:
            driver: Primary `BaseDriver` instance.
            aux_driver: Optional secondary `BaseDriver` for auxiliary corrections. Defaults to `None`.
            label: Label for energy reporting. Defaults to 'E'.
        """
        self.main_driver = driver
        self.aux_driver = aux_driver
        self.label = label

    @staticmethod
    def _get_components(
            driver: BaseDriver, geometry, use_solv: bool, use_gtrv: bool, T: float, output: TextIO
    ) -> tuple[float, float, float]:
        """Extract energy components from driver call.
        
        Routes calls to `get_energy()` or `get_gibbs_free_energy()` based on requested
        corrections and unpacks results into (electronic_energy, solvation_energy, thermal_rrho_energy).
        
        Args:
            driver: `BaseDriver` to query.
            geometry: `Molecule` to evaluate.
            use_solv: If `True`, compute solvation corrections.
            use_gtrv: If `True`, compute Gibbs free energy and thermal corrections.
            T: Temperature in Kelvin.
            output: File object for driver output.
            
        Returns:
            Tuple of (electronic_energy, g_solvation, g_thermal_rrho) in Hartree.
            Components not requested are returned as 0.0.
        """
        if use_gtrv:
            res = driver.get_gibbs_free_energy(geometry, add_solvent=use_solv, T=T, output=output)
            e_elec = res[0]
            g_solv = (res[1] - res[0]) if use_solv else 0.0
            g_gtrv = (res[2] - res[1]) if use_solv else (res[1] - res[0])
        elif use_solv:
            res = driver.get_energy(geometry, add_solvent=True, output=output)
            e_elec, g_solv, g_gtrv = res[0], res[1] - res[0], 0.0
        else:
            e_elec = driver.get_energy(geometry, add_solvent=False, output=output)
            g_solv, g_gtrv = 0.0, 0.0

        return e_elec, g_solv, g_gtrv

    def _compute_total_energy(
            self, geometry: Molecule, output: TextIO, T: float,
            gsolv: SelectDriver, gtrv: SelectDriver, skip_main: bool = False
    ) -> float:
        """Assemble total energy from main and/or auxiliary drivers.
        
        Computes electronic energy plus selected solvation and thermal corrections.
        Allows flexible routing of corrections to different drivers.
        
        Args:
            geometry: `Molecule` to evaluate.
            output: File object for driver output.
            T: Temperature in Kelvin.
            gsolv: Which driver provides solvation correction (`NONE`, `MAIN`, or `AUX`).
            gtrv: Which driver provides thermal correction (`NONE`, `MAIN`, or `AUX`).
            skip_main: If `True`, use cached `geometry.energy` instead of main driver. Defaults to `False`.
            
        Returns:
            Total energy in Hartree as `E_elec + G_solv + G_rrho`.
            
        Raises:
            RuntimeError: If `AUX` driver is required but not provided.
        """
        if not skip_main:
            m_e, m_s, m_t = self._get_components(
                self.main_driver, geometry, gsolv == SelectDriver.MAIN, gtrv == SelectDriver.MAIN, T, output)
        else:
            m_e, m_s, m_t = geometry.energy, 0, 0

        # 2. Aux Driver
        a_s, a_t = 0.0, 0.0
        if SelectDriver.AUX in (gsolv, gtrv):
            if not self.aux_driver:
                raise RuntimeError('AUX driver required but not provided!')

            _, a_s, a_t = self._get_components(
                self.aux_driver, geometry,
                gsolv == SelectDriver.AUX, gtrv == SelectDriver.AUX, T, output)

        return m_e + (m_s if gsolv == SelectDriver.MAIN else a_s) + (m_t if gtrv == SelectDriver.MAIN else a_t)

    def _report_results(
            self, old_ensemble: Ensemble, final_ensemble: Ensemble, ethr: float, check_convergence: bool = False):
        """Print standardized filtering results and RMSD matrix.
        
        Outputs relative energies of all input geometries (marking retained ones with '*'),
        number of retained conformers, and pairwise RMSD matrix.
        
        Args:
            old_ensemble: Original ensemble before filtering.
            final_ensemble: Filtered ensemble (retained geometries).
            ethr: Energy threshold used for filtering (in Hartree).
            check_convergence: If `True`, mark non-converged structures as excluded. Defaults to `False`.
        """
        if not final_ensemble.elements:
            print('! No conformers retained.')
            return

        min_e = min(x.energy for x in old_ensemble.elements)
        print('\n* Final values:')
        print(f"{'':15} {self.label:^14} {f'Δ{self.label}':^12}")
        for i, geom in enumerate(old_ensemble.elements):
            rel = geom.energy - min_e
            mark = '*' if rel < ethr else ''
            if check_convergence and not geom.converged:
                mark = ''
            print(f'{geom.name:15} {geom.energy:14.8f} {rel:12.8f} {mark}')

        print(f'\n* Done, retained {len(final_ensemble)} conformer(s)', flush=True)

        print(f'\n* RMSD matrix (Å) for {len(final_ensemble)} structures:')
        header = ' '.join(f'..{final_ensemble.elements[i].name[-4:]:4}' for i in range(len(final_ensemble)))
        print(f"{' ':16}{header}")
        for i, g1 in enumerate(final_ensemble.elements):
            row = [f'{rmsd.kabsch_rmsd(g1.positions, g2.positions):6.3f}' for j, g2 in
                   enumerate(final_ensemble.elements[:i + 1])]
            print(f"{g1.name:15} {' '.join(row)}")

    def filter(self, ensemble: Ensemble, output: TextIO = sys.stdout, T: float = 298.15) -> Ensemble:
        """Filter ensemble by a criterion implemented in subclasses.
        
        Args:
            ensemble: Input `Ensemble` to filter.
            output: File object for driver output. Defaults to `sys.stdout`.
            T: Temperature in Kelvin. Defaults to 298.15.
            
        Returns:
            Filtered `Ensemble` with retained geometries.
            
        Raises:
            NotImplementedError: Must be implemented by subclasses.
        """
        raise NotImplementedError()


class EnergyFilter(BaseFilter):
    """Filter ensemble by energy threshold using single-point calculations.
    
    Evaluates all geometries with specified drivers and removes those above
    energy threshold relative to the minimum.
    
    Attributes:
        ethr: Energy threshold in Hartree.
        gsolv: `SelectDriver` for solvation corrections.
        gtrv: `SelectDriver` for thermal corrections.
    """

    def __init__(
            self, driver: BaseDriver, ethr: float,
            gsolv_component: SelectDriver = SelectDriver.NONE, gtrv_component: SelectDriver = SelectDriver.NONE,
            aux_driver=None,
            label='E'
    ):
        """Initialize an EnergyFilter.
        
        Args:
            driver: Primary `BaseDriver` for energy calculations.
            ethr: Energy threshold in Hartree above minimum.
            gsolv_component: Which driver provides solvation energy. Defaults to `SelectDriver.NONE`.
            gtrv_component: Which driver provides thermal energy. Defaults to `SelectDriver.NONE`.
            aux_driver: Optional auxiliary `BaseDriver` for correction terms. Defaults to `None`.
            label: Label for reporting. Defaults to 'E'.
        """
        super().__init__(driver, aux_driver, label)
        self.ethr, self.gsolv, self.gtrv = ethr, gsolv_component, gtrv_component

    def filter(self, ensemble: Ensemble, output: TextIO = sys.stdout, T: float = 298.15) -> Ensemble:
        print(f'* Filtering on {self.label} (threshold: {self.ethr:.6f} a.u.)')

        print('* Setup Summary:')
        print(f'  - Components: elec=MAIN, gsolv={self.gsolv.name}, gtrv={self.gtrv.name}')
        print(f'  - Main driver: {self.main_driver}')
        if self.aux_driver:
            print(f'  - Aux driver:  {self.aux_driver}', flush=True)

        for i, geom in enumerate(ensemble.elements, 1):
            geom.energy = self._compute_total_energy(geom, output, T, self.gsolv, self.gtrv)
            geom.converged = True
            print(f'> {geom.name}: {geom.energy:.8f} a.u.', flush=True)

        min_energy = min(x.energy for x in ensemble.elements)
        filtered = ensemble.filter(lambda x: (x.energy - min_energy) < self.ethr)
        self._report_results(ensemble, filtered, self.ethr)
        return filtered


class OptFilter(BaseFilter):
    """Filter ensemble by geometry optimization and energy threshold.
    
    Optimizes all geometries and retains only those that converged and are
    below energy threshold.
    
    Attributes:
        ethr: Energy threshold in Hartree.
        use_solvent: If `True`, optimize with solvation model.
        gtrv: `SelectDriver` for thermal corrections.
        maxcycles: Maximum optimization cycles.
        optlevel: Optimization level ('loose', 'normal', 'tight').
    """

    def __init__(
            self, driver: BaseDriver, ethr: float, use_solvent: bool = True,
            gtrv_component: SelectDriver = SelectDriver.NONE,
            aux_driver=None,
            maxcycles: int = -1, optlevel: str = 'normal',
            label='E'
    ):
        """Initialize an OptFilter.
        
        Args:
            driver: Primary `BaseDriver` for optimization and energy.
            ethr: Energy threshold in Hartree.
            use_solvent: If `True`, optimize with solvation. Defaults to `True`.
            gtrv_component: Which driver provides thermal corrections. Defaults to `SelectDriver.NONE`.
            aux_driver: Optional auxiliary `BaseDriver`. Defaults to `None`.
            maxcycles: Maximum optimization cycles (-1 = default). Defaults to -1.
            optlevel: Optimization level. Defaults to 'normal'.
            label: Label for reporting. Defaults to 'E'.
        """
        super().__init__(driver, aux_driver, label)
        self.ethr, self.use_solvent, self.gtrv, self.maxcycles = ethr, use_solvent, gtrv_component, maxcycles
        self.optlevel = optlevel

    def filter(self, ensemble: Ensemble, output: TextIO = sys.stdout, T: float = 298.15) -> Ensemble:
        new_ensemble = Ensemble([x.copy() for x in ensemble.elements])

        print(f'* Optimization + Filter on {self.label} (threshold: {self.ethr:.6f} a.u.)')
        s = SelectDriver.MAIN if self.use_solvent else SelectDriver.NONE
        print('* Setup Summary:')
        print(f'  - Components: elec=MAIN, gsolv={s.name}, gtrv={self.gtrv.name}')
        print(f'  - Main driver: {self.main_driver}')
        if self.aux_driver:
            print(f'  - Aux driver:  {self.aux_driver}', flush=True)

        for i, geom in enumerate(new_ensemble.elements, 1):
            print(f'> Optimizing {geom.name}...', end=' ', flush=True)
            optimized = self.main_driver.optimize_geometry(
                geom, self.use_solvent, output, maxcycle=self.maxcycles, optlevel=self.optlevel)
            optimized.energy = self._compute_total_energy(
                optimized, output, T,
                SelectDriver.MAIN if self.use_solvent else SelectDriver.NONE, self.gtrv, skip_main=True)

            new_ensemble.elements[i - 1] = optimized
            print(f"Done. E = {optimized.energy:.8f} a.u. {'[FAILED]' if not optimized.converged else ''}")

        min_e = min(x.energy for x in new_ensemble.elements)
        filtered = new_ensemble.filter(lambda x: (x.energy - min_e) < self.ethr and x.converged)
        self._report_results(new_ensemble, filtered, self.ethr, check_convergence=True)
        return filtered


class MacroOptFilter(BaseFilter):
    """Filter ensemble by iterative macro-optimization.
    
    Performs multiple optimization cycles per structure with early discard logic.
    Useful for difficult optimizations or exhaustive refinement.
    
    Attributes:
        ethr: Energy threshold in Hartree for discard.
        use_solvent: If `True`, optimize with solvation.
        gtrv: `SelectDriver` for thermal corrections.
        maxcycles: Maximum total optimization iterations.
        optcycles: Optimization steps per macro-cycle.
        gradthr: Gradient threshold for early discard in Hartree/Bohr.
        optlevel: Optimization level.
    """

    def __init__(
            self, driver: BaseDriver, ethr: float, use_solvent: bool = True,
            gtrv_component: SelectDriver = SelectDriver.NONE, aux_driver=None,
            maxcycles: int = -1, optcycles: int = 10, gradthr: float = 1e-2, optlevel: str = 'normal',
            label: str = 'E'
    ):
        """Initialize a MacroOptFilter.
        
        Args:
            driver: Primary `BaseDriver` for optimization.
            ethr: Energy threshold in Hartree.
            use_solvent: If `True`, optimize with solvation. Defaults to `True`.
            gtrv_component: Which driver provides thermal corrections. Defaults to `SelectDriver.NONE`.
            aux_driver: Optional auxiliary `BaseDriver`. Defaults to `None`.
            maxcycles: Maximum total optimization cycles (-1 = unlimited). Defaults to -1.
            optcycles: Optimization steps per macro-cycle. Defaults to 10.
            gradthr: Gradient norm threshold for early discard in Hartree/Bohr. Defaults to 1e-2.
            optlevel: Optimization level. Defaults to 'normal'.
            label: Label for reporting. Defaults to 'E'.
        """
        super().__init__(driver, aux_driver, label)
        self.ethr, self.use_solvent, self.gtrv = ethr, use_solvent, gtrv_component
        self.maxcycles, self.optcycles, self.gradthr = maxcycles, optcycles, gradthr
        self.optlevel = optlevel

    def filter(self, ensemble: Ensemble, output: TextIO = sys.stdout, T: float = 298.15) -> Ensemble:
        print(f'* Macro-Optimization on {self.label} (threshold: {self.ethr:.6f} a.u.)')

        s = SelectDriver.MAIN if self.use_solvent else SelectDriver.NONE
        print('* Setup Summary:')
        print(f'  - Components: elec=MAIN, gsolv={s.name}, gtrv={self.gtrv.name}')
        print(f'  - Main driver: {self.main_driver}')
        if self.aux_driver:
            print(f'  - Aux driver:  {self.aux_driver}', flush=True)

        new_elements = [x.copy() for x in ensemble.elements]

        for geometry in new_elements:
            geometry.converged = False

        # Status: 0=active, 1=converged, 2=discarded
        status = [0] * len(new_elements)
        iteration = 0

        while 0 in status and (self.maxcycles < 0 or iteration < self.maxcycles):
            print(f'\n> Macrocycle {iteration // self.optcycles + 1}')
            for i, geom in enumerate(new_elements):
                if status[i] != 0:
                    continue

                print(f'  - Otimizing {geom.name}...', end=' ', flush=True)
                opt_geom = self.main_driver.optimize_geometry(
                    geom, self.use_solvent, output, maxcycle=self.optcycles, optlevel=self.optlevel)

                opt_geom.energy = self._compute_total_energy(
                    opt_geom, output, T,
                    SelectDriver.MAIN if self.use_solvent else SelectDriver.NONE, self.gtrv, skip_main=True)

                new_elements[i] = opt_geom
                print(f'E={opt_geom.energy:.8f}, grad={opt_geom.gnorm:.6f}', end=' ')

                if opt_geom.converged:
                    status[i] = 1
                    print('[CONVERGED]', flush=True)
                else:
                    print('[NOT CONVERGED]', flush=True)

            # Early discard logic
            min_e = min(x.energy for x in new_elements if status[list(new_elements).index(x)] != 2)
            for i, geom in enumerate(new_elements):
                if status[i] != 2 and geom.gnorm < self.gradthr and (geom.energy - min_e) > self.ethr:
                    status[i] = 2

            print(f'→ retained {len(list(filter(lambda x: x != 2, status)))} conformer(s)', flush=True)
            iteration += self.optcycles

        final_ensemble = Ensemble([new_elements[i] for i, s in enumerate(status) if s == 1])
        self._report_results(Ensemble(new_elements), final_ensemble, self.ethr, check_convergence=True)
        return final_ensemble


# Boltzmann constant in atomic units (Eh/K)
BOLTZMANN_CONSTANT_IN_AU = 3.166811563e-6


class BoltzmannFilter(BaseFilter):
    """Filter ensemble by Boltzmann population threshold.
    
    Computes Boltzmann populations at given temperature and retains conformers
    accounting for specified cumulative population threshold.
    
    Attributes:
        pthr: Cumulative population threshold (0-1).
        gsolv: `SelectDriver` for solvation corrections.
        gtrv: `SelectDriver` for thermal corrections.
    """

    def __init__(
            self,
            driver: BaseDriver,
            pthr: float,
            gsolv_component: Optional[SelectDriver] = None,
            gtrv_component: Optional[SelectDriver] = None,
            aux_driver: BaseDriver = None,
            label: str = 'E'
    ):
        """Initialize a BoltzmannFilter.
        
        Args:
            driver: Primary `BaseDriver` for energy calculations.
            pthr: Cumulative population threshold (e.g., 0.95 for 95%).
            gsolv_component: Which driver provides solvation energy. Defaults to `None`.
            gtrv_component: Which driver provides thermal energy. Defaults to `None`.
            aux_driver: Optional auxiliary `BaseDriver`. Defaults to `None`.
            label: Label for reporting. Defaults to 'E'.
        """
        super().__init__(driver, aux_driver, label)
        self.pthr = pthr
        self.gsolv = gsolv_component
        self.gtrv = gtrv_component

    def filter(self, ensemble: Ensemble, output: TextIO = sys.stdout, T: float = 298.15) -> Ensemble:
        print(f'* Boltzmann population filtering on Δ{self.label} (threshold: {self.pthr:.3f})')
        print('* Setup Summary:')
        print(f'  - Components:  elec=MAIN, gsolv={self.gsolv.name}, gtrv={self.gtrv.name}')
        print(f'  - Main driver: {self.main_driver}')
        if self.aux_driver:
            print(f'  - Aux driver:  {self.aux_driver}')

        # 1. Compute energies for all elements
        for i, geom in enumerate(ensemble.elements, 1):
            geom.energy = self._compute_total_energy(geom, output, T, self.gsolv, self.gtrv)
            print(f'> {geom.name}: {geom.energy:.8f} a.u.', flush=True)

        # 2. Calculate Boltzmann populations
        energies = np.array([x.energy for x in ensemble.elements])
        rel_energies = energies - np.min(energies)
        exp_factors = np.exp(-rel_energies / (BOLTZMANN_CONSTANT_IN_AU * T))
        populations = exp_factors / np.sum(exp_factors)

        # 3. Determine selection based on cumulative sum of sorted populations
        sorted_indices = np.argsort(populations)[::-1]
        cumulative_pop = np.cumsum(populations[sorted_indices])

        # Find how many conformers we need to keep to reach pthr
        cutoff_idx = np.searchsorted(cumulative_pop, self.pthr)
        keep_indices = sorted_indices[:cutoff_idx + 1]

        # 4. Final Reporting
        print('\n* Final Boltzmann population:')
        print(f"{' ':15} {self.label:^14} {f'Δ{self.label}':^12} {'Pop %':^8}")

        for i, (geom, pop) in enumerate(zip(ensemble.elements, populations)):
            rel_e = rel_energies[i]
            mark = '*' if i in keep_indices else ''
            print(f'{geom.name:15} {geom.energy:14.8f} {rel_e:12.8f} {pop * 100:6.1f} {mark}')

        # 5. Filter the ensemble using the set of valid indices
        keep_set = set(keep_indices)
        filtered_elements = [
            geom for idx, geom in enumerate(ensemble.elements) if idx in keep_set
        ]

        new_ensemble = Ensemble(filtered_elements)
        print(f'* Done, retained {len(new_ensemble)} conformer(s)')

        return new_ensemble
