import sys
from enum import Enum, auto
from typing import TextIO, Optional, Tuple

from pymerk.ensemble import Ensemble
from pymerk.driver import BaseDriver


class BaseFilter:
    """Base class for filters"""

    def __init__(self, driver: BaseDriver):
        self.main_driver = driver

    def filter(self, ensemble: Ensemble, output: TextIO = sys.stdout) -> Ensemble:
        raise NotImplementedError()


class SelectDriver(Enum):
    NONE = auto()
    MAIN = auto()
    AUX = auto()


class EnergyFilter(BaseFilter):
    """Filter on energy"""

    def __init__(
            self, driver: BaseDriver, ethr: float,
            gsolv_component: SelectDriver = SelectDriver.NONE,
            gtrv_component: SelectDriver = SelectDriver.NONE,
            aux_driver: Optional[BaseDriver] = None,
            label: str = 'E'
    ):
        super().__init__(driver)
        self.ethr = ethr
        self.gsolv_component = gsolv_component
        self.gtrv_component = gtrv_component
        self.aux_driver = aux_driver
        self.label = label

    def _get_components(
            self, driver: BaseDriver, geometry, use_solv: bool, use_gtrv: bool, T: float, output: TextIO
    ) -> Tuple[float, float, float]:
        """
        Helper to fetch E, Gsolv, and Gtrv from a driver based on requirements.
        Handles the logic of tuple unpacking for different driver methods.
        """
        if use_gtrv:
            # Returns (E, E+Gsolv, E+Gsolv+Gtrv) if use_solv=True else (E, E+Gtrv)
            res = driver.get_gibbs_free_energy(geometry, add_solvent=use_solv, T=T, output=output)
            e_elec = res[0]
            g_solv = (res[1] - res[0]) if use_solv else 0.0
            g_gtrv = (res[2] - res[1]) if use_solv else (res[1] - res[0])
        elif use_solv:
            # Returns (E, E+Gsolv)
            res = driver.get_energy(geometry, add_solvent=True, output=output)
            e_elec = res[0]
            g_solv = res[1] - res[0]
            g_gtrv = 0.0
        else:
            # Returns float (E)
            e_elec = driver.get_energy(geometry, add_solvent=False, output=output)
            g_solv = 0.0
            g_gtrv = 0.0

        return e_elec, g_solv, g_gtrv

    def _compute_total_energy(self, geometry, output: TextIO, T: float = 298.15) -> float:
        # 1. Determine requirements for the MAIN driver
        main_needs_solv = (self.gsolv_component == SelectDriver.MAIN)
        main_needs_gtrv = (self.gtrv_component == SelectDriver.MAIN)

        # MAIN always provides the Electronic Energy
        m_elec, m_gsolv, m_gtrv = self._get_components(
            self.main_driver, geometry, main_needs_solv, main_needs_gtrv, T, output
        )

        # 2. Determine requirements for the AUX driver (if needed)
        a_gsolv, a_gtrv = 0.0, 0.0
        aux_needs_solv = (self.gsolv_component == SelectDriver.AUX)
        aux_needs_gtrv = (self.gtrv_component == SelectDriver.AUX)

        if aux_needs_solv or aux_needs_gtrv:
            if self.aux_driver is None:
                raise RuntimeError('AUX driver required!')

            _, a_gsolv, a_gtrv = self._get_components(
                self.aux_driver, geometry, aux_needs_solv, aux_needs_gtrv, T, output
            )

        # 3. Sum the components based on the selected drivers
        total = m_elec

        # Add Solvation Correction
        if self.gsolv_component == SelectDriver.MAIN:
            total += m_gsolv
        elif self.gsolv_component == SelectDriver.AUX:
            total += a_gsolv

        # Add Thermal Correction (Gtrv)
        if self.gtrv_component == SelectDriver.MAIN:
            total += m_gtrv
        elif self.gtrv_component == SelectDriver.AUX:
            total += a_gtrv

        return total

    def filter(self, ensemble: Ensemble, output: TextIO = sys.stdout) -> Ensemble:
        print(f'* Filtering on {self.label} (threshold is {self.ethr:.6f} a.u.)')
        print(f'* elec={SelectDriver.MAIN}; gsolv={self.gsolv_component}; & gtrv={self.gtrv_component}')
        print(f'  Using MAIN={self.main_driver}')
        if self.aux_driver is not None:
            print(f'       & AUX={self.aux_driver}')

        # 1. Calculate and assign energies
        for i, geometry in enumerate(ensemble.elements, 1):
            print(f'> Computing energy of molecule #{i}', flush=True)
            energy = self._compute_total_energy(geometry, output)
            geometry.energy = energy
            print(f'  .. {energy:.8f} a.u.')

        # 2. Perform the threshold filtering
        min_energy = min(x.energy for x in ensemble.elements)

        print(f'* Final Δ{self.label} (w.r.t more stable) of conformer(s):')
        for i, geometry in enumerate(ensemble.elements):
            rel_e = geometry.energy - min_energy
            mark = '*' if rel_e < self.ethr else ''
            print(f'{i:5} {rel_e:.8f} {mark}')

        filtered = ensemble.filter(lambda x: x.energy - min_energy < self.ethr)
        print(f'* Done, retained {len(filtered)} conformer(s)', flush=True)
        return filtered
