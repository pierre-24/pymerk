import sys
from typing import TextIO
from pymerk.ensemble import Ensemble
from pymerk.driver import BaseDriver, XtbDriver


class BaseFilter:
    """Base class for filters"""

    def __init__(self, driver: BaseDriver):
        self.driver = driver

    def filter(self, ensemble: Ensemble, output: TextIO = sys.stdout) -> Ensemble:
        raise NotImplementedError()


class BaseEnergyFilter(BaseFilter):
    def __init__(self, driver: BaseDriver, ethr: float, label: str = 'ΔE'):
        super().__init__(driver)
        self.ethr = ethr
        self.label = label

    def _compute_total_energy(self, geometry, output: TextIO) -> float:
        """Subclasses implement specific energy/correction logic here."""
        raise NotImplementedError

    def filter(self, ensemble: Ensemble, output: TextIO = sys.stdout) -> Ensemble:
        print(f'* Filtering ({self.label} < {self.ethr} a.u.)')

        # 1. Calculate and assign energies
        for i, geometry in enumerate(ensemble.elements, 1):
            print(f'> Computing energy of molecule #{i}')
            energy = self._compute_total_energy(geometry, output)
            geometry.energy = energy
            print(f'  .. {energy:.8f} a.u.')

        # 2. Perform the threshold filtering
        min_energy = min(x.energy for x in ensemble.elements)

        print(f'* Final relative {self.label} of conformer(s):')
        for i, geometry in enumerate(ensemble.elements):
            rel_e = geometry.energy - min_energy
            mark = '*' if rel_e < self.ethr else ''
            print(f'{i:5} {rel_e:.8f} {mark}')

        filtered = ensemble.filter(lambda x: x.energy - min_energy < self.ethr)
        print(f'* Done, retained {len(filtered)} conformer(s)')
        return filtered


class EnergyFilter(BaseEnergyFilter):
    """Filter on the energy (E)"""
    def _compute_total_energy(self, geometry, output: TextIO) -> float:
        return self.driver.get_energy(geometry, False, output=output)


class EnergyWithXtbGsolvFilter(BaseEnergyFilter):
    """Filter on the energy + gsolv (g*), the latter computed with `xtb`"""
    def __init__(self, driver: BaseDriver, xtb_driver: XtbDriver, ethr: float):
        super().__init__(driver, ethr, label='Δg*')
        self.xtb_driver = xtb_driver

    def _compute_total_energy(self, geometry, output: TextIO) -> float:
        # If the main driver is already XTB, don't do double work
        if self.driver is self.xtb_driver:
            return self.xtb_driver.get_energy(geometry, False, output)

        energy = self.driver.get_energy(geometry, False, output)
        xtb_gsolv = self.xtb_driver.get_gsolv(geometry, output)
        return energy + xtb_gsolv


class GibbsFreeEnergyWithXtbFilter(BaseEnergyFilter):
    def __init__(self, driver: BaseDriver, xtb_driver: XtbDriver, ethr: float, T: float = 298.15):
        super().__init__(driver, ethr, label=f'ΔG* @ {T}K')
        self.xtb_driver = xtb_driver
        self.T = T

    def _compute_total_energy(self, geometry, output: TextIO) -> float:
        if self.driver is self.xtb_driver:
            _, xtb_gibbs = self.xtb_driver.get_gibbs_free_energy(geometry, self.T, output)
            return xtb_gibbs

        energy = self.driver.get_energy(geometry, False, output)
        xtb_e, xtb_gibbs = self.xtb_driver.get_gibbs_free_energy(geometry, self.T, output)
        # Use XTB as a delta correction to the main driver energy
        return energy + (xtb_gibbs - xtb_e)


class GibbsFreeEnergyFilter(BaseEnergyFilter):
    def __init__(self, driver: BaseDriver, ethr: float, T: float = 298.15):
        super().__init__(driver, ethr, label=f'ΔG* @ {T}K')
        self.T = T

    def _compute_total_energy(self, geometry, output: TextIO) -> float:
        _, gibbs = self.driver.get_gibbs_free_energy(geometry, self.T, output)
        return gibbs
