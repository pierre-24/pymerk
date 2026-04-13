import sys

from pymerk.ensemble import Ensemble
from pymerk.driver import BaseDriver, XtbDriver

from typing import TextIO


class BaseFilter:
    """Base class for filters"""

    def __init__(self, driver: BaseDriver):
        self.driver = driver

    def filter(self, ensemble: Ensemble, output: TextIO = sys.stdout) -> Ensemble:
        raise NotImplementedError()


class EnergyFilter(BaseFilter):
    def __init__(self, driver: BaseDriver, ethr: float):
        super().__init__(driver)
        self.ethr = ethr

    def filter(self, ensemble: Ensemble, output: TextIO = sys.stdout) -> Ensemble:
        output.write('> Filtering (ΔE < {} a.u.)\n'.format(self.ethr))
        filtered_ensemble = Ensemble([])

        i = 0
        for geometry, _ in ensemble.elements:
            output.write('> Computing electronic energy of molecule #{}\n'.format(i + 1))
            energy = self.driver.get_energy(geometry, output=output)
            filtered_ensemble.elements.append((geometry, energy))
            i += 1

        min_energy = min(x[1] for x in filtered_ensemble.elements)
        filtered_ensemble.elements = list(filter(lambda x: x[1] - min_energy < self.ethr, filtered_ensemble.elements))

        return filtered_ensemble


class GibbsFreeEnergyWithXtbFilter(BaseFilter):
    def __init__(self, driver: BaseDriver, xtb_driver: XtbDriver, ethr: float, T: float = 298.15):
        super().__init__(driver)
        self.xtb_driver = xtb_driver
        self.ethr = ethr
        self.T = T

    def filter(self, ensemble: Ensemble, output: TextIO = sys.stdout) -> Ensemble:
        output.write('> Filtering (ΔG* < {} a.u.)\n'.format(self.ethr))

        filtered_ensemble = Ensemble([])

        i = 0
        for geometry, _ in ensemble.elements:
            output.write('> Computing Gibbs free energy of molecule #{}\n'.format(i + 1))

            if isinstance(self.driver, XtbDriver):
                energy, xtb_gibbs_free_energy = self.xtb_driver.get_gibbs_free_energy(geometry, self.T, output)
                xtb_energy = energy
            else:
                energy = self.driver.get_energy(geometry, output)
                xtb_energy, xtb_gibbs_free_energy = self.xtb_driver.get_gibbs_free_energy(geometry, self.T, output)

            filtered_ensemble.elements.append((geometry, energy + xtb_gibbs_free_energy - xtb_energy))
            i += 1

        min_energy = min(x[1] for x in filtered_ensemble.elements)
        print([x[1] - min_energy for x in filtered_ensemble.elements])
        filtered_ensemble.elements = list(filter(lambda x: x[1] - min_energy < self.ethr, filtered_ensemble.elements))

        return filtered_ensemble
