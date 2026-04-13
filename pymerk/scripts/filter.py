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


class BaseEnergyFilter(BaseFilter):
    def __init__(self, driver: BaseDriver, ethr: float):
        super().__init__(driver)
        self.ethr = ethr

    def _filter(self, ensemble: Ensemble) -> Ensemble:
        min_energy = min(x[1] for x in ensemble.elements)

        print('* Final relative energy of conformer(s):')
        for i, (_, energy) in enumerate(ensemble.elements):
            print('{:5} {:.8f} {}'.format(i, energy - min_energy, '*' if (energy - min_energy) < self.ethr else ''))

        return ensemble.filter(lambda x: x[1] - min_energy < self.ethr)


class EnergyFilter(BaseEnergyFilter):
    def filter(self, ensemble: Ensemble, output: TextIO = sys.stdout) -> Ensemble:
        print('* Filtering (ΔE < {} a.u.) with {}'.format(self.ethr, self.driver))
        filtered_ensemble = Ensemble([])

        i = 0
        for geometry, _ in ensemble.elements:
            print('> Computing electronic energy of molecule #{}'.format(i + 1))
            energy = self.driver.get_energy(geometry, output=output)
            filtered_ensemble.elements.append((geometry, energy))

            print('  .. {} a.u.'.format(energy))

            i += 1

        filtered_ensemble = self._filter(filtered_ensemble)
        print('* Done, retained {} conformer(s)'.format(len(filtered_ensemble)))

        return filtered_ensemble


class EnergyWithXtbGsolvFilter(BaseEnergyFilter):
    def __init__(self, driver: BaseDriver, xtb_driver: XtbDriver, ethr: float):
        super().__init__(driver, ethr)
        self.xtb_driver = xtb_driver

    def filter(self, ensemble: Ensemble, output: TextIO = sys.stdout) -> Ensemble:
        print('* Filtering (ΔE* < {} a.u.)'.format(self.ethr))

        filtered_ensemble = Ensemble([])

        i = 0
        for geometry, _ in ensemble.elements:
            print('> Computing g* of molecule #{}'.format(i + 1))

            if isinstance(self.driver, XtbDriver):
                total_energy = self.xtb_driver.get_energy(geometry, output)
            else:
                energy = self.driver.get_energy(geometry, output)
                xtb_gsolv = self.xtb_driver.get_gsolv(geometry, output)
                total_energy = energy + xtb_gsolv

            print('  .. {} a.u.'.format(total_energy))

            filtered_ensemble.elements.append((geometry, total_energy))
            i += 1

        filtered_ensemble = self._filter(filtered_ensemble)
        print('* Done, retained {} conformer(s)'.format(len(filtered_ensemble)))

        return filtered_ensemble


class GibbsFreeEnergyWithXtbFilter(BaseEnergyFilter):
    def __init__(self, driver: BaseDriver, xtb_driver: XtbDriver, ethr: float, T: float = 298.15):
        super().__init__(driver, ethr)
        self.xtb_driver = xtb_driver
        self.T = T

    def filter(self, ensemble: Ensemble, output: TextIO = sys.stdout) -> Ensemble:
        print('* Filtering (ΔG* < {} a.u. at T={})'.format(self.ethr, self.T))

        filtered_ensemble = Ensemble([])

        i = 0
        for geometry, _ in ensemble.elements:
            print('> Computing Gibbs free energy of molecule #{}'.format(i + 1))

            if isinstance(self.driver, XtbDriver):
                _, xtb_gibbs_free_energy = self.xtb_driver.get_gibbs_free_energy(geometry, self.T, output)
                total_energy = xtb_gibbs_free_energy
            else:
                energy = self.driver.get_energy(geometry, output)
                xtb_energy, xtb_gibbs_free_energy = self.xtb_driver.get_gibbs_free_energy(geometry, self.T, output)
                total_energy = energy + xtb_gibbs_free_energy - xtb_energy

            print('  .. {} a.u.'.format(total_energy))

            filtered_ensemble.elements.append((geometry, total_energy))
            i += 1

        filtered_ensemble = self._filter(filtered_ensemble)
        print('* Done, retained {} conformer(s)'.format(len(filtered_ensemble)))

        return filtered_ensemble
