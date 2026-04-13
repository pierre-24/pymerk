import pytest
import shutil

from pymerk.driver import XtbDriver
from pymerk.filter import EnergyFilter, GibbsFreeEnergyWithXtbFilter


@pytest.fixture
def xtb_driver(tmpdir):
    XTB_DRIVER = XtbDriver(tmpdir, shutil.which('xtb'))
    XTB_DRIVER.solvatation_model = 'alpb'
    XTB_DRIVER.solvent = 'thf'

    return XTB_DRIVER


AU_TO_KCAL = 6.275030e2


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_filter_energy(Ca_THF2_ensemble, xtb_driver):
    # remove all conformers above 1 kcal/mol
    new_ensemble = EnergyFilter(xtb_driver, 1 / AU_TO_KCAL).filter(Ca_THF2_ensemble)
    assert len(new_ensemble) == len(Ca_THF2_ensemble) - 2


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_filter_gibbs_energy_with_xtb(Ca_THF2_ensemble, xtb_driver):
    # remove all conformers above 1 kcal/mol
    new_ensemble = GibbsFreeEnergyWithXtbFilter(xtb_driver, xtb_driver, 1 / AU_TO_KCAL).filter(Ca_THF2_ensemble)
    assert len(new_ensemble) == len(Ca_THF2_ensemble) - 1
