import pytest
import shutil

from pymerk.driver import XtbDriver
from pymerk.scripts.filter import EnergyFilter, SelectDriver, OptFilter, MacroOptFilter, BoltzmannFilter


@pytest.fixture
def xtb_driver(tmpdir):
    XTB_DRIVER = XtbDriver(tmpdir, shutil.which('xtb'))
    XTB_DRIVER.solvatation_model = 'alpb'
    XTB_DRIVER.solvent = 'thf'

    return XTB_DRIVER


AU_TO_KCAL = 6.275030e2


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_filter_energy_with_xtb(Ca_THF2_ensemble, xtb_driver):
    # remove all conformers above 1 kcal/mol
    with (xtb_driver.workdir / 'output.log').open('w') as f:
        filt = EnergyFilter(xtb_driver, 1 / AU_TO_KCAL, label='E')
        new_ensemble = filt.filter(Ca_THF2_ensemble, f)
        assert len(new_ensemble) == len(Ca_THF2_ensemble) - 2


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_filter_gstar_with_xtb(Ca_THF2_ensemble, xtb_driver):
    # remove all conformers above 1 kcal/mol
    with (xtb_driver.workdir / 'output.log').open('w') as f:
        filt = EnergyFilter(xtb_driver, 1 / AU_TO_KCAL, gsolv_component=SelectDriver.MAIN, label='g*')
        new_ensemble = filt.filter(Ca_THF2_ensemble, f)
        assert len(new_ensemble) == len(Ca_THF2_ensemble) - 2


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_filter_gibbs_energy_with_xtb(Ca_THF2_ensemble, xtb_driver):
    # remove all conformers above 1 kcal/mol
    with (xtb_driver.workdir / 'output.log').open('w') as f:
        filt = EnergyFilter(
            xtb_driver, 1 / AU_TO_KCAL,
            gsolv_component=SelectDriver.MAIN, gtrv_component=SelectDriver.MAIN, label='G*')
        new_ensemble = filt.filter(Ca_THF2_ensemble, f)
        assert len(new_ensemble) == len(Ca_THF2_ensemble) - 1


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_filter_gibbs_energy_with_xtb_and_aux(Ca_THF2_ensemble, xtb_driver):
    # remove all conformers above 1 kcal/mol
    with (xtb_driver.workdir / 'output.log').open('w') as f:
        filt = EnergyFilter(
            xtb_driver, 1 / AU_TO_KCAL,
            gsolv_component=SelectDriver.MAIN, gtrv_component=SelectDriver.AUX, label='G*',
            aux_driver=xtb_driver
        )
        new_ensemble = filt.filter(Ca_THF2_ensemble, f)
        assert len(new_ensemble) == len(Ca_THF2_ensemble) - 2


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_filter_opt_with_xtb(Ca_THF2_ensemble, xtb_driver):
    # remove all conformers above 1.0 kcal/mol, if any
    with (xtb_driver.workdir / 'output.log').open('w') as f:
        xtb_driver.solvent = 'water'

        filt = OptFilter(
            xtb_driver, 1.0 / AU_TO_KCAL,
            True, gtrv_component=SelectDriver.AUX, aux_driver=xtb_driver, label='G*')
        new_ensemble = filt.filter(Ca_THF2_ensemble, f)
        assert len(new_ensemble) == len(Ca_THF2_ensemble)


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_filter_macroopt_with_xtb(Ca_THF2_ensemble, xtb_driver):
    # remove all conformers above 1.0 kcal/mol, if any
    with (xtb_driver.workdir / 'output.log').open('w') as f:
        xtb_driver.solvent = 'water'

        filt = MacroOptFilter(
            xtb_driver, 1.0 / AU_TO_KCAL,
            True, gtrv_component=SelectDriver.AUX, aux_driver=xtb_driver, label='G*')
        new_ensemble = filt.filter(Ca_THF2_ensemble, f)
        assert len(new_ensemble) < 4


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_filter_population_with_xtb(Ca_THF2_ensemble, xtb_driver):
    with (xtb_driver.workdir / 'output.log').open('w') as f:
        filt = BoltzmannFilter(
            xtb_driver, .95,
            gsolv_component=SelectDriver.MAIN, gtrv_component=SelectDriver.MAIN, label='G*')
        new_ensemble = filt.filter(Ca_THF2_ensemble, f)
        assert len(new_ensemble) == len(Ca_THF2_ensemble) - 1
