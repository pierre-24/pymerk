import pytest
from pymerk.driver import XtbDriver
import shutil

XTB_DRIVER = None

if shutil.which('xtb'):
    XTB_DRIVER = XtbDriver('xtb')
    XTB_DRIVER.solvatation_model = 'alpb'
    XTB_DRIVER.solvent = 'thf'


@pytest.mark.skipif(XTB_DRIVER is None, reason='xtb driver not available')
def test_get_energy(Ca_THF3_ensemble):
    assert XTB_DRIVER.get_energy(Ca_THF3_ensemble.elements[0][0]) == pytest.approx(-50.024, abs=1e-3)
