import pytest
import numpy
from pymerk.driver import XtbDriver
import shutil

XTB_DRIVER = None

if shutil.which('xtb'):
    XTB_DRIVER = XtbDriver('xtb')
    XTB_DRIVER.solvatation_model = 'alpb'
    XTB_DRIVER.solvent = 'thf'


@pytest.mark.skipif(XTB_DRIVER is None, reason='xtb driver not available')
def test_get_energy(Ca_THF3_ensemble):
    assert XTB_DRIVER.get_energy(Ca_THF3_ensemble.elements[0][0]) == pytest.approx(-50.02, abs=1e-2)


@pytest.mark.skipif(XTB_DRIVER is None, reason='xtb driver not available')
def test_get_gibbs(Ca_THF3_ensemble):
    assert XTB_DRIVER.get_gibbs_free_energy(Ca_THF3_ensemble.elements[0][0])[1] == pytest.approx(-49.72, abs=1e-2)


@pytest.mark.skipif(XTB_DRIVER is None, reason='xtb driver not available')
def test_opt(Ca_THF3_ensemble):
    old_geometry, old_energy = Ca_THF3_ensemble.elements[0]
    new_geometry, new_energy = XTB_DRIVER.optimize_geometry(old_geometry)

    assert new_geometry.charge == 2
    assert numpy.allclose(new_geometry.positions, old_geometry.positions, atol=1e-2)
    assert new_energy == pytest.approx(-50.02, abs=1e-2)
