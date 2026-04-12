import pytest
from pymerk.driver import XtbDriver
import shutil
import rmsd


@pytest.fixture
def xtb_driver(tmpdir):
    XTB_DRIVER = XtbDriver(tmpdir, shutil.which('xtb'))
    XTB_DRIVER.solvatation_model = 'alpb'
    XTB_DRIVER.solvent = 'thf'

    return XTB_DRIVER


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_get_energy(Ca_THF3_ensemble, xtb_driver):
    output_file = xtb_driver.workdir / 'output.log'
    with output_file.open('w') as f:
        assert xtb_driver.get_energy(Ca_THF3_ensemble.elements[0][0], f) == pytest.approx(-50.02, abs=1e-2)


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_get_gibbs(Ca_THF3_ensemble, xtb_driver):
    output_file = xtb_driver.workdir / 'output.log'
    with output_file.open('w') as f:
        energy, gibbs_energy = xtb_driver.get_gibbs_free_energy(Ca_THF3_ensemble.elements[0][0], output=f)
        assert energy == pytest.approx(-50.02, abs=1e-2)
        assert gibbs_energy == pytest.approx(-49.72, abs=1e-2)


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_opt(Ca_THF3_ensemble, xtb_driver):
    disp = .1
    old_geometry, old_energy = Ca_THF3_ensemble.elements[0]
    output_file = xtb_driver.workdir / 'output.log'

    # only one cycle
    old_geometry.positions[4, :] += disp
    with output_file.open('w') as f:
        new_geometry, new_energy = xtb_driver.optimize_geometry(old_geometry, f, maxcycle=1)

    old_geometry.positions[4, :] -= disp
    assert rmsd.rmsd(new_geometry.positions, old_geometry.positions) > .01

    # all cycles
    old_geometry.positions[4, :] += disp
    with output_file.open('w') as f:
        new_geometry, new_energy = xtb_driver.optimize_geometry(old_geometry, f)

    old_geometry.positions[4, :] -= disp
    assert rmsd.rmsd(new_geometry.positions, old_geometry.positions) < .01

    # final check
    assert new_geometry.charge == 2
    assert new_energy == pytest.approx(-50.02, abs=1e-2)
