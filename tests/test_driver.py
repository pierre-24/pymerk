import pytest
from pymerk.driver import XtbDriver, VlxDriver
import shutil
import rmsd


@pytest.fixture
def xtb_driver(tmpdir):
    XTB_DRIVER = XtbDriver(tmpdir, shutil.which('xtb'))
    XTB_DRIVER.solvatation_model = 'alpb'
    XTB_DRIVER.solvent = 'thf'

    return XTB_DRIVER


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_xtb_get_energy(Ca_THF2_ensemble, xtb_driver):
    output_file = xtb_driver.workdir / 'output.log'
    with output_file.open('w') as f:
        assert xtb_driver.get_energy(Ca_THF2_ensemble.elements[0][0], f) == pytest.approx(-33.247, abs=1e-2)


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_xtb_get_gsolv(Ca_THF2_ensemble, xtb_driver):
    output_file = xtb_driver.workdir / 'output.log'
    with output_file.open('w') as f:
        assert xtb_driver.get_gsolv(Ca_THF2_ensemble.elements[0][0], f) == pytest.approx(-0.244, abs=1e-2)


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_xtb_get_gibbs(Ca_THF2_ensemble, xtb_driver):
    output_file = xtb_driver.workdir / 'output.log'
    with output_file.open('w') as f:
        energy, gibbs_energy = xtb_driver.get_gibbs_free_energy(Ca_THF2_ensemble.elements[0][0], output=f)
        assert energy == pytest.approx(-33.247, abs=1e-2)
        assert gibbs_energy == pytest.approx(-33.061, abs=1e-2)


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_xtb_opt(Ca_THF2_ensemble, xtb_driver):
    disp = .05
    atom_index = 0

    old_geometry, old_energy = Ca_THF2_ensemble.elements[0]
    old_positions = old_geometry.positions.copy()
    output_file = xtb_driver.workdir / 'output.log'

    # only one cycle
    old_geometry.positions[atom_index, :] += disp
    modified_positions = old_geometry.positions.copy()

    with output_file.open('w') as f:
        new_geometry, new_energy = xtb_driver.optimize_geometry(old_geometry, f, maxcycle=1)

    assert rmsd.kabsch_rmsd(new_geometry.positions, old_positions) < rmsd.kabsch_rmsd(modified_positions, old_positions)

    # all cycles
    with output_file.open('w') as f:
        new_geometry, new_energy = xtb_driver.optimize_geometry(old_geometry, f)

    assert rmsd.kabsch_rmsd(new_geometry.positions, old_positions) < .01

    # final check
    assert new_geometry.charge == 2
    assert new_energy == pytest.approx(-33.24, abs=1e-2)


@pytest.fixture
def vlx_driver(tmpdir):
    VLX_DRIVER = VlxDriver(tmpdir, shutil.which('vlx'), 'b3lyp', 'sto-3g')
    VLX_DRIVER.solvatation_model = 'cpcm'
    VLX_DRIVER.solvent = 7.54

    return VLX_DRIVER


@pytest.mark.skipif(not shutil.which('vlx'), reason='vlx driver not available')
def test_vlx_get_energy(Ca_THF2_ensemble, vlx_driver):
    output_file = vlx_driver.workdir / 'output.log'
    with output_file.open('w') as f:
        assert vlx_driver.get_energy(Ca_THF2_ensemble.elements[0][0], f) == pytest.approx(-1129.626, abs=1e-2)


@pytest.mark.skipif(not shutil.which('vlx'), reason='vlx driver not available')
def test_vlx_opt(Ca_THF2_ensemble, vlx_driver):
    old_geometry, old_energy = Ca_THF2_ensemble.elements[0]
    output_file = vlx_driver.workdir / 'output.log'

    # only one cycle
    with output_file.open('w') as f:
        new_geometry, new_energy = vlx_driver.optimize_geometry(old_geometry, f, maxcycle=1)

    modified_positions = new_geometry.positions.copy()
    modfied_energy = new_energy

    assert rmsd.kabsch_rmsd(old_geometry.positions, modified_positions) > .01

    assert old_energy > new_energy

    # 2 cycles
    with output_file.open('w') as f:
        new_geometry, new_energy = vlx_driver.optimize_geometry(old_geometry, f, maxcycle=2)

    # final check
    assert rmsd.kabsch_rmsd(new_geometry.positions, old_geometry.positions) > .01
    assert rmsd.kabsch_rmsd(new_geometry.positions, modified_positions) > .01

    assert old_energy > new_energy
    assert modfied_energy > new_energy

    assert new_geometry.charge == 2
