import pytest
import shutil
import rmsd
import pathlib
import numpy as np

from pymerk.molecule import Molecule
from pymerk.driver import XtbDriver, VlxDriver, OrcaDriver


@pytest.fixture
def xtb_driver(tmpdir):
    XTB_DRIVER = XtbDriver(tmpdir, shutil.which('xtb'))
    XTB_DRIVER.solvatation_model = 'alpb'
    XTB_DRIVER.solvent = 'thf'

    return XTB_DRIVER


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_xtb_get_energy_wo_solvent(Ca_THF2_ensemble, xtb_driver):
    output_file = xtb_driver.workdir / 'output.log'
    with output_file.open('w') as f:
        assert xtb_driver.get_energy(Ca_THF2_ensemble.elements[0], False, f) == pytest.approx(-33.002, abs=1e-2)


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_xtb_get_energy_w_solvent(Ca_THF2_ensemble, xtb_driver):
    output_file = xtb_driver.workdir / 'output.log'
    with output_file.open('w') as f:
        assert xtb_driver.get_energy(Ca_THF2_ensemble.elements[0], True, f) == (
            pytest.approx(-33.002, abs=1e-2),
            pytest.approx(-33.247, abs=1e-2),
        )


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_xtb_get_gibbs(Ca_THF2_ensemble, xtb_driver):
    output_file = xtb_driver.workdir / 'output.log'
    with output_file.open('w') as f:
        elec_energy, elec_solv_energy, gibbs_energy = xtb_driver.get_gibbs_free_energy(
            Ca_THF2_ensemble.elements[0], add_solvent=True, output=f)

        assert elec_energy == pytest.approx(-33.002, abs=1e-2)
        assert elec_solv_energy == pytest.approx(-33.247, abs=1e-2)
        assert gibbs_energy == pytest.approx(-33.061, abs=1e-2)


def test_xtb_get_gibbs_from_log(tmpdir):
    cwd = pathlib.Path(__file__).parent
    drv = XtbDriver(tmpdir, '{} {}'.format(cwd / 'assets/run-dummy.sh', cwd / 'assets/gibbs.xtb.log'))

    total_energy, total_free_energy = drv.get_gibbs_free_energy(Molecule([], np.zeros((0, 3))))

    assert total_energy == pytest.approx(-114.63718369861, abs=1e-2)
    assert total_free_energy == pytest.approx(-113.902279812806, abs=1e-2)


@pytest.mark.skipif(not shutil.which('xtb'), reason='xtb driver not available')
def test_xtb_opt(Ca_THF2_ensemble, xtb_driver):
    disp = .05
    atom_index = 0

    old_geometry = Ca_THF2_ensemble.elements[0]
    old_positions = old_geometry.positions.copy()
    output_file = xtb_driver.workdir / 'output.log'

    # only one cycle
    old_geometry.positions[atom_index, :] += disp
    modified_positions = old_geometry.positions.copy()

    with output_file.open('w') as f:
        new_geometry = xtb_driver.optimize_geometry(old_geometry, True, f, maxcycle=1)

    assert rmsd.kabsch_rmsd(new_geometry.positions, old_positions) < rmsd.kabsch_rmsd(modified_positions, old_positions)
    assert new_geometry.gnorm > 5e-4
    assert not new_geometry.converged

    # all cycles
    with output_file.open('w') as f:
        new_geometry = xtb_driver.optimize_geometry(old_geometry, True, f)

    assert rmsd.kabsch_rmsd(new_geometry.positions, old_positions) < .01

    # final check
    assert new_geometry.charge == 2
    assert new_geometry.energy == pytest.approx(-33.24, abs=1e-2)
    assert new_geometry.gnorm < 5e-4
    assert new_geometry.converged


@pytest.fixture
def vlx_driver(tmpdir):
    VLX_DRIVER = VlxDriver(tmpdir, shutil.which('vlx'), 'b3lyp', 'sto-3g')
    VLX_DRIVER.solvatation_model = 'cpcm'
    VLX_DRIVER.solvent = 7.54

    return VLX_DRIVER


@pytest.mark.skipif(not shutil.which('vlx'), reason='vlx driver not available')
def test_vlx_get_energy(_2H2O_ensemble, vlx_driver):
    output_file = vlx_driver.workdir / 'output.log'
    with output_file.open('w') as f:
        assert vlx_driver.get_energy(_2H2O_ensemble.elements[0], True, f) == (
            pytest.approx(-150.628, abs=1e-2),
            pytest.approx(-150.636, abs=1e-2)
        )


@pytest.mark.skipif(not shutil.which('vlx'), reason='vlx driver not available')
def test_vlx_opt(_2H2O_ensemble, vlx_driver):
    old_geometry = _2H2O_ensemble.elements[0]
    output_file = vlx_driver.workdir / 'output.log'

    # only one cycle
    with output_file.open('w') as f:
        modified_geom_1 = vlx_driver.optimize_geometry(old_geometry, True, f, maxcycle=1)

    assert rmsd.kabsch_rmsd(old_geometry.positions, modified_geom_1.positions) > .01
    assert not modified_geom_1.converged

    # 2 cycles
    with output_file.open('w') as f:
        modified_geom_3 = vlx_driver.optimize_geometry(old_geometry, True, f, maxcycle=2)

    # final check
    assert rmsd.kabsch_rmsd(modified_geom_3.positions, old_geometry.positions) > .01
    assert rmsd.kabsch_rmsd(modified_geom_3.positions, modified_geom_1.positions) > .01

    assert not modified_geom_1.converged
    assert modified_geom_1.energy > modified_geom_3.energy
    assert modified_geom_1.gnorm > modified_geom_3.gnorm


@pytest.fixture
def orca_driver(tmpdir):
    ORCA_DRIVER = OrcaDriver(tmpdir, shutil.which('orca'), 'b3lyp', 'sto-3g')
    ORCA_DRIVER.solvatation_model = 'cpcm'
    ORCA_DRIVER.solvent = 'tetrahydrofuran'

    return ORCA_DRIVER


@pytest.mark.skipif(not shutil.which('orca'), reason='Orca driver not available')
def test_orca_get_energy_cpcm(Ca_THF2_ensemble, orca_driver):
    output_file = orca_driver.workdir / 'output.log'
    orca_driver.nprocs = 2
    with output_file.open('w') as f:
        assert orca_driver.get_energy(Ca_THF2_ensemble.elements[0], True, f) == (
            pytest.approx(-1129.026, abs=1e-2),
            pytest.approx(-1129.283, abs=1e-2)
        )


@pytest.mark.skipif(not shutil.which('orca'), reason='Orca driver not available')
def test_orca_get_energy_smd(Ca_THF2_ensemble, orca_driver):
    orca_driver.solvatation_model = 'smd'
    output_file = orca_driver.workdir / 'output.log'
    orca_driver.nprocs = 2
    with output_file.open('w') as f:
        assert orca_driver.get_energy(Ca_THF2_ensemble.elements[0], True, f) == (
            pytest.approx(-1129.023, abs=1e-2),
            pytest.approx(-1129.317, abs=1e-2)
        )


def test_orca_get_energy_from_log(tmpdir):
    cwd = pathlib.Path(__file__).parent
    drv = OrcaDriver(
        tmpdir, '{} {}'.format(cwd / 'assets/run-dummy.sh', cwd / 'assets/energy.orca.log'), 'hf', 'sto-3g')

    total_energy = drv.get_energy(Molecule([], np.zeros((0, 3))))

    assert total_energy == pytest.approx(-2859.479849944195, abs=1e-2)


@pytest.mark.skipif(not shutil.which('orca'), reason='Orca driver not available')
def test_orca_opt(Ca_THF2_ensemble, orca_driver):
    old_geometry = Ca_THF2_ensemble.elements[0]
    output_file = orca_driver.workdir / 'output.log'

    # only one cycle
    with output_file.open('w') as f:
        modified_geom_1 = orca_driver.optimize_geometry(old_geometry, True, f, maxcycle=1)
    assert rmsd.kabsch_rmsd(old_geometry.positions, modified_geom_1.positions) > .01
    assert not modified_geom_1.converged

    # 2 cycles
    with output_file.open('w') as f:
        modified_geom_3 = orca_driver.optimize_geometry(old_geometry, True, f, maxcycle=2)

    # final check
    assert rmsd.kabsch_rmsd(modified_geom_3.positions, old_geometry.positions) > .01
    assert rmsd.kabsch_rmsd(modified_geom_3.positions, modified_geom_1.positions) > .01

    assert not modified_geom_1.converged
    assert modified_geom_1.energy > modified_geom_3.energy
    assert modified_geom_1.gnorm > modified_geom_3.gnorm

    assert modified_geom_3.charge == 2
