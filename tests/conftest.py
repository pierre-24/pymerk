import pytest
import pathlib

from pymerk.ensemble import Ensemble


@pytest.fixture(scope='module')
def Ca_THF2_ensemble():
    with (pathlib.Path(__file__).parent / 'assets/Ca_THF2_ensemble.xyz').open() as f:
        return Ensemble.from_multi_xyz(f, charge=2)
