import pytest
import pathlib

from pymerk.ensemble import Ensemble


@pytest.fixture(scope='module')
def Ca_THF3_ensemble():
    with (pathlib.Path(__file__).parent / 'assets/Ca_THF3_ensemble.xyz').open() as f:
        return Ensemble.from_multi_xyz(f)
