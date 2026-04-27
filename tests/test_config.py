import pathlib
import pytest

from pymerk.scripts import Config


@pytest.fixture
def config_x():
    with (pathlib.Path(__file__).parent / 'assets/input_x.toml').open('rb') as f:
        return Config.from_toml(f)


def test_config(config_x):
    assert config_x.paths.vlx == 'vlx'
    assert config_x.general.solvent == 'thf'
