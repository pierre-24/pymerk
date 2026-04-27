import pytest
import pathlib
import shutil

from pymerk.scripts import Config
from pymerk.scripts.run import DefaultWorkflow


@pytest.fixture
def config_2H2O():
    with (pathlib.Path(__file__).parent / 'assets/input_2H2O.toml').open('rb') as f:
        return Config.from_toml(f)


@pytest.mark.skipif(not shutil.which('vlx'), reason='vlx driver not available')
def test_default_workflow_vlx(config_2H2O, _2H2O_ensemble, tmpdir):
    # Execute Workflow
    workflow = DefaultWorkflow(tmpdir, config_2H2O)
    final_ensemble = workflow.filter(_2H2O_ensemble)

    assert len(final_ensemble) == 1
