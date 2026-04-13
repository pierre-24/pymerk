from pymerk.driver import BaseDriver
from pymerk.filter import BaseFilter


def test_filter(Ca_THF2_ensemble, tmpdir):
    BaseFilter(Ca_THF2_ensemble, BaseDriver(tmpdir))
