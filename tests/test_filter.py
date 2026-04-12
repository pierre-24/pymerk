from pymerk.driver import BaseDriver
from pymerk.filter import BaseFilter


def test_filter(Ca_THF3_ensemble, tmpdir):
    BaseFilter(Ca_THF3_ensemble, BaseDriver(tmpdir))
