from pymerk.ensemble import Ensemble
from pymerk.filter import BaseFilter


def test_filter():
    BaseFilter(Ensemble([]))
