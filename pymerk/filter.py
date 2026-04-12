from pymerk.ensemble import Ensemble
from pymerk.driver import BaseDriver


class BaseFilter:
    """Base class for filters"""

    def __init__(self, ensemble: Ensemble, driver: BaseDriver):
        self.ensemble = ensemble
        self.driver = driver

    def filter(self) -> Ensemble:
        raise NotImplementedError()
