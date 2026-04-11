from pymerk.ensemble import Ensemble


class BaseFilter:
    """Base class for filters"""

    def __init__(self, ensemble: Ensemble):
        self.ensemble = ensemble

    def filter(self) -> Ensemble:
        raise NotImplementedError()
