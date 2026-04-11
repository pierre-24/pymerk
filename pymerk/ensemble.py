from typing import TextIO

from pymerk.geometry import Geometry


class Ensemble:
    """Ensemble of geometries, with associated energies"""

    def __init__(self, elements: list[tuple[Geometry, float]]):

        for i in range(1, len(elements)):
            if len(elements[i][0]) != len(elements[i - 1][0]):
                raise ValueError('Ensemble geometries must have same size')

        self.elements = elements

    def __len__(self) -> int:
        return len(self.elements)

    @classmethod
    def from_multi_xyz(cls, f: TextIO) -> 'Ensemble':
        geometries = Geometry.from_multi_xyz(f)
        return cls([(g, .0) for g in geometries])
