from typing import TextIO

from pymerk.geometry import Geometry


class Ensemble:
    """Ensemble of geometries, with associated energies"""

    def __init__(self, elements: list[tuple[Geometry, float]]):

        for i in range(1, len(elements)):
            if len(elements[i][0]) != len(elements[i - 1][0]):
                raise ValueError('Ensemble geometries must have same size')

        self.elements: list[tuple[[Geometry, float]]] = elements

    def __len__(self) -> int:
        return len(self.elements)

    @classmethod
    def from_multi_xyz(cls, f: TextIO, charge: int = 0, multiplicity: int = 1) -> 'Ensemble':
        geometries = Geometry.from_multi_xyz(f, charge, multiplicity)
        return cls([(g, .0) for g in geometries])

    def as_multi_xyz(self, f: TextIO):
        for i, (geometry, energy) in enumerate(self.elements):
            if i > 0:
                f.write('\n')
            f.write(geometry.to_xyz(title='Conformer #{}, E={}'.format(i + 1, energy)))

    def filter(self, predicate) -> 'Ensemble':
        """Filter elements, keep them if `predicate` is `True`
        """
        return Ensemble(list(filter(lambda x: predicate(x), self.elements)))
