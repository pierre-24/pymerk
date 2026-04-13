from typing import TextIO

from pymerk.molecule import Molecule


class Ensemble:
    """Ensemble of geometries, with associated energies"""

    def __init__(self, elements: list[Molecule]):

        for i in range(1, len(elements)):
            if len(elements[i]) != len(elements[i - 1]):
                raise ValueError('Ensemble geometries must have same size')

        self.elements: list[Molecule] = elements

    def __len__(self) -> int:
        return len(self.elements)

    @classmethod
    def from_multi_xyz(cls, f: TextIO, charge: int = 0, multiplicity: int = 1) -> 'Ensemble':
        geometries = Molecule.from_multi_xyz(f, charge, multiplicity)
        return cls(geometries)

    def as_multi_xyz(self, f: TextIO):
        for i, geometry in enumerate(self.elements):
            if i > 0:
                f.write('\n')
            f.write(geometry.to_xyz(title='Conformer #{}, E={}'.format(i + 1, geometry.energy)))

    def filter(self, predicate) -> 'Ensemble':
        """Filter elements, keep them if `predicate` is `True`
        """
        return Ensemble(list(filter(lambda x: predicate(x), self.elements)))
