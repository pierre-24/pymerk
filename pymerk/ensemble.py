from typing import TextIO

from pymerk.molecule import Molecule


class Ensemble:
    """Container for a collection of `Molecule` instances with associated properties.
    
    An ensemble represents multiple conformations or isomers of the same molecular
    structure. All geometries must have the same number of atoms. This class provides
    operations for filtering, I/O, and bulk manipulation of conformers.
    
    Attributes:
        elements: List of `Molecule` instances representing each geometry.
    """

    def __init__(self, elements: list[Molecule]):
        """Initialize an Ensemble.
        
        Args:
            elements: List of `Molecule` instances to include in the ensemble.
            
        Raises:
            ValueError: If molecules have different numbers of atoms.
        """

        for i in range(1, len(elements)):
            if len(elements[i]) != len(elements[i - 1]):
                raise ValueError('Ensemble geometries must have same size')

        self.elements: list[Molecule] = elements

    def __len__(self) -> int:
        """Return the number of geometries in the ensemble.
        
        Returns:
            Number of molecules in `elements`.
        """
        return len(self.elements)

    @classmethod
    def from_multi_xyz(cls, f: TextIO, charge: int = 0, multiplicity: int = 1) -> 'Ensemble':
        """Create an ensemble by reading multiple geometries from an XYZ file.
        
        Reads consecutive XYZ blocks from the file and creates one `Molecule` per block.
        Each geometry is automatically named 'Conformer #N' where `N` is the 1-based index.
        
        Args:
            f: File object opened in read mode containing XYZ format data.
            charge: Total molecular charge for all geometries. Defaults to 0.
            multiplicity: Spin multiplicity for all geometries. Defaults to 1.
            
        Returns:
            A new `Ensemble` instance containing all parsed geometries.
        """
        geometries = Molecule.from_multi_xyz(f, charge, multiplicity, names=lambda i: 'Conformer #{}'.format(i + 1))
        return cls(geometries)

    def as_multi_xyz(self, f: TextIO):
        """Write the ensemble to a file in multi-XYZ format.
        
        Writes each geometry as an XYZ block with title containing the name and energy.
        Consecutive blocks are separated by a blank line.
        
        Args:
            f: File object opened in write mode.
        """
        for i, geometry in enumerate(self.elements):
            if i > 0:
                f.write('\n')
            f.write(geometry.to_xyz(title='{}, E={}'.format(geometry.name, geometry.energy)))

    def filter(self, predicate) -> 'Ensemble':
        """Filter geometries by a predicate function.
        
        Creates a new ensemble containing only molecules for which `predicate` returns `True`.
        
        Args:
            predicate: Callable that takes a `Molecule` and returns `bool`.
            
        Returns:
            A new `Ensemble` with filtered molecules. Returns empty ensemble if no molecules match.
        """
        return Ensemble(list(filter(lambda x: predicate(x), self.elements)))
