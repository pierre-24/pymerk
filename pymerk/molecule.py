import numpy
from typing import TextIO, Callable, Optional
from numpy.typing import NDArray


class Molecule:
    """Represents a molecular geometry with atomic positions and associated properties.
    
    This class stores atomic symbols, their 3D positions, and quantum chemical properties
    such as energy, convergence status, and gradient norm.
    
    Attributes:
        symbols: List of atomic symbols (e.g., `['C', 'H', 'O']`).
        positions: `(N, 3)` array of atomic coordinates in Angstroms.
        charge: Total molecular charge in elementary charge units. Defaults to 0.
        multiplicity: Spin multiplicity (2S+1). Defaults to 1 (singlet).
        energy: Electronic energy in Hartree. Defaults to 0.0.
        gnorm: Gradient norm (max. atomic force) in Hartree/Bohr. Defaults to 0.
        converged: Whether geometry optimization has converged. Defaults to False.
        name: Identifier for this geometry (e.g., 'Conformer #1'). Defaults to empty string.
    """
    
    def __init__(
            self, symbols: list[str], positions: NDArray, charge: int = 0, multiplicity: int = 1,
            energy: float = .0, gnorm: float = 0, converged: bool = False, name: str = ''
    ):
        """Initialize a `Molecule` instance.
        
        Args:
            symbols: List of atomic element symbols.
            positions: `(N, 3)` numpy array of atomic coordinates.
            charge: Total molecular charge. Defaults to 0.
            multiplicity: Spin multiplicity. Defaults to 1.
            energy: Electronic energy in Hartree. Defaults to 0.0.
            gnorm: Gradient norm in Hartree/Bohr. Defaults to 0.
            converged: Convergence status flag. Defaults to False.
            name: Geometry identifier string. Defaults to empty string.
            
        Raises:
            AssertionError: If positions shape is not (len(symbols), 3).
        """
        assert positions.shape == (len(symbols), 3)

        self.symbols = symbols
        self.positions = positions
        self.charge = charge
        self.multiplicity = multiplicity
        self.energy = energy
        self.gnorm = gnorm
        self.converged = converged
        self.name = name

    def __len__(self) -> int:
        """Return the number of atoms in the molecule.
        
        Returns:
            Number of atoms.
        """
        return len(self.symbols)

    def copy(self) -> 'Molecule':
        """Create a deep copy of this molecule.
        
        Returns a new Molecule instance with independent copies of positions and symbols
        to ensure modifications do not affect the original.
        
        Returns:
            A new Molecule instance with copied data.
        """

        return Molecule(
            self.symbols.copy(),
            self.positions.copy(),
            self.charge,
            self.multiplicity,
            self.energy,
            self.gnorm,
            self.converged,
            self.name
        )

    @classmethod
    def from_xyz(
            cls, f: TextIO, charge: int = 0, multiplicity: int = 1,
            energy: float = .0, gnorm: float = 0, converged: bool = False, name: str = ''
    ) -> 'Molecule':
        """Read a single geometry from an XYZ format file, and add a bunch of properties.
        
        Args:
            f: File object opened in read mode.
            charge: Total molecular charge. Defaults to 0.
            multiplicity: Spin multiplicity. Defaults to 1.
            energy: Electronic energy in Hartree. Defaults to 0.0.
            gnorm: Gradient norm. Defaults to 0.
            converged: Convergence flag. Defaults to False.
            name: Geometry identifier. Defaults to empty string.
            
        Returns:
            A new `Molecule` instance parsed from the file.
            
        Raises:
            EOFError: If no data is found in the file.
            RuntimeError: If XYZ format is invalid (each atom line must have 4 chunks).
        """

        symbols = []
        positions = []

        data = f.readline()
        if data == '':  # nothing else to read?
            raise EOFError

        n = int(data)
        f.readline()

        for i in range(n):
            chunks = f.readline().split()
            if len(chunks) != 4:
                raise RuntimeError('Invalid XYZ format, each line must contain 4 chunks')

            symbols.append(chunks[0])
            positions.append([float(x) for x in chunks[1:]])

        return cls(symbols, numpy.array(positions), charge, multiplicity, energy, gnorm, converged, name)

    @staticmethod
    def from_multi_xyz(
            f: TextIO, charge: int = 0, multiplicity: int = 1,
            names: Callable[[int], str] = lambda g: str(g)
    ) -> list['Molecule']:
        """Read multiple geometries from an XYZ file.
        
        Reads consecutive XYZ blocks from a file, creating one `Molecule` per block.
        Stops at EOF.
        
        Args:
            f: File object opened in read mode.
            charge: Total molecular charge. Defaults to 0.
            multiplicity: Spin multiplicity. Defaults to 1.
            names: Function mapping geometry index to name. Defaults to string representation of index.
            
        Returns:
            List of `Molecule` instances parsed from the file.
        """
        geometries = []
        i = 0
        while True:
            try:
                geometries.append(Molecule.from_xyz(f, charge, multiplicity, name=names(i)))
                i += 1
            except EOFError:
                break

        return geometries

    def to_string(self) -> str:
        """Generate XYZ atom block (symbols and coordinates) as a string.
        
        Returns:
            Formatted string with atom coordinates.
        """

        r = ''
        for i in range(len(self)):
            if i > 0:
                r += '\n'
            r += '{:2} {: .7f} {: .7f} {: .7f}'.format(self.symbols[i], *self.positions[i])

        return r

    def to_xyz(self, title: Optional[str] = None) -> str:
        """Generate complete XYZ format representation of this `Molecule`.
        
        Args:
            title: Optional title/comment line. If None, uses `self.name`.
            
        Returns:
            Complete XYZ formatted string.
        """

        r = '{}\n{}\n'.format(len(self), title if title is not None else self.name)
        r += self.to_string()

        return r
