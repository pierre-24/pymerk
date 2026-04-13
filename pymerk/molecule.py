import numpy
from typing import TextIO
from numpy.typing import NDArray


class Molecule:
    def __init__(
            self, symbols: list[str], positions: NDArray, charge: int = 0, multiplicity: int = 1, energy: float = .0):
        assert positions.shape == (len(symbols), 3)

        self.symbols = symbols
        self.positions = positions
        self.charge = charge
        self.multiplicity = multiplicity
        self.energy = energy

    def __len__(self) -> int:
        return len(self.symbols)

    def copy(self) -> 'Molecule':
        """Copy itself. Involves a copy of positions and symbols.
        """

        return Molecule(
            self.symbols.copy(),
            self.positions.copy()
        )

    @classmethod
    def from_xyz(cls, f: TextIO, charge: int = 0, multiplicity: int = 1, energy: float = .0) -> 'Molecule':
        """Read geometry from a XYZ file
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

        return cls(symbols, numpy.array(positions), charge, multiplicity, energy)

    @staticmethod
    def from_multi_xyz(f: TextIO, charge: int = 0, multiplicity: int = 1) -> list['Molecule']:
        geometries = []
        while True:
            try:
                geometries.append(Molecule.from_xyz(f, charge, multiplicity))
            except EOFError:
                break

        return geometries

    def to_string(self) -> str:
        """Get the positions"""

        r = ''
        for i in range(len(self)):
            if i > 0:
                r += '\n'
            r += '{:2} {: .7f} {: .7f} {: .7f}'.format(self.symbols[i], *self.positions[i])

        return r

    def to_xyz(self, title: str = '') -> str:
        """Get XYZ representation of this geometry"""

        r = '{}\n{}'.format(len(self), title)
        for i in range(len(self)):
            r += '\n{:2} {: .7f} {: .7f} {: .7f}'.format(self.symbols[i], *self.positions[i])

        return r
