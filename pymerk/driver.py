import pathlib
import tempfile
import subprocess

from pymerk.ensemble import Geometry


class BaseDriver:
    def get_energy(self, geometry: Geometry) -> float:
        raise NotImplementedError()

    def get_gibbs_free_energy(self, geometry: Geometry) -> tuple[float, float]:
        raise NotImplementedError()

    def optimize(self, geometry: Geometry) -> Geometry:
        raise NotImplementedError()


class QMDriver(BaseDriver):
    def __init__(self, method: str, basis: str):
        self.method = method
        self.basis = basis


class XtbDriver(BaseDriver):
    def __init__(self, path: str, version: str = 'gfn2'):
        self.path = path
        self.version = version
        self.solvatation_model = None
        self.solvent = None

    def get_energy(self, geometry: Geometry) -> float:
        tdir = tempfile.mkdtemp()
        xyz_path = pathlib.Path(tdir) / 'input.xyz'
        with xyz_path.open('w') as f:
            f.write(geometry.to_xyz())

        command_line = []
        if geometry.charge != 0:
            command_line.extend(['-c', str(geometry.charge)])

        if self.solvatation_model is not None:
            command_line.extend(['-' + self.solvatation_model, self.solvent])

        result = subprocess.run(
            [self.path, xyz_path, *command_line], stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        if result.returncode != 0:
            raise RuntimeError('error while running xtb: {}'.format(result.stderr.decode('utf-8')))

        stdout = result.stdout.decode('utf-8')
        position = stdout.rfind('TOTAL ENERGY')

        if position < 0:
            raise RuntimeError('error while running xtb: unable to find TOTAL ENERGY in output')

        return float(stdout[position + 26: position + 43])
