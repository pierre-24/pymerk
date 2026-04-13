import pathlib
import subprocess
import io
import sys
import h5py
import select
from typing import TextIO

from pymerk.ensemble import Geometry


class BaseDriver:
    def __init__(self, workdir: pathlib.Path):
        self.workdir = workdir

    @property
    def workdir(self) -> pathlib.Path:
        return self._workdir

    @workdir.setter
    def workdir(self, workdir: str | pathlib.Path):
        self._workdir = pathlib.Path(workdir)

    def get_energy(self, geometry: Geometry, output: TextIO = sys.stdout) -> float:
        """Get the electronic energy of the given geometry"""
        raise NotImplementedError()

    def get_gibbs_free_energy(
            self, geometry: Geometry, T: float = 298.15, output: TextIO = sys.stdout) -> tuple[float, float]:
        """Get the electronic energy and the Gibbs free energy (at `T`) of the given geometry"""
        raise NotImplementedError()

    def optimize_geometry(
            self, geometry: Geometry, output: TextIO = sys.stdout, maxcyle: int = -1) -> tuple[Geometry, float]:
        """Optimize the given geometry, and get the optimized geometry (and its electronic energy)"""
        raise NotImplementedError()


def _make_temp_xyz(workdir: pathlib.Path, geometry: Geometry, file_name: str = 'input.xyz') -> pathlib.Path:
    """Make a temporary xyz file"""

    xyz_path = workdir / file_name
    with xyz_path.open('w') as f:
        f.write(geometry.to_xyz())

    return xyz_path


def _run_and_capture(
    cmd: list[str], cwd: str | pathlib.Path, output: TextIO = sys.stdout, err_output: TextIO = sys.stderr
) -> tuple[int, str, str]:
    """Run `cmd` and capture stdout and stderr, while also writing to `output`

    From https://me.micahrl.com/blog/magicrun/
    """

    process = subprocess.Popen(  # type: ignore
        cmd,
        bufsize=1,  # Output is line buffered, required to print output in real time
        universal_newlines=True,  # Required for line buffering
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )

    stdoutbuf = io.StringIO()
    stderrbuf = io.StringIO()
    stdout_fileno = process.stdout.fileno()  # type: ignore
    stderr_fileno = process.stderr.fileno()  # type: ignore

    while process.poll() is None:
        # select() waits until there is data to read (or an "exceptional case") on any of the streams
        readready, writeready, exceptionready = select.select(
            [process.stdout, process.stderr],
            [],
            [process.stdout, process.stderr],
            0.5,
        )

        # Check if what is ready is a stream, and if so, which stream.
        # Copy the stream to the buffer so we can use it, and print it to stdout/stderr in real time
        for stream in readready:
            if stream.fileno() == stdout_fileno:
                line = process.stdout.readline()  # type: ignore
                stdoutbuf.write(line)
                output.write(line)
            elif stream.fileno() == stderr_fileno:
                line = process.stderr.readline()  # type: ignore
                stderrbuf.write(line)
                err_output.write(line)
            else:
                raise RuntimeError(f'Unknown file descriptor in select result. Fileno: {stream.fileno()}')

    # Check for any remaining output after the process has exited.
    # Without this, the last line of output may not be printed, if output is buffered (very normal)
    # and the process doesn't explicitly flush upon exit
    # (also very normal, and will definitely happen if the process crashes or gets KILLed).
    for stream in [process.stdout, process.stderr]:
        for line in stream.readlines():
            if stream.fileno() == stdout_fileno:
                stdoutbuf.write(line)
                output.write(line)
            elif stream.fileno() == stderr_fileno:
                stderrbuf.write(line)
                err_output.write(line)

    return process.wait(), stdoutbuf.getvalue(), stderrbuf.getvalue()


class QMDriver(BaseDriver):
    def __init__(self, workdir: pathlib.Path, method: str, basis: str):
        super().__init__(workdir)

        self.method = method
        self.basis = basis


class XtbDriver(BaseDriver):
    def __init__(self, workdir: pathlib.Path, exe_path: str | pathlib.Path, version: str = 'gfn2'):
        super().__init__(workdir)

        self.exe_path = exe_path
        self.version = version
        self.solvatation_model = None
        self.solvent = None

    def _make_command_line(self, geometry: Geometry) -> list[str]:
        command_line = []
        if geometry.charge != 0:
            command_line.extend(['-c', str(geometry.charge)])

        if self.solvatation_model is not None:
            command_line.extend(['--' + self.solvatation_model, self.solvent])

        if self.version == 'gfn1':
            command_line.extend(['--gfn', '1'])
        elif self.version == 'gfnff':
            command_line.extend(['--gfnff'])
        elif self.version == 'gfn2':
            pass
        else:
            raise RuntimeError('unrecognized version: {}'.format(self.version))

        if geometry.multiplicity > 1:
            command_line.extend(['--uhf', str(geometry.multiplicity - 1)])

        return command_line

    def get_energy(self, geometry: Geometry, output: TextIO = sys.stdout) -> float:
        xyz_path = _make_temp_xyz(self.workdir, geometry)
        command_line = self._make_command_line(geometry)

        returncode, stdout, stderr = _run_and_capture(
            [self.exe_path, xyz_path, *command_line], self.workdir, output)

        if returncode != 0:
            raise RuntimeError('error while running xtb: {}'.format(stderr))

        position = stdout.rfind('TOTAL ENERGY')

        if position < 0:
            raise RuntimeError('error while running xtb: unable to find TOTAL ENERGY in output')

        return float(stdout[position + 26: position + 43])

    def get_gibbs_free_energy(
        self, geometry: Geometry, T: float = 298.15, output: TextIO = sys.stdout,
        use_bhess: bool = True, imagthr: float = -100, sthr: float = 50, scale: float = 1.0
    ) -> tuple[float, float]:
        xyz_path = _make_temp_xyz(self.workdir, geometry)
        command_line = self._make_command_line(geometry)

        input_path = self.workdir / 'input.xtb'

        with input_path.open('w') as f:
            f.write('$thermo\n  temp={}'.format(T))

            if use_bhess:
                f.write('  imagthr={}\n  scale={}\n  sthr={}\n'.format(imagthr, scale, sthr))

            f.write('$end')

        returncode, stdout, stderr = _run_and_capture(
            [
                self.exe_path, xyz_path,
                *command_line,
                '--bhess' if use_bhess else '--ohess',
                '-I', str(input_path)
            ], self.workdir, output)

        if returncode != 0:
            raise RuntimeError('error while running xtb: {}'.format(stderr))

        position = stdout.rfind('TOTAL ENERGY')

        if position < 0:
            raise RuntimeError('error while running xtb: unable to find TOTAL ENERGY in output')

        total_energy = float(stdout[position + 26: position + 43])
        position = stdout.find('TOTAL FREE ENERGY', position)

        if position < 0:
            raise RuntimeError('error while running xtb: unable to find TOTAL FREE ENERGY in output')

        total_free_energy = float(stdout[position + 26: position + 43])

        return total_energy, total_free_energy

    def optimize_geometry(
            self, geometry: Geometry, output: TextIO = sys.stdout, maxcycle: int = -1, optlevel: int = 0
    ) -> tuple[Geometry, float]:
        xyz_path = _make_temp_xyz(self.workdir, geometry)
        command_line = self._make_command_line(geometry)

        input_path = self.workdir / 'input.xtb'

        with input_path.open('w') as f:
            f.write('$opt\n  optlevel={}\n'.format(optlevel))

            if maxcycle > 0:
                f.write('  maxcycle={}\n'.format(maxcycle))

            f.write('$end')

        returncode, stdout, stderr = _run_and_capture(
            [self.exe_path, xyz_path, *command_line, '--opt', '-I', str(input_path)], self.workdir, output)

        if returncode != 0:
            raise RuntimeError('error while running xtb: {}'.format(stderr))

        position = stdout.rfind('TOTAL ENERGY')

        if position < 0:
            raise RuntimeError('error while running xtb: unable to find TOTAL ENERGY in output')

        total_energy = float(stdout[position + 26: position + 43])

        with (self.workdir / 'xtbopt.xyz').open() as f:
            new_geometry = Geometry.from_xyz(f, geometry.charge, geometry.multiplicity)

        return new_geometry, total_energy


BORH_TO_ANG = 5.29177210544e-1


class VlxDriver(QMDriver):
    def __init__(self, workdir: pathlib.Path, exe_path: str | pathlib.Path, method: str, basis: str):
        super().__init__(workdir, method, basis)

        self.exe_path = exe_path
        self.solvatation_model = None
        self.solvent = None

    def _write_input(self, f: TextIO, geometry: Geometry):
        f.write('@method settings\nxcfun: {}\nbasis: {}\n@end\n'.format(self.method, self.basis))

        f.write('@molecule\ncharge: {}\nmultiplicity: {}\nxyz:\n{}\n@end\n'.format(
            geometry.charge, geometry.multiplicity, geometry.to_string()))

    def get_energy(self, geometry: Geometry, output: TextIO = sys.stdout) -> float:
        input_path = self.workdir / 'input.vlx'

        with input_path.open('w') as f:
            f.write('@jobs\ntask: scf\n@end\n')
            self._write_input(f, geometry)

        returncode, stdout, stderr = _run_and_capture(
            [self.exe_path, str(input_path)], self.workdir, output)

        if returncode != 0:
            raise RuntimeError('error while running vlx: {}'.format(stderr))

        position = stdout.rfind('Total Energy')

        if position < 0:
            raise RuntimeError('error while running vlx: unable to find `Total energy` in output')

        total_energy = float(stdout[position + 36: position + 56])

        return total_energy

    def optimize_geometry(
        self, geometry: Geometry, output: TextIO = sys.stdout, maxcycle: int = -1,
        conv_energy: float = 1e-6, conv_grms: float = 3e-4, conv_gmax: float = 4.5e-4, conv_drms: float = 1.2e-3,
        conv_dmax: float = 1.8e-3
    ) -> tuple[Geometry, float]:
        input_path = self.workdir / 'input.vlx'

        with input_path.open('w') as f:
            f.write('@jobs\ntask: optimize\n@end\n')

            f.write('@optimize\nconv_energy: {}\nconv_grms: {}\nconv_gmax: {}\nconv_drms: {}\nconv_dmax: {}\n'.format(
                conv_energy, conv_grms, conv_gmax, conv_drms, conv_dmax
            ))

            if maxcycle > 0:
                f.write('max_iter: {}\nconv_maxiter: true\n'.format(maxcycle))

            f.write('@end\n')

            self._write_input(f, geometry)

        returncode, stdout, stderr = _run_and_capture(
            [self.exe_path, str(input_path)], self.workdir, output)

        if returncode != 0:
            raise RuntimeError('error while running vlx: {}'.format(stderr))

        position = stdout.rfind('Geometry optimization completed.')

        if position < 0:
            raise RuntimeError('error while running vlx: unable to find `Geometry optimization completed` in output')

        position = stdout.rfind('* Info *   Energy')
        position_end = stdout.find('a.u.', position)

        total_energy = float(stdout[position + 22: position_end])

        with h5py.File(self.workdir / 'input.h5') as f:
            new_position = f['atom_coordinates'][:] * BORH_TO_ANG
            new_geometry = Geometry(geometry.symbols, new_position, geometry.charge, geometry.multiplicity)

        return new_geometry, total_energy
