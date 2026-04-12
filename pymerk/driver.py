import pathlib
import subprocess
import io
import sys
import select
from typing import TextIO

from pymerk.ensemble import Geometry


class BaseDriver:
    def __init__(self, workdir: pathlib.Path):
        self.workdir = workdir

    def get_energy(self, geometry: Geometry, output: TextIO = sys.stdout) -> float:
        raise NotImplementedError()

    def get_gibbs_free_energy(
            self, geometry: Geometry, T: float = 298.15, output: TextIO = sys.stdout) -> tuple[float, float]:
        raise NotImplementedError()

    def optimize_geometry(self, geometry: Geometry, output: TextIO = sys.stdout) -> tuple[Geometry, float]:
        raise NotImplementedError()


def _make_temp_xyz(workdir: str | pathlib.Path, geometry: Geometry, file_name: str = 'input.xyz') -> pathlib.Path:
    """Make a temporary xyz file"""

    xyz_path = pathlib.Path(workdir) / file_name
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

        input_path = pathlib.Path(self.workdir) / 'input.xtb'

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
            self, geometry: Geometry, output: TextIO = sys.stdout, optlevel: int = 0) -> tuple[Geometry, float]:
        xyz_path = _make_temp_xyz(self.workdir, geometry)
        command_line = self._make_command_line(geometry)

        input_path = pathlib.Path(self.workdir) / 'input.xtb'

        with input_path.open('w') as f:
            f.write('$thermo\n  optlevel={}\n$end'.format(optlevel))

        returncode, stdout, stderr = _run_and_capture(
            [self.exe_path, xyz_path, *command_line, '--opt', '-I', str(input_path)], self.workdir, output)

        if returncode != 0:
            raise RuntimeError('error while running xtb: {}'.format(stderr))

        position = stdout.rfind('TOTAL ENERGY')

        if position < 0:
            raise RuntimeError('error while running xtb: unable to find TOTAL ENERGY in output')

        total_energy = float(stdout[position + 26: position + 43])

        with (pathlib.Path(self.workdir) / 'xtbopt.xyz').open() as f:
            new_geometry = Geometry.from_xyz(f, geometry.charge, geometry.multiplicity)

        return new_geometry, total_energy
