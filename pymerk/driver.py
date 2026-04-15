import pathlib
import shutil
import subprocess
import io
import sys
import h5py
import select
from typing import TextIO, Optional

from pymerk.molecule import Molecule


class BaseDriver:
    def __init__(self, workdir: pathlib.Path):
        self.workdir = workdir

    @property
    def workdir(self) -> pathlib.Path:
        return self._workdir

    @workdir.setter
    def workdir(self, workdir: str | pathlib.Path):
        self._workdir = pathlib.Path(workdir)

    def clear_workdir(self):
        for i in self.workdir.iterdir():
            if i.is_file():
                i.unlink()
            else:
                shutil.rmtree(str(i))

    def get_energy(
            self, geometry: Molecule, add_solvent: bool = False, output: TextIO = sys.stdout
    ) -> float | tuple[float, float]:
        """Get the electronic energy of the given geometry"""
        raise NotImplementedError()

    def get_gibbs_free_energy(
            self, geometry: Molecule, T: float = 298.15, add_solvent: bool = False, output: TextIO = sys.stdout
    ) -> tuple[float, float] | tuple[float, float, float]:
        """Get the electronic energy and the Gibbs free energy (at `T`) of the given geometry"""
        raise NotImplementedError()

    def optimize_geometry(
            self, geometry: Molecule, add_solvent: bool = False, output: TextIO = sys.stdout, maxcycle: int = -1
    ) -> Molecule:
        """Optimize the given geometry, and get the optimized geometry (and its electronic energy)"""
        raise NotImplementedError()

    def __str__(self):
        return 'BaseDriver'


def _make_temp_xyz(workdir: pathlib.Path, geometry: Molecule, file_name: str = 'input.xyz') -> pathlib.Path:
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
                output.flush()
            elif stream.fileno() == stderr_fileno:
                line = process.stderr.readline()  # type: ignore
                stderrbuf.write(line)
                err_output.write(line)
                err_output.flush()
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


def _find_float(s: str, out: str, pstart: int, pend: int, label: str = 'prog') -> float:
    position = out.rfind(s)

    if position < 0:
        raise RuntimeError('error while running {}: unable to find `{}` in output'.format(label, s))

    return float(out[position + pstart: position + pend])


class QMDriver(BaseDriver):
    def __init__(self, workdir: pathlib.Path, method: str, basis: str):
        super().__init__(workdir)

        self.method = method
        self.basis = basis


class XtbDriver(BaseDriver):
    def __init__(
        self, workdir: pathlib.Path, exe_path: str | pathlib.Path, version: str = 'gfn2',
        use_bhess: bool = True, imagthr: float = -100, sthr: float = 50, scale: float = 1.0,
        optlevel: int = 0, solvatation_model: str = None, solvent: str = None
    ):
        super().__init__(workdir)

        self.exe_path = exe_path
        self.version = version

        self.solvatation_model = solvatation_model
        self.solvent = solvent

        self.use_bhess = use_bhess
        self.imagthr = imagthr
        self.sthr = sthr
        self.scale = scale

        self.optlevel = optlevel

    def __str__(self):
        return 'XtbDriver[{}{}]'.format(
            self.version,
            '' if self.solvatation_model is None else ',{}({})'.format(self.solvatation_model, self.solvent)
        )

    def _write_input(self, geometry: Molecule, add_solvent: bool, f: TextIO):
        f.write('$chrg {}\n'.format(geometry.charge))

        if geometry.multiplicity > 1:
            f.write('$spin {}\n'.format(geometry.multiplicity - 1))

        if self.version == 'gfn1':
            f.write('$gfn\n  method=1\n')
        elif self.version == 'gfn2':
            f.write('$gfn\n  method=2\n')
        elif self.version == 'gfnff':
            pass
        else:
            raise RuntimeError('unrecognized version: {}'.format(self.version))

    def _make_command_line(self, geometry: Molecule, add_solvent: bool) -> list[str]:
        command_line = []

        if add_solvent:
            if self.solvatation_model is None:
                raise RuntimeError('solvatation model is not set')

            command_line.extend(['--' + self.solvatation_model, self.solvent])

        if self.version == 'gfnff':
            command_line.extend(['--gfnff'])

        return command_line

    def get_energy(
            self, geometry: Molecule, add_solvent: bool = False, output: TextIO = sys.stdout
    ) -> float | tuple[float, float]:
        xyz_path = _make_temp_xyz(self.workdir, geometry)
        command_line = self._make_command_line(geometry, add_solvent)

        input_path = self.workdir / 'input.xtb'

        with input_path.open('w') as f:
            self._write_input(geometry, add_solvent, f)
            f.write('$end')

        returncode, stdout, stderr = _run_and_capture(
            [self.exe_path, xyz_path, '-I', str(input_path), *command_line], self.workdir, output, output)

        if returncode != 0:
            raise RuntimeError('error while running xtb: {}'.format(stderr))

        self.clear_workdir()

        total_energy = _find_float('TOTAL ENERGY', stdout, 26, 43, 'xtb')

        if add_solvent:
            gsolv = _find_float('-> Gsolv', stdout, 9, 42, 'xtb')

            return total_energy - gsolv, total_energy
        else:
            return total_energy

    def get_gibbs_free_energy(
            self, geometry: Molecule, T: float = 298.15, add_solvent: bool = False, output: TextIO = sys.stdout
    ) -> tuple[float, float] | tuple[float, float, float]:
        xyz_path = _make_temp_xyz(self.workdir, geometry)
        command_line = self._make_command_line(geometry, add_solvent)

        input_path = self.workdir / 'input.xtb'

        with input_path.open('w') as f:
            f.write('$thermo\n  temp={}'.format(T))
            if self.use_bhess:
                f.write('  imagthr={}\n  scale={}\n  sthr={}\n'.format(self.imagthr, self.scale, self.sthr))

            self._write_input(geometry, add_solvent, f)

            f.write('$end')

        returncode, stdout, stderr = _run_and_capture(
            [
                self.exe_path, xyz_path,
                *command_line,
                '--bhess' if self.use_bhess else '--ohess',
                '-I', str(input_path)
            ], self.workdir, output, output)

        if returncode != 0:
            raise RuntimeError('error while running xtb: {}'.format(stderr))

        self.clear_workdir()

        total_energy = _find_float('TOTAL ENERGY', stdout, 26, 43, 'xtb')
        total_free_energy = _find_float('TOTAL FREE ENERGY', stdout, 26, 43, 'xtb')

        if add_solvent:
            gsolv = _find_float('-> Gsolv', stdout, 9, 42, 'xtb')

            return total_energy - gsolv, total_energy, total_free_energy
        else:
            return total_energy, total_free_energy

    def optimize_geometry(
            self, geometry: Molecule, add_solvent: bool = False, output: TextIO = sys.stdout, maxcycle: int = -1
    ) -> Molecule:
        xyz_path = _make_temp_xyz(self.workdir, geometry)
        command_line = self._make_command_line(geometry, add_solvent)

        input_path = self.workdir / 'input.xtb'

        with input_path.open('w') as f:
            f.write('$opt\n  optlevel={}\n'.format(self.optlevel))
            if maxcycle > 0:
                f.write('  maxcycle={}\n'.format(maxcycle))

            self._write_input(geometry, add_solvent, f)

            f.write('$end')

        returncode, stdout, stderr = _run_and_capture(
            [self.exe_path, xyz_path, *command_line, '--opt', '-I', str(input_path)], self.workdir, output, output)

        if returncode != 0:
            raise RuntimeError('error while running xtb: {}'.format(stderr))

        total_energy = _find_float('TOTAL ENERGY', stdout, 26, 43, 'xtb')
        gnorm = _find_float('GRADIENT NORM', stdout, 26, 43, 'xtb')

        with (self.workdir / 'xtbopt.xyz').open() as f:
            new_geometry = Molecule.from_xyz(f, geometry.charge, geometry.multiplicity, total_energy, gnorm)

        position = stdout.rfind('GEOMETRY OPTIMIZATION CONVERGED')
        if position > 0:
            new_geometry.converged = True

        self.clear_workdir()
        return new_geometry


BORH_TO_ANG = 5.29177210544e-1


class VlxDriver(QMDriver):
    def __init__(
        self, workdir: pathlib.Path, exe_path: str | pathlib.Path, method: str, basis: str,
        solvatation_model: Optional[str] = None, solvent: Optional[str | float] = None,
        conv_energy: float = 1e-6, conv_grms: float = 3e-4, conv_gmax: float = 4.5e-4, conv_drms: float = 1.2e-3,
        conv_dmax: float = 1.8e-3
    ):
        super().__init__(workdir, method, basis)

        self.exe_path = exe_path
        self.solvatation_model = solvatation_model
        self.solvent = solvent

        self.conv_energy = conv_energy
        self.conv_grms = conv_grms
        self.conv_gmax = conv_gmax
        self.conv_drms = conv_drms
        self.conv_dmax = conv_dmax

    def __str__(self):
        return 'VlxDriver[{}/{}{}]'.format(
            self.method, self.basis,
            '' if self.solvatation_model is None else ',{}({})'.format(self.solvatation_model, self.solvent)
        )

    def _write_input(self, geometry: Molecule, add_solvent: bool, f: TextIO):
        f.write('@method settings\nxcfun: {}\nbasis: {}\n'.format(self.method, self.basis))

        if add_solvent:
            if self.solvatation_model is None:
                raise RuntimeError('solvatation model is not set')

            f.write('solvation model: {}\n'.format(self.solvatation_model))
            if self.solvatation_model == 'cpcm':
                f.write('cpcm epsilon: {}\n'.format(self.solvent))
            elif self.solvatation_model == 'smd':
                f.write('smd solvent: {}\n'.format(self.solvent))
            else:
                raise RuntimeError('unknown solvation model for vlx `{}`'.format(self.solvatation_model))

        f.write('@end\n')

        f.write('@molecule\ncharge: {}\nmultiplicity: {}\nxyz:\n{}\n@end\n'.format(
            geometry.charge, geometry.multiplicity, geometry.to_string()))

    def get_energy(
            self, geometry: Molecule, add_solvent: bool = False, output: TextIO = sys.stdout
    ) -> float | tuple[float, float]:
        input_path = self.workdir / 'input.vlx'

        with input_path.open('w') as f:
            f.write('@jobs\ntask: scf\n@end\n')
            self._write_input(geometry, add_solvent, f)

        returncode, stdout, stderr = _run_and_capture(
            [self.exe_path, str(input_path)], self.workdir, output)

        if returncode != 0:
            raise RuntimeError('error while running vlx: {}'.format(stderr))

        self.clear_workdir()

        total_energy = _find_float('Total Energy', stdout, 36, 56, 'vlx')

        if add_solvent:
            if self.solvatation_model.lower() == 'cpcm':
                gsolv = _find_float('Solvation Energy', stdout, 22, 43, 'vlx')
            else:
                gsolv = _find_float('Solvation Energy', stdout, 32, 53, 'vlx')

            return total_energy - gsolv, total_energy

        else:
            return total_energy

    def optimize_geometry(
            self, geometry: Molecule, add_solvent: bool = False, output: TextIO = sys.stdout, maxcycle: int = -1
    ) -> Molecule:
        input_path = self.workdir / 'input.vlx'

        with input_path.open('w') as f:
            f.write('@jobs\ntask: optimize\n@end\n')

            f.write('@optimize\nconv_energy: {}\nconv_grms: {}\nconv_gmax: {}\nconv_drms: {}\nconv_dmax: {}\n'.format(
                self.conv_energy, self.conv_grms, self.conv_gmax, self.conv_drms, self.conv_dmax
            ))

            if maxcycle > 0:
                f.write('max_iter: {}\nconv_maxiter: true\n'.format(maxcycle))

            f.write('@end\n')

            self._write_input(geometry, add_solvent, f)

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

        # attempt to check whether the geometry was optimized
        position = stdout.rfind('* Info *   Gradient')
        position_end = stdout.find('a.u.', position)

        grms = float(stdout[position + 22: position_end])

        position = stdout.find('* Info * ', position_end)
        position_end = stdout.find('a.u.', position)
        gmax = float(stdout[position + 22: position_end])

        position = stdout.find(' Statistical Deviation between', position_end)

        last_opt_line = stdout[position - 47 - 123 * 3: position - 47 - 123 * 2].split()
        de = float(last_opt_line[2])
        drms = float(last_opt_line[3])
        dmax = float(last_opt_line[4])

        is_converged = (
            de < self.conv_energy
            and drms < self.conv_drms
            and dmax < self.conv_dmax
            and grms < self.conv_grms
            and gmax < self.conv_gmax
        )

        # extract last geometry
        with h5py.File(self.workdir / 'input.h5') as f:
            new_position = f['atom_coordinates'][:] * BORH_TO_ANG
            new_geometry = Molecule(
                geometry.symbols, new_position,
                geometry.charge, geometry.multiplicity,
                total_energy, grms, is_converged
            )

        self.clear_workdir()
        return new_geometry
