import pathlib
import shutil
import subprocess
import io
import sys
import h5py
import select
import numpy

from typing import TextIO, Optional

from pymerk.molecule import Molecule


class BaseDriver:
    """Abstract base class for quantum chemistry program interfaces.

    Defines the interface for computing electronic energies, Gibbs free energies,
    and optimized geometries using external QM programs.

    Attributes:
        workdir: Working directory for temporary files and program output.
    """

    def __init__(self, workdir: pathlib.Path):
        """Initialize a BaseDriver.

        Args:
            workdir: Path to working directory.
        """
        self.workdir = workdir

    @property
    def workdir(self) -> pathlib.Path:
        """Get the working directory.

        Returns:
            Path object of the working directory.
        """
        return self._workdir

    @workdir.setter
    def workdir(self, workdir: str | pathlib.Path):
        """Set the working directory.

        Args:
            workdir: Path as string or `pathlib.Path` object.
        """
        self._workdir = pathlib.Path(workdir)

    def clear_workdir(self):
        """Remove all files and directories from the working directory."""
        for i in self.workdir.iterdir():
            if i.is_file():
                i.unlink()
            else:
                shutil.rmtree(str(i))

    def get_energy(
            self, geometry: Molecule, add_solvent: bool = False, output: TextIO = sys.stdout
    ) -> float | tuple[float, float]:
        """Compute the electronic energy of a geometry.

        Args:
            geometry: `Molecule` to evaluate.
            add_solvent: If `True`, compute solvated energy. Defaults to `False`.
            output: File object for writing program output. Defaults to `sys.stdout`.

        Returns:
            Electronic energy in Hartree if `add_solvent=False`.
            Tuple of `(electronic_energy, solvated_energy)` if `add_solvent=True`.

        Raises:
            NotImplementedError: Must be implemented by subclasses.
        """
        raise NotImplementedError()

    def get_gibbs_free_energy(
            self, geometry: Molecule, T: float = 298.15, add_solvent: bool = False, output: TextIO = sys.stdout
    ) -> tuple[float, float] | tuple[float, float, float]:
        """Compute electronic and Gibbs free energy of a geometry.

        Args:
            geometry: `Molecule` to evaluate.
            T: Temperature in Kelvin. Defaults to 298.15.
            add_solvent: If `True`, include solvation corrections. Defaults to `False`.
            output: File object for writing program output. Defaults to `sys.stdout`.

        Returns:
            Tuple of `(electronic_energy, gibbs_free_energy)` if `add_solvent=False`.
            Tuple of `(electronic_energy, solvated_energy, gibbs_free_energy)` if `add_solvent=True`.
            All energies in Hartree.

        Raises:
            NotImplementedError: Must be implemented by subclasses.
        """
        raise NotImplementedError()

    def optimize_geometry(
            self, geometry: Molecule, add_solvent: bool = False, output: TextIO = sys.stdout, maxcycle: int = -1,
            optlevel: str = 'normal'
    ) -> Molecule:
        """Optimize the geometry of a molecule.

        Performs geometry optimization and returns the optimized geometry with its energy
        and convergence status.

        Args:
            geometry: `Molecule` to optimize.
            add_solvent: If `True`, optimize with solvation model. Defaults to `False`.
            output: File object for writing program output. Defaults to `sys.stdout`.
            maxcycle: Maximum optimization cycles. If -1, use default. Defaults to -1.
            optlevel: Optimization level ('loose', 'normal', 'tight', etc.). Defaults to 'normal'.

        Returns:
            Optimized `Molecule` with converged flag set and gradient norm computed.

        Raises:
            NotImplementedError: Must be implemented by subclasses.
        """
        raise NotImplementedError()

    def __str__(self):
        """Return driver name.

        Returns:
            String representation of the driver type.
        """
        return 'BaseDriver'


def _make_temp_xyz(workdir: pathlib.Path, geometry: Molecule, file_name: str = 'input.xyz') -> pathlib.Path:
    """Create a temporary XYZ file from a geometry.

    Args:
        workdir: Working directory where the file will be created.
        geometry: `Molecule` to write to file.
        file_name: Name of output file. Defaults to 'input.xyz'.

    Returns:
        Path to the created XYZ file.
    """

    xyz_path = workdir / file_name
    with xyz_path.open('w') as f:
        f.write(geometry.to_xyz())

    return xyz_path


def _run_and_capture(
        cmd: list[str], cwd: str | pathlib.Path, output: TextIO = sys.stdout, err_output: TextIO = sys.stderr
) -> tuple[int, str, str]:
    """Execute a command while capturing and streaming output.

    Runs a subprocess command, capturing both stdout and stderr while simultaneously
    writing to provided output streams for real-time visibility.

    Note: Based on https://me.micahrl.com/blog/magicrun/

    Args:
        cmd: List of command and arguments to execute.
        cwd: Working directory for command execution.
        output: File object for stdout. Defaults to `sys.stdout`.
        err_output: File object for stderr. Defaults to `sys.stderr`.

    Returns:
        Tuple of `(return_code, stdout_string, stderr_string)`.
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
    """Extract a float from program output using substring search.

    Searches for a substring in output text and extracts a float value from a specific
    character range relative to the substring position.

    Args:
        s: Substring to search for in output.
        out: Program output string to search in.
        pstart: Start position offset relative to substring for float extraction.
        pend: End position offset relative to substring for float extraction.
        label: Program name for error messages. Defaults to `'prog'`.

    Returns:
        Extracted float value.

    Raises:
        RuntimeError: If substring is not found in output.
    """
    position = out.rfind(s)

    if position < 0:
        raise RuntimeError('error while running {}: unable to find `{}` in output'.format(label, s))

    return float(out[position + pstart: position + pend])


class QMDriver(BaseDriver):
    """Abstract base class for quantum chemistry drivers.

    Extends `BaseDriver` with method and basis set information common to all QM programs.

    Attributes:
        method: Exchange-correlation functional or QM method name.
        basis: Basis set specification.
    """

    def __init__(self, workdir: pathlib.Path, method: str, basis: str):
        """Initialize a QMDriver.

        Args:
            workdir: Path to working directory.
            method: QM method or functional name.
            basis: Basis set name.
        """
        super().__init__(workdir)

        self.method = method
        self.basis = basis


class XtbDriver(BaseDriver):
    """Interface for the xTB semiempirical QM program.

    Supports GFN1, GFN2, and GFNFF methods with optional solvation corrections.

    Attributes:
        exe_path: Path to xtb executable.
        version: xTB method version (`'gfn1'`, `'gfn2'`, `'gfnff'`).
        solvatation_model: Solvation model name (e.g., `'gbsa'`). `None` for gas-phase.
        solvent: Solvent name or identifier.
        use_bhess: If `True`, use Bhess for Hessian. Defaults to `True`.
        imagthr: Imaginary frequency threshold in cm⁻¹. Defaults to -100.
        sthr: Wave number threshold in cm⁻¹. Defaults to 50.
        scale: Frequency scaling factor. Defaults to 1.0.
    """

    def __init__(
            self, workdir: pathlib.Path, exe_path: str | pathlib.Path, version: str = 'gfn2',
            use_bhess: bool = True, imagthr: float = -100, sthr: float = 50, scale: float = 1.0,
            solvatation_model: str = None, solvent: str = None
    ):
        """Initialize an XtbDriver.

        Args:
            workdir: Path to working directory.
            exe_path: Path to xtb executable.
            version: xTB version to use. Defaults to 'gfn2'.
            use_bhess: Use Bhess for Hessian computation. Defaults to `True`.
            imagthr: Imaginary frequency threshold. Defaults to -100.
            sthr: Wave number threshold. Defaults to 50.
            scale: Frequency scaling factor. Defaults to 1.0.
            solvatation_model: Solvation model name. Defaults to `None`.
            solvent: Solvent identifier. Defaults to `None`.
        """
        super().__init__(workdir)

        self.exe_path = exe_path
        self.version = version

        self.solvatation_model = solvatation_model
        self.solvent = solvent

        self.use_bhess = use_bhess
        self.imagthr = imagthr
        self.sthr = sthr
        self.scale = scale

    def __str__(self):
        return 'XtbDriver[{}{}]'.format(
            self.version,
            '' if self.solvatation_model is None else ',{}({})'.format(self.solvatation_model, self.solvent)
        )

    def _write_input(self, geometry: Molecule, add_solvent: bool, f: TextIO):
        """Write xTB input file control block.

        Args:
            geometry: `Molecule` to write charge and multiplicity for.
            add_solvent: If `True`, include solvation directives (not used here).
            f: File object to write to.
        """
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
        """Build xTB command-line arguments.

        Args:
            geometry: `Molecule` being evaluated (not directly used).
            add_solvent: If `True`, add solvation flags.

        Returns:
            List of command-line arguments.
        """
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

        returncode, stdout, stderr = _run_and_capture([
            *self.exe_path.split(), xyz_path,
            '-I', str(input_path), *command_line
        ], self.workdir, output, output)

        if returncode != 0:
            raise RuntimeError('error while running xtb: {}'.format(stderr))

        self.clear_workdir()

        total_energy = _find_float('TOTAL ENERGY', stdout, 22, 43, 'xtb')

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
            f.write('$thermo\n  temp={}\n'.format(T))
            if self.use_bhess:
                f.write('  imagthr={}\n  scale={}\n  sthr={}\n'.format(self.imagthr, self.scale, self.sthr))

            self._write_input(geometry, add_solvent, f)

            f.write('$end')

        returncode, stdout, stderr = _run_and_capture([
            *self.exe_path.split(), xyz_path,
            *command_line,
            '--bhess' if self.use_bhess else '--ohess',
            '-I', str(input_path)
        ], self.workdir, output, output)

        if returncode != 0:
            raise RuntimeError('error while running xtb: {}'.format(stderr))

        self.clear_workdir()

        total_energy = _find_float('TOTAL ENERGY', stdout, 22, 43, 'xtb')
        total_free_energy = _find_float('TOTAL FREE ENERGY', stdout, 23, 43, 'xtb')

        if add_solvent:
            gsolv = _find_float('-> Gsolv', stdout, 9, 42, 'xtb')

            return total_energy - gsolv, total_energy, total_free_energy
        else:
            return total_energy, total_free_energy

    def optimize_geometry(
            self, geometry: Molecule, add_solvent: bool = False, output: TextIO = sys.stdout, maxcycle: int = -1,
            optlevel: str = 'normal'
    ) -> Molecule:
        xyz_path = _make_temp_xyz(self.workdir, geometry)
        command_line = self._make_command_line(geometry, add_solvent)

        input_path = self.workdir / 'input.xtb'

        with input_path.open('w') as f:
            f.write('$opt\n  optlevel={}\n'.format(optlevel))
            if maxcycle > 0:
                f.write('  maxcycle={}\n'.format(maxcycle))

            self._write_input(geometry, add_solvent, f)

            f.write('$end')

        returncode, stdout, stderr = _run_and_capture([
            *self.exe_path.split(), xyz_path, *command_line,
            '--opt', '-I', str(input_path)
        ], self.workdir, output, output)

        if returncode != 0:
            raise RuntimeError('error while running xtb: {}'.format(stderr))

        total_energy = _find_float('TOTAL ENERGY', stdout, 26, 43, 'xtb')
        gnorm = _find_float('GRADIENT NORM', stdout, 26, 43, 'xtb')

        with (self.workdir / 'xtbopt.xyz').open() as f:
            new_geometry = Molecule.from_xyz(
                f, geometry.charge, geometry.multiplicity, total_energy, gnorm, name=geometry.name)

        position = stdout.rfind('GEOMETRY OPTIMIZATION CONVERGED')
        if position > 0:
            new_geometry.converged = True

        self.clear_workdir()
        return new_geometry


BORH_TO_ANG = 5.29177210544e-1

#: Convergence thresholds: (E, grms, gmax, drms, dmax)
OPTLEVELS = {
    'loose': (.5e-4, .4e-2, .6e-2, .8e-2, 1.2e-2),
    'normal': (.5e-5, .1e-2, .15e-2, .5e-2, .75e-2),
    'tight': (.1e-5, .8e-3, 1.2e-3, .1e-2, .15e-2),
}


class VlxDriver(QMDriver):
    """Interface for the VeloxChem program.

    Supports arbitrary functionals and basis sets with optional solvation (CPCM or SMD).
    Provides support for prescreening, screening, and full optimizations.

    Attributes:
        exe_path: Path to VeloxChem executable.
        solvatation_model: Solvation model ('cpcm' or 'smd'). `None` for gas-phase.
        solvent: Solvent parameter (epsilon for CPCM, solvent name for SMD).
    """

    def __init__(
            self, workdir: pathlib.Path, exe_path: str | pathlib.Path, method: str, basis: str,
            solvatation_model: Optional[str] = None, solvent: Optional[str | float] = None,
    ):
        """Initialize a VlxDriver.

        Args:
            workdir: Path to working directory.
            exe_path: Path to VeloxChem executable.
            method: Functional name (e.g., 'pbe0', 'rcam-b3lyp').
            basis: Basis set (e.g., 'def2-svp', 'def2-tzvpd').
            solvatation_model: Solvation model ('cpcm' or 'smd'). Defaults to `None`.
            solvent: Solvent specification. Defaults to `None`.
        """
        super().__init__(workdir, method, basis)

        self.exe_path = exe_path
        self.solvatation_model = solvatation_model
        self.solvent = solvent

    def __str__(self):
        """Return driver name with method, basis, and solvation info.

        Returns:
            String like 'VlxDriver[pbe0/def2-svp]' or 'VlxDriver[rcam-b3lyp/def2-tzvpd,cpcm(18.5)]'.
        """
        return 'VlxDriver[{}/{}{}]'.format(
            self.method, self.basis,
            '' if self.solvatation_model is None else ',{}({})'.format(self.solvatation_model, self.solvent)
        )

    def _write_input(self, geometry: Molecule, add_solvent: bool, f: TextIO):
        """Write VeloxChem input file.

        Args:
            geometry: `Molecule` with charge and multiplicity.
            add_solvent: If `True`, include solvation settings.
            f: File object to write to.
        """
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
            [*self.exe_path.split(), str(input_path)], self.workdir, output)

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
            self, geometry: Molecule, add_solvent: bool = False, output: TextIO = sys.stdout, maxcycle: int = -1,
            optlevel: str = 'normal'
    ) -> Molecule:
        input_path = self.workdir / 'input.vlx'

        conv_energy, conv_grms, conv_gmax, conv_drms, conv_dmax = OPTLEVELS[optlevel]

        with input_path.open('w') as f:
            f.write('@jobs\ntask: optimize\n@end\n')

            f.write('@optimize\nconv_energy: {}\nconv_grms: {}\nconv_gmax: {}\nconv_drms: {}\nconv_dmax: {}\n'.format(
                conv_energy, conv_grms, conv_gmax, conv_drms, conv_dmax
            ))

            if maxcycle > 0:
                f.write('max_iter: {}\nconv_maxiter: true\n'.format(maxcycle))

            f.write('@end\n')

            self._write_input(geometry, add_solvent, f)

        returncode, stdout, stderr = _run_and_capture(
            [*self.exe_path.split(), str(input_path)], self.workdir, output)

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
            de < conv_energy
            and drms < conv_drms
            and dmax < conv_dmax
            and grms < conv_grms
            and gmax < conv_gmax
        )

        # extract last geometry
        with h5py.File(self.workdir / 'input.h5') as f:
            new_position = f['atom_coordinates'][:] * BORH_TO_ANG
            new_geometry = Molecule(
                geometry.symbols, new_position,
                geometry.charge, geometry.multiplicity,
                total_energy, grms * numpy.sqrt(len(geometry)), is_converged, name=geometry.name
            )

        self.clear_workdir()
        return new_geometry


class OrcaDriver(QMDriver):
    """Interface for the Orca quantum chemistry program.

    Supports arbitrary functionals and basis sets with optional solvation (CPCM or SMD).
    Provides support for single-point energy calculations and geometry optimizations.

    Attributes:
        exe_path: Path to Orca executable.
        solvatation_model: Solvation model ('cpcm' or 'smd'). `None` for gas-phase.
        solvent: Solvent parameter (epsilon for CPCM, solvent name for SMD).
        nprocs: Number of processors for parallel calculation. Defaults to 1.
    """

    def __init__(
            self, workdir: pathlib.Path, exe_path: str | pathlib.Path, method: str, basis: str,
            solvatation_model: Optional[str] = None, solvent: Optional[str | float] = None, nprocs: int = 1
    ):
        """Initialize a OrcaDriver.
        """
        super().__init__(workdir, method, basis)

        self.exe_path = exe_path
        self.solvatation_model = solvatation_model
        self.solvent = solvent
        self.nprocs = nprocs

    def __str__(self):
        """Return driver name with method, basis, and solvation info.

        Returns:
            String like 'OrcaDriver[pbe0/def2-svp]'.
        """
        return 'OrcaDriver[{}/{}{}]'.format(
            self.method, self.basis,
            '' if self.solvatation_model is None else ',{}({})'.format(self.solvatation_model, self.solvent)
        )

    def _write_input(self, geometry: Molecule, add_solvent: bool, extra_keywords: Optional[str], f: TextIO):
        """Write Orca input file with method, basis, and solvation settings.

        Generates Orca input specification including the functional, basis set, charge,
        multiplicity, geometry, and optional solvation model.

        Args:
            geometry: `Molecule` with charge and multiplicity.
            add_solvent: If `True`, include solvation model in input.
            extra_keywords: Additional keywords for the Orca input file (e.g., 'opt' for optimization).
            f: File object to write to.

        Raises:
            RuntimeError: If solvation model is requested but not set, or if unknown solvation model is specified.
        """

        if add_solvent:
            if self.solvatation_model is None:
                raise RuntimeError('solvatation model is not set')

            if self.solvatation_model.lower() == 'cpcm':
                solvent = 'cpcm({})'.format(self.solvent)
            elif self.solvatation_model.lower() == 'smd':
                solvent = 'smd({})'.format(self.solvent)
            else:
                raise RuntimeError('unknown solvation model for orca `{}`'.format(self.solvatation_model))

            extra_keywords = (extra_keywords + ' ' if extra_keywords else '') + solvent

        f.write('! {} {} {}\n'.format(self.method, self.basis, extra_keywords if extra_keywords else ''))

        if self.nprocs > 1:
            f.write('%pal nprocs {}\nend\n'.format(self.nprocs))

        f.write('*xyzfile {} {} input.xyz\n'.format(geometry.charge, geometry.multiplicity))

    def get_energy(
            self, geometry: Molecule, add_solvent: bool = False, output: TextIO = sys.stdout
    ) -> float | tuple[float, float]:
        _make_temp_xyz(self.workdir, geometry)
        input_path = self.workdir / 'input.orca'

        with input_path.open('w') as f:
            self._write_input(geometry, add_solvent, None, f)

        returncode, stdout, stderr = _run_and_capture(
            [*self.exe_path.split(), str(input_path)], self.workdir, output)

        if returncode != 0:
            raise RuntimeError('error while running orca: {}'.format(stderr))

        self.clear_workdir()

        total_energy = _find_float('FINAL SINGLE POINT ENERGY', stdout, 26, 49, 'orca')
        if add_solvent:
            gsolv = _find_float('CPCM Dielectric    :', stdout, 22, 47, 'orca')
            if self.solvatation_model.lower() == 'smd':
                gsolv += _find_float('Free-energy (cav+disp)  :', stdout, 26, 52, 'orca')

            return total_energy - gsolv, total_energy

        else:
            return total_energy

    def optimize_geometry(
            self, geometry: Molecule, add_solvent: bool = False, output: TextIO = sys.stdout, maxcycle: int = -1,
            optlevel: str = 'normal'
    ) -> Molecule:
        _make_temp_xyz(self.workdir, geometry)
        input_path = self.workdir / 'input.orca'

        conv_energy, conv_grms, conv_gmax, conv_drms, conv_dmax = OPTLEVELS[optlevel]

        with input_path.open('w') as f:
            self._write_input(geometry, add_solvent, 'opt', f)

            f.write('%geom\n')
            if maxcycle > 0:
                f.write('  MaxIter {}\n'.format(maxcycle + 1))
            f.write('  TolE {}\n'.format(conv_energy))
            f.write('  TolRMSG {}\n'.format(conv_grms))
            f.write('  TolMaxG {}\n'.format(conv_gmax))
            f.write('  TolRMSD {}\n'.format(conv_drms))
            f.write('  TolMaxD {}\n'.format(conv_dmax))
            f.write('end\n')

        returncode, stdout, stderr = _run_and_capture(
            [*self.exe_path.split(), str(input_path)], self.workdir, output)

        if returncode != 0:
            raise RuntimeError('error while running orca: {}'.format(stderr))

        total_energy = _find_float('FINAL SINGLE POINT ENERGY', stdout, 26, 49, 'orca')
        gnorm = _find_float('Norm of the Cartesian gradient', stdout, 39, 55, 'orca')

        new_positions = geometry.positions.copy()

        position = stdout.rfind('CARTESIAN COORDINATES (ANGSTROEM)')

        if position < 0:
            raise RuntimeError('error while running orca: unable to find `CARTESIAN COORDINATES (ANGSTROEM)` in output')

        for i in range(len(geometry)):
            line = stdout[position + 68 + 42 * i: position + 68 + 42 * (i + 1)].split()
            new_positions[i] = [float(line[1]), float(line[2]), float(line[3])]

        new_geometry = Molecule(
            geometry.symbols, new_positions,
            geometry.charge, geometry.multiplicity,
            total_energy, gnorm, False, name=geometry.name
        )

        position = stdout.rfind('THE OPTIMIZATION HAS CONVERGED')
        if position > 0:
            new_geometry.converged = True

        self.clear_workdir()

        return new_geometry
