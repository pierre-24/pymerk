Installation & usage
====================

**Note:** At the moment, only `VeloxChem <https://veloxchem.org/>`_ and `Orca <https://orca-manual.mpi-muelheim.mpg.de/>`_ can be used as the QM driver.

Installation
------------

To install this package, you need:

- Python 3.11 or higher (with ``pip`` and ``venv``, but those are generally default), and
- `xTB <https://xtb-docs.readthedocs.io/en/latest/>`_ (mandatory), and, optionally,
- `VeloxChem <https://veloxchem.org/>`_ and/or `Orca <https://orca-manual.mpi-muelheim.mpg.de/>`_.

Use:

.. code-block:: bash

    # for the stable version
    pip3 install git+https://github.com/pierre-24/pymerk.git@v0.2.1
    # for the latest version
    pip3 install git+https://github.com/pierre-24/pymerk.git@dev

Note: as this script install programs, you might need to add them (such as ``$HOME/.local/bin``, if you use ``--user``) to your ``$PATH``.

Verify the installation:

.. code-block:: bash

    pymerk_run --help


Usage
-----

PyMERK is controlled entirely through TOML configuration files.
The main entry point is ``pymerk_run`` (see below).

Configuration File
__________________

Create a TOML configuration file to define the workflow parameters and override default values where needed.

A template file containing all default settings can be generated with:

.. code:: bash

   pymerk_config > input.toml

This file can then be edited to suit your system and computational requirements.

The overall structure closely follows that of the `CENSO .censo2rc <https://xtb-docs.readthedocs.io/en/latest/CENSO_docs/censorc.html>`_ configuration file.
Some keywords have been adapted, added, or removed to reflect pyMERK-specific features.

General Settings
~~~~~~~~~~~~~~~~

.. pymkwdef:: temperature
   :section: general
   :type: float
   :default: 298.15

   Temperature in Kelvin used for all thermochemical and Boltzmann population calculations.

.. pymkwdef:: evaluate_rrho
   :section: general
   :type: bool
   :default: true

   Enable or disable the evaluation of RRHO thermochemical contributions.

.. pymkwdef:: sm_rrho
   :section: general
   :type: str
   :default: "gbsa"

   Solvation model used for RRHO corrections, typically ``gbsa`` or ``alpb``.

.. pymkwdef:: imagthr
   :section: general
   :type: float
   :default: -100.0

   Threshold for imaginary frequencies in cm⁻¹ below which modes are treated as invalid or ignored.

.. pymkwdef:: sthr
   :section: general
   :type: float
   :default: 50.0

   Low-frequency threshold in cm⁻¹ used in the treatment of vibrational modes.

.. pymkwdef:: solvent
   :section: general
   :type: str
   :default: "h2o"

   Identifier of the solvent used in the calculations, for example ``h2o`` or ``dmso``.

.. pymkwdef:: gas_phase
   :section: general
   :type: bool
   :default: false

   If set to true, all solvation contributions are ignored and calculations are performed in the gas phase.

.. note::

    If the value of :pymkw:`general.solvent` is not the same in ``xtb`` and the QM driver, you can add ``alternate_solvent`` (:pymkw:`screening.alternate_solvent`, :pymkw:`optimization.alternate_solvent`, and :pymkw:`refinement.alternate_solvent`) to latter stages to define the equivalent solvent.
    For example, for ``"thf"``, you need to set ``alternate_solvent = 'tetrahydrofuran'`` with VeloxChem.

Program Paths
~~~~~~~~~~~~~

.. pymkwdef:: xtb
   :section: paths
   :type: str
   :default: ""

   Path to the xTB executable.

.. pymkwdef:: vlx
   :section: paths
   :type: str
   :default: ""

   Path to the VeloxChem executable.
   This setting is used when ``prog = "vlx"`` in one of the stage.

.. pymkwdef:: orca
   :section: paths
   :type: str
   :default: ""

   Path to the Orca executable.
   Note that in order to use multiple processes (:pymkw:`paths.orca_nprocs`), you need to provide the full path (see `Parallel instructions for ORCA <https://orca-manual.mpi-muelheim.mpg.de/contents/essentialelements/parallel.html>`_).
   This setting is used when ``prog = "orca"`` in one of the stage.

.. pymkwdef:: orca_nprocs
   :section: paths
   :type: int
   :default: 1

   Number of processes used for Orca calculations, with a maximum of 64.
   This setting is used when ``prog = "orca"`` in one of the stage.


.. note::

    You can also set runners and default options, by using, *e.g.*,

    .. code-block:: toml

        [paths]
        xtb = "xtb -v"    # more verbose output with xTB
        vlx = "srun vlx"  # run VeloxChem via srun


Prescreening Stage (fast single-points)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. pymkwdef:: prog
   :section: prescreening
   :type: str
   :default: "orca"

   Program used for electronic structure calculations, typically ``vlx`` (VeloxChem) or ``orca``.

.. pymkwdef:: func
   :section: prescreening
   :type: str
   :default: "pbe d3"

   Exchange-correlation functional used for DFT calculations.

.. pymkwdef:: basis
   :section: prescreening
   :type: str
   :default: "def2-sv(p)"

   Basis set used during the prescreening stage.

.. pymkwdef:: gfnv
   :section: prescreening
   :type: str
   :default: "gfn2"

   xTB variant used for auxiliary contributions, such as ``gfn1``, ``gfn2``, or ``gfnff``.

.. pymkwdef:: threshold
   :section: prescreening
   :type: float
   :default: 4.0

   Energy threshold in kcal/mol used to retain conformers relative to the lowest-energy structure.

Screening Stage (refined calculations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. pymkwdef:: prog
   :section: screening
   :type: str
   :default: "orca"

   Program used for electronic structure calculations.

.. pymkwdef:: func
   :section: screening
   :type: str
   :default: "r2scan-3c"

   Exchange-correlation functional used for the screening stage.

.. pymkwdef:: basis
   :section: screening
   :type: str
   :default: "def2-mTZVPP"

   Basis set used during the screening stage.

.. pymkwdef:: sm
   :section: screening
   :type: str
   :default: "smd"

   Solvation model used in this stage, such as ``smd``, ``cpcm``, or ``gbsa``.

.. pymkwdef:: alternate_solvent
   :section: screening
   :type: Optional[str | int]
   :default: null

   Alternate solvent name to be used by the QM driver, to be provided if it does not match :pymkw:`general.solvent`.

.. pymkwdef:: gfnv
   :section: screening
   :type: str
   :default: "gfn2"

   xTB variant used for auxiliary energy corrections.

.. pymkwdef:: threshold
   :section: screening
   :type: float
   :default: 3.5

   Relative energy cutoff in kcal/mol used to retain conformers.

.. pymkwdef:: gsolv_included
   :section: screening
   :type: bool
   :default: false

   If set to true, solvation effects are included directly in the QM driver energies, otherwise they are computed separately using xTB.

Optimization Stage (full geometry optimization)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. note::

    At this stage, the solvent must be included directly by the QM driver, as ``xTB``-gradient correction is not yet possible.

.. pymkwdef:: prog
   :section: optimization
   :type: str
   :default: "orca"

   Program used to perform geometry optimizations.

.. pymkwdef:: func
   :section: optimization
   :type: str
   :default: "r2scan-3c"

   Exchange-correlation functional used during optimization.

.. pymkwdef:: basis
   :section: optimization
   :type: str
   :default: "def2-mTZVPP"

   Basis set used for geometry optimizations.

.. pymkwdef:: sm
   :section: optimization
   :type: str
   :default: "cpcm"

   Solvation model used during optimization.

.. pymkwdef:: alternate_solvent
   :section: optimization
   :type: Optional[str | int]
   :default: null

   Alternate solvent name to be used by the QM driver, to be provided if it does not match :pymkw:`general.solvent`.

.. pymkwdef:: optlevel
   :section: optimization
   :type: str
   :default: "normal"

   Optimization convergence level. At the moment, only ``loose``, ``normal``, or ``tight`` are supported.

.. pymkwdef:: gfnv
   :section: optimization
   :type: str
   :default: "gfn2"

   xTB variant used for RRHO thermochemical corrections.

.. pymkwdef:: threshold
   :section: optimization
   :type: float
   :default: 3.0

   Energy threshold in kcal/mol used to discard conformers after optimization.

.. pymkwdef:: macrocycles
   :section: optimization
   :type: bool
   :default: true

   Enable or disable the macrocycle optimization protocol.

.. pymkwdef:: gradthr
   :section: optimization
   :type: float
   :default: 0.01

   Gradient norm threshold in atomic units below which conformers are compared and filtered.

.. pymkwdef:: maxcyc
   :section: optimization
   :type: int
   :default: 200

   Maximum number of optimization iterations per conformer.

.. pymkwdef:: optcycles
   :section: optimization
   :type: int
   :default: 8

   Number of microcycles per macrocycle when the macrocycle protocol is enabled.


Refinement Stage (Boltzmann population filtering)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. pymkwdef:: prog
   :section: refinement
   :type: str
   :default: "orca"

   Program used for final single-point energy calculations.

.. pymkwdef:: func
   :section: refinement
   :type: str
   :default: "wb97m-v"

   Exchange-correlation functional used for high-accuracy refinement.

.. pymkwdef:: basis
   :section: refinement
   :type: str
   :default: "def2-TZVP"

   Basis set used during the refinement stage.

.. pymkwdef:: sm
   :section: refinement
   :type: str
   :default: "smd"

   Solvation model used in the refinement calculations.

.. pymkwdef:: alternate_solvent
   :section: refinement
   :type: Optional[str | int]
   :default: null

   Alternate solvent name to be used by the QM driver, to be provided if it does not match :pymkw:`general.solvent`.

.. pymkwdef:: gfnv
   :section: refinement
   :type: str
   :default: "gfn2"

   xTB variant used for auxiliary corrections.

.. pymkwdef:: threshold
   :section: refinement
   :type: float
   :default: 0.95

   Cumulative Boltzmann population threshold (between 0 and 1) used to select the final set of conformers.

.. pymkwdef:: gsolv_included
   :section: refinement
   :type: bool
   :default: false

   If set to true, solvation effects are included directly in the QM driver energies, otherwise they are computed separately using xTB.

Running PyMERK
______________

Execute the workflow with:

.. code-block:: bash

   pymerk_run input_ensemble.xyz -i config.toml -o output_ensemble.xyz \
       -c <charge> -m <multiplicity> -w <workdir>

Where:

- ``input_ensemble.xyz`` is a multi-structure XYZ file containing the initial conformers.
- ``config.toml`` is the TOML configuration file defining the workflow parameters.
- ``output_ensemble.xyz`` is the output XYZ file containing the final selected conformers.
- ``<charge>`` is the total molecular charge (default: 0).
- ``<multiplicity>`` is the spin multiplicity (default: 1).
- ``<workdir>`` is the *working directory* for the files generated by the different programs as well as intermediate results (default: .).

Only the input structure file is required, and all other arguments are optional.

.. tip::

   If no output file is specified, a default name will be generated automatically.

Output
______

The workflow processes the conformer ensemble through all enabled stages and produces, in the standard output:

- **Progress**;
- **Energy summaries** reporting relative energies for all evaluated structures; and
- **RMSD matrices** (in Ångströms) for the conformers that pass each stage.

In the *working directory*, you will also find:

- **Log files** containing the output of each program, and
- **Filtered ensembles** at each stage as ``.xyz`` files containing only retained conformers.

These outputs allow you to track how the conformer set is progressively reduced and refined.
