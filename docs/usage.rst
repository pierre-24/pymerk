Installation & usage
====================

Installation
------------

To install this package, you need:

- Python 3.11 or higher (with ``pip```), and
- xTB and/or VeloxChem quantum chemistry packages (optional but recommended).

Use:

.. code-block:: bash

    pip3 install git+https://github.com/pierre-24/pymerk.git


Note: as this script install programs, you might need to add them (such as ``$HOME/.local/bin``, if you use ``--user``) to your ``$PATH``.

Verify the installation:

.. code-block:: bash

    pymerk_run --help


Usage
-----

PyMERK is controlled entirely through TOML configuration files. The main entry point is ``pymerk_run``.

Configuration File
__________________

Create a TOML configuration file specifying the workflow parameters.

General Settings
~~~~~~~~~~~~~~~~

.. code-block:: toml

    [general]
    temperature = 298.15              # Temperature in Kelvin for calculations
    evaluate_rrho = true              # Calculate RRHO contributions
    sm_rrho = "gbsa"                  # Solvation model for RRHO (gbsa, cpcm, smd)
    imagthr = -100.0                  # Imaginary frequency threshold (cm⁻¹)
    sthr = 50.0                       # Wave number threshold (cm⁻¹)
    solvent = "h2o"                   # Default solvent identifier (h2o, dmso, etc)
    gas_phase = false                 # If true, ignore all solvation corrections

Prescreening Stage (fast single-points)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: toml

    [prescreening]
    prog = "vlx"                      # Program: vlx (VeloxChem) or xtb
    func = "pbe0"                     # Functional for DFT
    basis = "def2-svp"                # Basis set for prescreening
    gfnv = "gfn2"                     # xTB variant for gsolv contributions (gfn1, gfn2, gfnff)
    threshold = 4.0                   # Energy threshold (kcal/mol) to retain candidates

Screening Stage (refined calculations)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: toml

    [screening]
    prog = "vlx"                      # Program: vlx or xtb
    func = "rcam-b3lyp"               # Functional for screening
    basis = "def2-tzvpd"              # Basis set for screening
    sm = "smd"                        # Solvation model (smd, cpcm, gbsa)
    alternate_solvent = null          # Optional: override default solvent
    gfnv = "gfn2"                     # xTB variant for corrections
    threshold = 3.5                   # Relative energy cutoff (kcal/mol)
    gsolv_included = false            # Whether solvation is included in energies or calculated separately

Optimization Stage (full geometry optimization)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: toml

    [optimization]
    prog = "vlx"                      # Program for optimization
    func = "rcam-b3lyp"               # Functional for optimization
    basis = "def2-tzvpd"              # Basis set for optimization
    sm = "cpcm"                       # Solvation model (cpcm, smd, gbsa)
    alternate_solvent = null          # Optional: override solvent (epsilon value for CPCM, name for SMD)
    gfnv = "gfn2"                     # xTB variant for RRHO corrections
    optcycles = 8                     # Microcycles per macrocycle (if using macrocycle protocol)
    maxcyc = 200                      # Maximum optimization iterations
    optlevel = "normal"               # Optimization thoroughness (loose, normal, tight)
    threshold = 3.0                   # Energy threshold (kcal/mol) for filtering after optimization
    gradthr = 0.01                    # Gradient threshold (a.u.) below which energy threshold applies
    macrocycles = true                # Enable macrocycle optimization protocol

Refinement Stage (Boltzmann population filtering)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: toml

    [refinement]
    prog = "vlx"                      # Program for refinement
    func = "wb97m-d4"                 # High-accuracy functional
    basis = "def2-tzvpd"              # Basis set for refinement
    sm = "smd"                        # Solvation model
    alternate_solvent = null          # Optional: override solvent
    gfnv = "gfn2"                     # xTB variant
    threshold = 0.95                  # Cumulative Boltzmann population threshold (0-1)
    gsolv_included = false            # How solvation energy is handled

Program Paths
~~~~~~~~~~~~~

.. code-block:: toml

    [paths]
    xtb = "/usr/local/bin/xtb"        # Path to xtb executable
    vlx = "/usr/local/bin/vlx"        # Path to VeloxChem executable

Example Configuration
_____________________

A minimal working example ``config.toml``:

.. code-block:: toml

    [general]
    temperature = 298.15
    solvent = "water"
    gas_phase = false

    [prescreening]
    prog = "vlx"
    func = "pbe0"
    basis = "def2-svp"
    threshold = 4.0

    [screening]
    prog = "vlx"
    func = "rcam-b3lyp"
    basis = "def2-tzvpd"
    sm = "smd"
    threshold = 3.5

    [optimization]
    prog = "vlx"
    func = "rcam-b3lyp"
    basis = "def2-tzvpd"
    sm = "cpcm"
    alternate_solvent = 80.  # Dielectric constant for CPCM (e.g., 80 for water)
    maxcyc = 200
    threshold = 3.0

    [refinement]
    prog = "vlx"
    func = "wb97m-d4"
    basis = "def2-tzvpd"
    sm = "smd"
    threshold = 0.95

    [paths]
    xtb = "/usr/local/bin/xtb"
    vlx = "/usr/local/bin/vlx"

Running PyMERK
______________

Execute the workflow with:

.. code-block:: bash

    pymerk_run input_ensemble.xyz -i config.toml -o output_ensemble.xyz -c <charge> -m <multiplicity>

Where:

- ``input_ensemble.xyz`` - Multi-structure XYZ file with initial conformers
- ``config.toml`` - TOML configuration file with workflow settings
- ``output_ensemble.xyz`` - Output XYZ file with refined conformers
- ``<charge>`` - Total charge of the molecule
- ``<multiplicity>`` - Spin multiplicity of the molecule

Output
______

The tool processes the ensemble through all enabled stages and generates:

- **Filtered ensembles** at each stage (as ``.xyz``  files containing only the retained conformers).
- **Energy reports** with relative energies for all geometries
- **RMSD matrices** (in Ångströms) for retained conformers
