Installation & usage
====================

**Note:** At the moment, only `VeloxChem <https://veloxchem.org/>`_ can be used as the QM driver.

Installation
------------

To install this package, you need:

- Python 3.11 or higher (with ``pip`` and ``venv``, but those are generally default), and
- `xTB <https://xtb-docs.readthedocs.io/en/latest/>`_ (mandatory) and `VeloxChem <https://veloxchem.org/>`_ (optional but recommended).

Use:

.. code-block:: bash

    pip3 install git+https://github.com/pierre-24/pymerk.git


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

Create a TOML configuration file specifying the workflow parameters to override their default values.
Its format follows closely the one of the `CENSO .censo2rc <https://xtb-docs.readthedocs.io/en/latest/CENSO_docs/censorc.html>`_ file, but some keywords are missing and others have been added.

General Settings
~~~~~~~~~~~~~~~~

.. code-block:: toml

    [general]
    temperature = 298.15              # Temperature in Kelvin for calculations
    evaluate_rrho = true              # Calculate RRHO contributions
    sm_rrho = "gbsa"                  # Solvation model for RRHO (gbsa, alpb)
    imagthr = -100.0                  # Imaginary frequency threshold (cm⁻¹)
    sthr = 50.0                       # Wave number threshold (cm⁻¹)
    solvent = "h2o"                   # Default solvent identifier (h2o, dmso, etc)
    gas_phase = false                 # If true, ignore all solvation corrections

.. note::

    If the value of ``solvent`` is not the same in ``xtb`` and the QM driver, you can add ``alternate_solvent`` to latter stages to define the equivalent solvent.
    For example, for ``"thf"``, you need to set ``alternate_solvent = 'tetrahydrofuran'`` with VeloxChem.

Program Paths
~~~~~~~~~~~~~

.. code-block:: toml

    [paths]
    xtb = "xtb"        # Path to xtb executable
    vlx = "vlx"        # Path to VeloxChem executable


.. note::

    You can also set runners and default options, by using, *e.g.*,

    .. code-block:: toml

        [paths]
        xtb = "xtb -v"    # more verbose output with xTB
        vlx = "srun vlx"  # run VeloxChem via srun


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
    gfnv = "gfn2"                     # xTB variant for corrections
    threshold = 3.5                   # Relative energy cutoff (kcal/mol)
    gsolv_included = false            # Whether solvation is included in QM driver energies or calculated separately with xtb

Optimization Stage (full geometry optimization)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: toml

    [optimization]
    prog = "vlx"                      # Program for optimization
    func = "rcam-b3lyp"               # Functional for optimization
    basis = "def2-tzvpd"              # Basis set for optimization
    sm = "cpcm"                       # Solvation model
    optlevel = "normal"               # Optimization thoroughness (loose, normal, tight)
    gfnv = "gfn2"                     # xTB variant for RRHO corrections
    threshold = 3.0                   # Energy threshold (kcal/mol) for filtering after optimization

    # for MACROCYLE protocol
    macrocycles = true                # Enable macrocycle optimization protocol
    gradthr = 0.01                    # Gradient threshold (a.u.) below which energy threshold applies
    maxcyc = 200                      # Maximum optimization iterations
    optcycles = 8                     # Microcycles per macrocycle (if using macrocycle protocol)

.. note::

    At this stage, the solvent must be included directly by the QM driver, as ``xTB``-gradient correction is not yet possible.


Refinement Stage (Boltzmann population filtering)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: toml

    [refinement]
    prog = "vlx"                      # Program for refinement
    func = "wb97m-d4"                 # High-accuracy functional
    basis = "def2-tzvpd"              # Basis set for refinement
    sm = "smd"                        # Solvation model
    gfnv = "gfn2"                     # xTB variant
    threshold = 0.95                  # Cumulative Boltzmann population threshold (0-1)
    gsolv_included = false            # Whether solvation is included in QM driver energies or calculated separately with xtb

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
