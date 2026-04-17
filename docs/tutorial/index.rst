Tutorial
========

The goal of this tutorial is to run the **pyMERK** program on a set of 12 conformers of :math:`\ce{(H_2O)_6}`.
You can download the initial structures here: :download:`6H2O_12.xyz`.

.. note::

   Running the full workflow may take several hours on a standard workstation.
   Precomputed output files are also provided if you prefer to follow the tutorial more quickly.

This tutorial assumes that you have `xTB <https://xtb-docs.readthedocs.io/en/latest/>`_ and `VeloxChem <https://veloxchem.org/>`_ properly installed and accessible from your environment.

Setting up & running
--------------------

To run pyMERK, you also need an input file in TOML format to override default parameters if needed.
For this tutorial, we will use: :download:`input.toml`.

The content of this file will be discussed progressively in the analysis section below.

Before running the program, **verify and update the paths** to xTB (``xtb``) and VeloxChem (``vlx``):

.. literalinclude:: input.toml
   :language: toml
   :lines: 19-21

.. tip::

   On HPC systems using SLURM or MPI, you may want to adapt the VeloxChem command, e.g.:

   - ``vlx = 'srun vlx'``
   - ``vlx = 'mpirun vlx'``

   Refer to `VeloxChem cluster documentation <https://veloxchem.org/docs/run-cluster/>`_ for details on parallel execution.

Run the workflow with:

.. code:: bash

   # set up 'OMP_STACKSIZE' for xtb
   export OMP_STACKSIZE=4G

   # define OpenMP parallelization
   export OMP_PLACES=cores
   export OMP_NUM_THREADS=16

   # run the workflow
   pymerk_run -i input.toml 6H2O_12.xyz

.. hint::

   Adjust ``OMP_NUM_THREADS`` according to the number of cores available on your machine.

Once the workflow starts, you can monitor progress in the generated log files.
If you do not want to run it yourself, you can grab :download:`run.log`.

Analyzing the output
--------------------

A `CENSO <https://xtb-docs.readthedocs.io/en/latest/CENSO_docs/censo.html>`_ workflow is divided into four stages:

1. Pre-screening
2. Screening
3. Optimization
4. Refinement

Each stage is performed at a specified level of theory (controlled via the input file).
Energies obtained from the main quantum chemistry program (here, VeloxChem) are corrected using xTB.

General settings are defined in the ``[general]`` section of the TOML input:

.. literalinclude:: input.toml
   :language: toml
   :lines: 1-4

In this example:

- The solvent is set to water.
- ``gas_phase = false`` (default), so solvent contributions are included in all reported energies (depending on the stage, they might added via xTB or VeloxChem).
- ``evaluate_rrho = true`` (default) enables thermochemical corrections via the `SPH correction <https://pubs.acs.org/doi/10.1021/acs.jctc.0c01306>`_ using xTB.

Further options are documented in :doc:`../usage`.

Pre-screening step
~~~~~~~~~~~~~~~~~~

The first stage is the **pre-screening**, a computationally inexpensive step designed to eliminate high-energy conformers early in the workflow.

The corresponding input section is:

.. literalinclude:: input.toml
   :language: toml
   :lines: 5-7

Here, the default basis set (def2-SVP) is replaced with a smaller one to reduce computational cost and accelerate this initial filtering.

The output begins with a summary of the level of theory used:

.. literalinclude:: run.log
   :lines: 1-15
   :append: ...

In this section, the notation ``elec=MAIN, gsolv=AUX, gtrv=NONE`` indicates how different energy contributions are evaluated:

- ``elec``: the electronic energy is computed by the *main* QM driver (VeloxChem), here at the PBE0/def2-SV(P) level.
- ``gsolv``: the solvation contribution is added using the *auxiliary* driver (xTB), at the GFN2-xTB/GBSA(water) level.
- ``gtrv``: no thermochemical correction is included at this stage.

Once all conformer energies have been computed, a summary table reports relative energies with respect to the lowest-energy structure:

.. literalinclude:: run.log
   :lines: 26-39

Based on these values, conformers above a given energy threshold are discarded.
In this example, 3 conformers are removed.
The energy cutoff can be adjusted via the ``threshold`` parameter (in kcal/mol) in the TOML input file.

Finally, the RMSD matrix between the remaining conformers is printed:

.. literalinclude:: run.log
   :lines: 43-53

These values provide a measure of structural similarity.
In this case, the relatively large RMSD values indicate that the retained conformers are structurally distinct, which is desirable before proceeding to the next stage.