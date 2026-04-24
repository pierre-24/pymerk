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

Before running the program, **verify and update the paths** to xTB (:pymkw:`paths.xtb`) and VeloxChem (:pymkw:`paths.vlx`):

.. literalinclude:: input.toml
   :language: toml
   :start-at: [paths]

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

A pyMERK workflow (loosely inspired by `CENSO <https://xtb-docs.readthedocs.io/en/latest/CENSO_docs/censo.html>`_) is organized into four successive stages.
These stages are pre-screening, screening, optimization, and refinement.

Each stage is carried out at a defined level of theory set in the input file.
Electronic energies are computed with the main QM engine (here, VeloxChem) and complemented by xTB contributions when required.

General settings are defined in the ``[general]`` section of the TOML input:

.. literalinclude:: input.toml
   :language: toml
   :start-at: [general]
   :end-before: [prescreening]

In this example, the solvent is set to water.
The option :pymkw:`general.gas_phase` is set to ``false`` (default) means that solvent effects are included in all reported energies.
Depending on the stage, these contributions may come from either VeloxChem or xTB.
The option :pymkw:`general.evaluate_rrho` is set to ``true`` (default) enables thermochemical corrections via the `SPH correction <https://pubs.acs.org/doi/10.1021/acs.jctc.0c01306>`_ using xTB.

Additional options are described in :doc:`../usage`.

Pre-screening step
~~~~~~~~~~~~~~~~~~

The **pre-screening** stage provides a fast and inexpensive way to discard clearly unfavorable conformers.

Input section:

.. literalinclude:: input.toml
   :language: toml
   :start-at: [prescreening]
   :end-before: [screening]

Thanks to :pymkw:`prescreening.basis`, a reduced basis set is used to lower computational cost.
Conformers within 4 kcal/mol of the lowest-energy structure (:pymkw:`prescreening.threshold`) are retained.

The output starts with a summary of the computational setup:

.. literalinclude:: run.log
   :language: text
   :lines: 1-15
   :append: ...

The notation ``elec=MAIN, gsolv=AUX, gtrv=NONE`` specifies how energy components are evaluated.
The ``elec`` term corresponds to the electronic energy computed by the main QM driver (VeloxChem) at the PBE0/def2-SV(P) level.
The ``gsolv`` term corresponds to the solvation correction computed by the auxiliary driver (xTB) using GFN2-xTB/GBSA(water).
The ``gtrv`` term indicates that no thermochemical contribution is included at this stage.

Relative energies are then reported:

.. literalinclude:: run.log
   :language: text
   :lines: 26-39

Conformers above the threshold are removed, and in this example 3 conformers are discarded.
The cutoff can be tuned via the :pymkw:`prescreening.threshold` parameter in the input file.

An RMSD matrix between the remaining conformers is also printed:

.. literalinclude:: run.log
   :language: text
   :lines: 43-53

These values provide a measure of structural diversity.
In this case, the conformers are all significantly different from one another.

Screening step
~~~~~~~~~~~~~~

The **screening** stage refines the selection using a higher level of theory and includes thermochemical corrections.

Input section:

.. literalinclude:: input.toml
   :language: toml
   :start-at: [screening]
   :end-before: [optimization]

A reduced basis set is still used to balance cost and accuracy.
The energy threshold (via :pymkw:`screening.threshold`) is slightly stricter and set to 3.5 kcal/mol by default.

The output is similar to the previous stage:

.. literalinclude:: run.log
   :language: text
   :start-at: 2_screening
   :end-at: * Done

The notation ``elec=MAIN, gsolv=AUX, gtrv=AUX`` indicates how energy contributions are evaluated.
The ``elec`` term corresponds to the electronic energy computed by VeloxChem at the rcam-b3lyp/def2-svpd level.
The ``gsolv`` term corresponds to the solvation correction computed by xTB at the GFN2-xTB/GBSA level.
The ``gtrv`` term corresponds to thermochemical corrections computed by xTB via the SPH scheme.

At this stage, one additional conformer is discarded.

Optimization step
~~~~~~~~~~~~~~~~~

The **optimization** stage performs geometry optimizations of the retained conformers.

Input section:

.. literalinclude:: input.toml
   :language: toml
   :start-at: [optimization]
   :end-before: [refinement]

The level of theory is rcam-b3lyp/def2-svpd.
Solvent effects must now be included directly in the QM calculation, for example via CPCM in VeloxChem.
The dielectric constant is set using the :pymkw:`optimization.alternate_solvent` parameter.

The *macrocycle* procedure is used by default.
Up to 8 optimization cycles are performed per conformer, controlled by :pymkw:`optimization.optcycles`.
After each macrocycle, converged structures based on :pymkw:`optimization.gradthr` are compared.
Conformers above the energy threshold of 3 kcal/mol (value of :pymkw:`optimization.threshold`) are discarded early.

The output for this stage is:

.. literalinclude:: run.log
   :language: text
   :start-at: Macrocycle 1
   :end-at: * Done

In this example, all 8 conformers converge after 5 macrocycles and are retained.
The optimized geometries are available in :download:`3_optimize.selected.xyz`.

Refinement step
~~~~~~~~~~~~~~~

The **refinement** stage evaluates the conformer ensemble at a higher level of theory and selects structures based on their Boltzmann population.

.. math::

   p_i = \frac{e^{-\Delta G_i / RT}}{\sum_j e^{-\Delta G_j / RT}}

Input section:

.. literalinclude:: input.toml
   :language: toml
   :start-at: [refinement]
   :end-before: [paths]

Solvation is handled directly by VeloxChem using the SMD model with :pymkw:`refinement.gsolv_included` set to ``true``.

In this stage, the :pymkw:`refinement.threshold` parameter defines a cumulative Boltzmann population cutoff of 95 percent rather than an energy difference.

The output for this stage is:

.. literalinclude:: run.log
   :language: text
   :start-at: 4_refinement
   :end-at: * Done

Only the conformers required to reach the target population are retained.
In this example, 7 conformers are sufficient.

The final structures are available in :download:`4_refinement.selected.xyz`.