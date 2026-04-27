Tutorial 2: Calcium solvation with Orca
=======================================

The goal of this tutorial is to run the **pyMERK** program on a set of 10 conformers of :math:`\ce{Ca_2(BH_4)_4(THF)_2}` in order to identify and optimize the most relevant structures.
All quantum chemistry calculations in this workflow are performed using `Orca <https://orca-manual.mpi-muelheim.mpg.de/>`_.

You can download the initial structures here: :download:`goat.finalensemble.xyz`.
These geometries were generated from a `GOAT calculation <https://orca-manual.mpi-muelheim.mpg.de/contents/structurereactivity/goat.html>`_ using rigid fragments.

.. note::

   Running the full workflow may take several **hours** on a standard workstation.
   The reference calculation presented here was performed on a supercomputer.
   Precomputed output files are provided if you prefer to follow the tutorial without running the calculations.

Setting up & running
--------------------

For this tutorial, we use the following input file: :download:`input.toml`.
Apart from the choice of solvent (tetrahydrofuran instead of water), most parameters are left at their default values as described in :doc:`../../usage`.

This tutorial assumes that `xTB <https://xtb-docs.readthedocs.io/en/latest/>`_ and Orca are properly installed and available in your environment.
Before running the workflow, **verify and update the paths** to xTB (:pymkw:`paths.xtb`) and Orca (:pymkw:`paths.orca`):

.. literalinclude:: input.toml
   :language: toml
   :start-at: [paths]
   :end-before: [general]

.. tip::

   To enable MPI parallelization in Orca, you must provide the full path to the Orca executable.
   See `Parallel instructions for ORCA <https://orca-manual.mpi-muelheim.mpg.de/contents/essentialelements/parallel.html>`_ for details.
   You can then adjust :pymkw:`paths.orca_nprocs` to control the number of processes.

Run the workflow with:

.. code:: bash

   pymerk_run -i input.toml goat.finalensemble.xyz -o final_ensemble.xyz

Once the workflow starts, progress can be monitored from the standard output.
If you prefer not to run the calculation, you can use the provided :download:`run.log`.

Analyzing the output
--------------------

Although this system is significantly larger than in the previous tutorial, the use of the `r²SCAN-3c approach <https://chemrxiv.org/doi/full/10.26434/chemrxiv.13333520.v2>`_ in Orca greatly accelerates the screening and optimization stages.

The input section corresponding to the optimization stage is:

.. literalinclude:: input.toml
   :language: toml
   :start-at: [optimization]
   :end-before: [refinement]

This setup uses the r²SCAN-3c/def2-mTZVPP method in tetrahydrofuran, with solvation described via SMD.

The corresponding output is:

.. literalinclude:: run.log
   :language: text
   :start-at: Macrocycle 1
   :end-at: * Done

A total of 11 macrocycles are required to reduce the set to the 5 lowest-energy conformers out of the 9 retained after the previous stages.
The selection effectively begins at macrocycle 5, when the gradient norms fall below :pymkw:`optimization.gradthr`.
This illustrates the efficiency of the macrocycle procedure in discarding high-energy conformers early.

The refinement stage further reduces the ensemble using a higher level of theory (wb97x-d3/def2-TZVPP):

.. literalinclude:: run.log
   :language: text
   :start-at: 4_refinement
   :end-at: * Done

At this stage, the conformer set is reduced to 3 structures, since the two first ones has negligible Boltzmann weight at 298 K.

The final geometries are available in :download:`4_refinement.selected.xyz`.