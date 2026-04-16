PyMERK - Python Molecular Ensembles RanKing
=============================================

Purpose
-------

**PyMERK** is a Python tool for ranking and filtering ensembles of molecular systems using a `CENSO <https://xtb-docs.readthedocs.io/en/latest/CENSO_docs/censo.html>`_-like workflow.
It processes multiple molecular conformers through progressive filtering stages:

1. **Prescreening**: Fast single-point calculations with cheap solvation correction (xTB) to eliminate high-energy candidates;
2. **Screening**: Refined DFT with larger basis sets, solvation and thermochemical corrections (the latter with xTB);
3. **Optimization**: Full geometry optimization with best methods;
4. **Refinement**: Boltzmann population filtering at desired temperature.

The tool supports:

- Multiple quantum chemistry programs;
- Flexible energy threshold filtering;
- Solvent and thermal corrections;
- Macro-cycle optimization protocol;
- Comprehensive conformer statistics.

About
-----

**PyMERK** is developed by `Pierre Beaujean <https://pierrebeaujean.net>`_, 
who holds a Ph.D. in quantum chemistry from the `University of Namur <https://unamur.be>`_ (Belgium).

The project was created to automate and streamline the ranking of molecular conformer ensembles, a common task in computational chemistry workflows. 
It implements a CENSO-like approach adapted for various quantum chemistry programs and is actively used in research on molecular systems and solvation effects.

For questions or discussions, please feel free to open an issue on GitHub or contact the author directly.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   usage
   api
   contributing