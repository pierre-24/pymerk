PyMERK - Python Molecular Ensembles RanKing
=============================================

Purpose
-------

**PyMERK** is a Python tool for ranking and filtering ensembles of molecular systems using a CENSO-like workflow.
It processes multiple molecular conformers through progressive filtering stages:

1. **Prescreening**: Fast single-point calculations using PBE0/def2-SVP to eliminate high-energy candidates
2. **Screening**: Refined DFT with larger basis sets (RCAM-B3LYP/def2-TZVPD) and solvation corrections
3. **Optimization**: Full geometry optimization with best methods (WB97M-D4 or custom)
4. **Refinement**: Boltzmann population filtering at desired temperature

The tool supports:

- Multiple quantum chemistry programs (xTB for semiempirical calculations, VeloxChem for DFT)
- Flexible energy threshold filtering
- Solvent and thermal corrections
- Macro-cycle optimization protocols
- Comprehensive conformer statistics (RMSD matrices, relative energies)


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   usage
   api
   contributing


About
-----

**PyMERK** is developed by `Pierre Beaujean <https://pierrebeaujean.net>`_, 
who holds a Ph.D. in quantum chemistry from the `University of Namur <https://unamur.be>`_ (Belgium).

The project was created to automate and streamline the ranking of molecular conformer ensembles, a common task in computational chemistry workflows. 
It implements a CENSO-like approach adapted for various quantum chemistry programs and is actively used in research on molecular systems and solvation effects.

For questions or discussions, please feel free to open an issue on GitHub or contact the author directly.
