# `pymerk` -- Python Molecular Ensembles RanKing

A [CENSO](https://xtb-docs.readthedocs.io/en/latest/CENSO_docs/censo.html)-like approach to rank ensembles of molecular systems using quantum chemistry calculations.

**Note:** early release, paint is still fresh!

**Note:** At the moment, only [VeloxChem](https://veloxchem.org/) can be used as the QM driver.

## Overview

**PyMERK** processes multiple molecular conformers through progressive filtering stages, following, in spirit, [the CENSO paper](https://pubs.acs.org/doi/10.1021/acs.jpca.1c00971):

1. **Prescreening** - Fast single-point calculations to eliminate high-energy candidates
2. **Screening** - Refined DFT calculations with solvation corrections
3. **Optimization** - Full geometry optimization
4. **Refinement** - Boltzmann population filtering

## Documentation

Full documentation is available at: [https://pierre-24.github.io/pymerk/](https://pierre-24.github.io/pymerk/)

- [Installation & Usage Guide](https://pierre-24.github.io/pymerk/usage.html)
- [API Reference](https://pierre-24.github.io/pymerk/api.html)
- [Contributing Guide](https://pierre-24.github.io/pymerk/contributing.html)

## License

MIT License - See [LICENSE](LICENSE) file for details.

