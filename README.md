# `pymerk` -- Python Molecular Ensembles RanKing

A CENSO-like approach to rank ensembles of molecular systems using quantum chemistry calculations.

## Overview

**PyMERK** processes multiple molecular conformers through progressive filtering stages:

1. **Prescreening** - Fast single-point calculations to eliminate high-energy candidates
2. **Screening** - Refined DFT calculations with solvation corrections
3. **Optimization** - Full geometry optimization
4. **Refinement** - Boltzmann population filtering

## Quick Start

```bash
# Install with development dependencies
python -m venv venv
source venv/bin/activate
make install-dev

# Run PyMERK
pymerk_run conformers.xyz config.toml ./results
```

## Documentation

Full documentation is available at: **[https://pierre-24.github.io/pymerk/](https://pierre-24.github.io/pymerk/)**

- [Installation & Usage Guide](https://pierre-24.github.io/pymerk/usage.html)
- [API Reference](https://pierre-24.github.io/pymerk/api.html)
- [Contributing Guide](https://pierre-24.github.io/pymerk/contributing.html)

## Requirements

- Python 3.11 or higher
- xTB and/or VeloxChem (optional but recommended)

## License

MIT License - See LICENSE file for details.

