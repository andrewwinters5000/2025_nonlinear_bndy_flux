# Numerical flux functions that give provable bounds for nonlinear initial boundary value problems with open boundaries

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/TODO/zenodo.TODO.svg)](https://doi.org/TODO/zenodo.TODO)


This repository contains information and code to reproduce the results presented in the
article
```bibtex
@online{winters2025numerical,
  title={Numerical flux functions that give provable bounds for nonlinear initial boundary value problems with open boundaries},
  author={Winters, Andrew R and Kopriva, David A and Nordström, Jan},
  year={2025},
  month={09},
  eprint={TODO},
  eprinttype={arxiv},
  eprintclass={math.NA}
}
```

If you find these results useful, please cite the article mentioned above. If you
use the implementations provided here, please **also** cite this repository as
```bibtex
@misc{winters2025numericalRepro,
  title={Reproducibility repository for
         "{N}umerical flux functions that give provable bounds for nonlinear initial boundary value problems with open boundaries"},
  author={Winters, Andrew R and Kopriva, David A and Nordström, Jan},
  year={2025},
  howpublished={\url{https://github.com/andrewwinters5000/2025_nonlinear_bndy_flux}},
  doi={TODO/zenodo.TODO}
}
```

## Abstract

We present a strategy to interpret nonlinear characteristic-type penalty terms as numerical flux functions that provide provable bounds for solutions of nonlinear hyperbolic initial boundary value problems with open boundaries.
This interpretation is rooted in the nonlinear energy method, which has recently enabled a systematic design of characteristic-based penalty terms for the weak imposition of boundary conditions.
The new boundary fluxes are directly compatible with high-order accurate split form discontinuous Galerkin spectral element methods and guarantee that
the solution is energy stable and bounded
solely by external data.
We derive inflow-outflow boundary fluxes specifically for the Burgers equation and the two-dimensional shallow water equations.
Numerical experiments demonstrate that the new nonlinear fluxes do not fail in situations where standard boundary treatments based on linear analysis do.

## Numerical experiments

To reproduce the numerical experiments presented in this article, you need
to install [Julia](https://julialang.org/). The numerical experiments presented
in this article were performed using Julia v1.11.3.

First, you need to download this repository, e.g., by cloning it with `git`
or by downloading an archive via the GitHub interface. Then, you need to start
Julia in the `code` directory of this repository and follow the instructions
described in the `README.md` file therein.


## Authors

- Andrew R. Winters (Linköping University, Sweden)
- David A. Kopriva (Florida State University, Florida, USA and San Diego State University, California, USA)
- Jan Nordström (Linköping University, Sweden and University of Johannesburg, South Africa)


## License

The code in this repository is published under the MIT license, see the
`LICENSE` file.


## Disclaimer

Everything is provided as is and without warranty. Use at your own risk!
