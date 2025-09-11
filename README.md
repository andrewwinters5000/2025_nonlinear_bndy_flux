# Numerical flux functions that impose nonlinearly stable open boundary conditions

[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/TODO/zenodo.TODO.svg)](https://doi.org/TODO/zenodo.TODO)


This repository contains information and code to reproduce the results presented in the
article
```bibtex
@online{winters2025numerical,
  title={Numerical flux functions that impose nonlinearly stable open boundary conditions},
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
         "{N}umerical flux functions that impose nonlinearly stable open boundary conditions"},
  author={Winters, Andrew R and Kopriva, David A and Nordström, Jan},
  year={2025},
  howpublished={\url{https://github.com/andrewwinters5000/2025_nonlinear_bndy_flux}},
  doi={TODO/zenodo.TODO}
}
```

## Abstract

We describe how provably stable nonlinear simultaneous approximation terms (SATs) that weakly impose boundary conditions can be reinterpreted in terms of a two-point numerical flux function.
Initially, this analysis is done for the one-dimensional Burgers equation, as its simplicity illustrates the fundamental ideas.
We then extend it into two-dimensional curvilinear coordinates with the shallow water equations.
The resulting numerical flux function inherits the provable stability bound of the SATs for nonlinear problems.
Additionally, this formulation can be used for any approximation strategy that uses a numerical flux to weakly impose boundary conditions, e.g. finite volume methods.
We provide numerical studies to demonstrate the performance of the novel boundary flux functions in the context of high-order discontinuous Galerkin spectral element methods.


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
