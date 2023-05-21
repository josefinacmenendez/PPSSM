# A point process state space model for estimating the instantaneous rate of the process given a time-series of event times.

This repository contains MATLAB routines for estimating spike rate functions via state space smoothing.

Algorithmic implementations stem from the following manuscript:

Smith, A.C., Scalon, J.D., Wirth, S., Yanike, M., Suzuki, W.A. and Brown, E.N., 2010. State-space algorithms for estimating spike rate functions. Computational Intelligence and Neuroscience, 2010.

The implementation in this repository improves the numerical stability of the Newton step in Smith et al., 2010, by using a damped newton algorithm for estimating the latent state.

If you use this code, please cite the repository:

@software{PPSSM,
  author = {Josefina Correa-Menendez},
  doi = {10.5281/zendo.7955060},
  month = {5},
  title = {{A numerically stable state-space model for spike rate estimation}},
  url = {https://github.com/josefinacmenendez/PPSSM},
  version = {1.0},
  year = {2023}
}




[![DOI](https://zenodo.org/badge/459431601.svg)](https://zenodo.org/badge/latestdoi/459431601)

