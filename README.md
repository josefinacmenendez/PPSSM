# SpikeRates

MATLAB routines for estimating spike rate functions via state space smoothing.

Algorithmic implementations stem from the following manuscript:

Smith, A.C., Scalon, J.D., Wirth, S., Yanike, M., Suzuki, W.A. and Brown, E.N., 2010. State-space algorithms for estimating spike rate functions. Computational Intelligence and Neuroscience, 2010.

The implementation in this repository improves the numerical stability of the Newton step in Smith et al., 2010, by using a damped newton algorithm for estimating the latent state.
