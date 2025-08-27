MGAL k-dimensional estimation (R)

Contents
--------
R/
  - rMGAL.R                   : k-dim generator
  - mgal_em_marginal.R        : marginal EM for (alpha, mu, eta, diag(Sigma))
  - sample_w_function_k.R     : k-dim Gibbs/MH sampler for latent W
  - covariance_estimation_k.R : covariance matrix EM via latent W draws

examples/
  - example_2x2.R             : simulate -> marginal EM -> covariance EM

Usage
-----
In R:
  source('R/rMGAL.R')
  source('R/mgal_em_marginal.R')
  source('R/sample_w_function_k.R')
  source('R/covariance_estimation_k.R')
  source('examples/example_2x2.R')

Dependencies
------------
Packages: MASS, GIGrvg, Bessel
(Your original scripts also used rmutil, rootSolve, nleqslv as general deps.)

Notes
-----
- Scripts are kept as plain .R modules (no Rmd), matching your request.
- Example uses rMGAL() for data generation and covariance_estimation_k().
