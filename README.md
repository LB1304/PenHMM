# PenHMM

`PenHMM` provides functions to simulate, estimate, and cross-validate hidden Markov models with covariates in the measurement model using a penalized maximum likelihood approach.

## Main functions

- `Draw.HMCovManifest()` simulates longitudinal binary data from a hidden Markov model with covariates.
- `Estimate.HMCovManifest()` estimates the model using penalized maximum likelihood.
- `Select.CV()` performs cross-validation to select the number of hidden states and the penalty strength.

## Installation from source

From the directory containing the package source:

```r
install.packages("PenHMM", repos = NULL, type = "source")
```

The package requires `Formula`, `LMest`, and `MASS`. Parallel cross-validation additionally requires `future` and `future.apply`.
