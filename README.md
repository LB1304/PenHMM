<h1 align="center">Penalized Hidden Markov Models with Covariates</h1>
<p align="center"> <span style="font-size: 14px;"><em><strong>Luca Brusa &middot; Fulvia Pennoni &middot; Francesco Bartolucci &middot; Romina Peruilh Bagolini</strong></em></span> </p>
<br>

`PenHMM` provides tools to simulate, estimate, and cross-validate discrete-time hidden Markov models with covariates in the measurement model using a penalized maximum likelihood approach. The penalty regularizes the dispersion of the state-specific support points, with the aim of addressing latent-state separation and preventing extreme estimates of the latent effects. The package also includes a cross-validation routine for selecting the number of hidden states and the penalty strength.

## Installation

To install the `PenHMM` package directly from GitHub:

```r
# install.packages("devtools")
require(devtools)
devtools::install_github("LB1304/PenHMM")
require(PenHMM)
```

To download the `.tar.gz` file for manual installation, use [this link](https://github.com/LB1304/PenHMM/archive/main.tar.gz).

The package requires `Formula`, `LMest` (version 3.2.7 or later), `MASS`, and `stats`. Parallel cross-validation additionally requires `future` and `future.apply`.

## Main functions

### `Draw.HMCovManifest()`

Simulates longitudinal binary data from a discrete-time hidden Markov model with covariates in the measurement model. The user can specify the latent-state support points, covariate effects, initial probabilities, transition probabilities, and whether the lagged response is included among the covariates.

### `Estimate.HMCovManifest()`

Estimates a hidden Markov model with covariates in the measurement model through maximum likelihood or penalized maximum likelihood. The penalized version regularizes the state-specific support points in order to reduce latent-state separation and improve the stability of the estimation procedure.

### `Select.CV()`

Performs cross-validation for model and tuning-parameter selection. The function evaluates candidate values for the number of hidden states and the penalty strength, and selects the preferred combination according to the chosen cross-validation rule. Both sequential and parallel computation are supported.

## Examples

Usage examples for each main function are available in the corresponding help pages. After installing the package, they can be opened in R with:

```r
?Draw.HMCovManifest
?Estimate.HMCovManifest
?Select.CV
```
Additional examples is available in [`inst/examples/PenHMM_examples.R`](inst/examples/PenHMM_examples.R). 
After installing the package, the example file can also be accessed directly from R using:

```r
example_file <- system.file("examples", "PenHMM_examples.R", package = "PenHMM")
source(example_file)
```
