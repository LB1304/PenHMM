## -----------------------------------------------------------------------------
## PenHMM: additional example
## -----------------------------------------------------------------------------
## This script illustrates a complete workflow:
##   1. simulation of longitudinal binary data;
##   2. cross-validation for selecting the penalty strength, with k fixed;
##   3. final estimation of the selected model.
##
## After installing PenHMM, this file can be accessed with:
## system.file("examples", "PenHMM_examples.R", package = "PenHMM")
## -----------------------------------------------------------------------------

library(PenHMM)

## Optional parallel computation for cross-validation. Set use_parallel <- TRUE
## and choose the number of workers if you want to run the CV step in parallel.
use_parallel <- TRUE
n_workers <- 5

## -----------------------------------------------------------------------------
## 1. Simulate data
## -----------------------------------------------------------------------------
set.seed(2801)

n_subjects <- 250
n_times <- 10
n_states <- 3

transition_matrix <- matrix(0.125, n_states, n_states) +
  diag(0.750 - 0.125, n_states)

support_points <- c(-20, -5, 5)
initial_prob <- rep(1 / n_states, n_states)
covariate_effects <- c(1, -1, 1, 1)

sim_data <- Draw.HMCovManifest(
  n = n_subjects,
  TT = n_times,
  al = support_points,
  be = covariate_effects,
  piv = initial_prob,
  Pi = transition_matrix
)

data <- sim_data$data

response_formula <- Y ~ X.1 + X.2 + X.3 + X.4
index_columns <- c("id", "time")

## -----------------------------------------------------------------------------
## 2. Select the penalty strength by cross-validation
## -----------------------------------------------------------------------------
## In this example, the number of hidden states is fixed at k = 3 and the
## cross-validation procedure selects the penalty strength. To jointly select k
## and the penalty strength, replace k_values = 3 with a vector such as 1:4.

k_candidates <- 3
gamma_candidates <- seq(0, 0.05, by = 0.01)

cv_result <- Select.CV(
  data = data,
  index = index_columns,
  responsesFormula = response_formula,
  k_values = k_candidates,
  gamma_values = gamma_candidates,
  n_folds = 5,
  start = 0,
  tol = 1e-4,
  maxit = 1e4,
  parallel = use_parallel,
  workers = n_workers,
  se_multiplier = 0.1,
  selection_rule = "smallest_k_largest_gamma",
  keep_fits = FALSE
)

cat("Selected number of hidden states:", cv_result$selected_k, "\n")
cat("Selected penalty strength:", cv_result$selected_gamma, "\n\n")

print(cv_result$cv_table)

## -----------------------------------------------------------------------------
## 3. Estimate the final selected model
## -----------------------------------------------------------------------------
final_fit <- Estimate.HMCovManifest(
  data = data,
  index = index_columns,
  responsesFormula = response_formula,
  k = cv_result$selected_k,
  start = 0,
  tol = 1e-4,
  maxit = 1e4,
  ga2 = cv_result$selected_gamma,
  out_se = TRUE,
  output = TRUE
)

## Main output components
final_fit$lk
final_fit$mu
final_fit$al
final_fit$be
final_fit$la
final_fit$PI
