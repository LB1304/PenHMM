\name{Select.CV}
\alias{Select.CV}
\title{Cross-Validation for Selecting the Number of Hidden States and Penalty Strength}
\description{
Performs cross-validation for hidden Markov models with covariates in the measurement model. The procedure can select both the number of hidden states and the penalty strength for the state-specific support points.
}
\usage{
Select.CV(data, index, responsesFormula, k_values, gamma_values,
  n_folds = 5, start = 0, tol = 1e-04, maxit = 1e+04,
  output = FALSE, out_se = FALSE, parallel = FALSE, workers = 1,
  future_seed = TRUE, se_multiplier = 0.1,
  selection_rule = c("smallest_k_largest_gamma",
    "largest_gamma_smallest_k", "max_mean"),
  shuffle_folds = FALSE, seed = NULL, keep_fits = TRUE)
}
\arguments{
  \item{data}{Data frame in long format.}
  \item{index}{Character vector of length two giving the names of the subject identifier and time variables in \code{data}.}
  \item{responsesFormula}{A formula specifying the response and covariates.}
  \item{k_values}{Integer vector containing the candidate numbers of hidden states.}
  \item{gamma_values}{Numeric vector containing the candidate penalty strengths. These values are passed to \code{ga2} in \code{Estimate.HMCovManifest}.}
  \item{n_folds}{Integer. Number of cross-validation folds.}
  \item{start}{Integer controlling the initialization strategy used in model estimation.}
  \item{tol}{Numeric. Convergence tolerance for the EM algorithm.}
  \item{maxit}{Integer. Maximum number of EM iterations.}
  \item{output}{Logical. Passed to \code{Estimate.HMCovManifest}.}
  \item{out_se}{Logical. Passed to \code{Estimate.HMCovManifest}.}
  \item{parallel}{Logical. If \code{TRUE}, models are fitted using \pkg{future.apply}.}
  \item{workers}{Integer. Number of workers used when \code{parallel = TRUE}.}
  \item{future_seed}{Logical or integer. Passed to \code{future.apply::future_lapply} for reproducible parallel random number generation.}
  \item{se_multiplier}{Numeric. Multiplier of the standard error used to define the selection threshold around the best cross-validated log-likelihood.}
  \item{selection_rule}{Character. Rule used to choose among candidate models within the threshold. \code{"smallest_k_largest_gamma"} prioritizes the smallest number of states and then the largest penalty; \code{"largest_gamma_smallest_k"} prioritizes the largest penalty and then the smallest number of states; \code{"max_mean"} selects the model with the maximum mean validation log-likelihood.}
  \item{shuffle_folds}{Logical. If \code{TRUE}, subject identifiers are randomly shuffled before assigning folds.}
  \item{seed}{Optional integer seed used when \code{shuffle_folds = TRUE}.}
  \item{keep_fits}{Logical. If \code{TRUE}, all fitted models are returned.}
}
\details{
Folds are assigned at the subject level, so all observations from the same subject are kept in the same fold. For each combination of \code{k_values}, \code{gamma_values}, and validation fold, the model is estimated on the training folds and the penalized log-likelihood is evaluated on the held-out fold.
}
\value{
A list with components:
\item{selected_k}{Selected number of hidden states.}
\item{selected_gamma}{Selected penalty strength.}
\item{selected_row}{The row of the cross-validation summary corresponding to the selected model.}
\item{cv_table}{Summary table with mean validation log-likelihood, standard error, threshold, and candidate indicator for each pair of \code{k} and \code{Gamma}.}
\item{fold_table}{Data frame mapping subject identifiers to cross-validation folds.}
\item{grid_results}{Fold-level validation log-likelihood values.}
\item{call}{The matched function call.}
\item{fits}{If \code{keep_fits = TRUE}, a list of fitted models for all grid and fold combinations.}
}
\examples{
\dontrun{
cv <- Select.CV(
  data = dat,
  index = c("id", "time"),
  responsesFormula = Y ~ X1 + X2,
  k_values = 1:4,
  gamma_values = seq(0, 0.05, by = 0.01),
  n_folds = 5,
  parallel = TRUE,
  workers = 4
)
cv$selected_k
cv$selected_gamma
}
}
\keyword{models}
\keyword{utilities}
