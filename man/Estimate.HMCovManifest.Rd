\name{Estimate.HMCovManifest}
\alias{Estimate.HMCovManifest}
\title{Estimate a Penalized Hidden Markov Model with Covariates in the Measurement Model}
\description{
Estimates a discrete-time hidden Markov model with covariates in the measurement model using maximum likelihood or penalized maximum likelihood. The penalty regularizes the state-specific support points to mitigate latent-state separation.
}
\usage{
Estimate.HMCovManifest(data, responsesFormula, index, k, start = 0,
  ga1 = 0, ga2 = 0, tol = 1e-08, maxit = 1e+04,
  out_se = FALSE, output = FALSE)
}
\arguments{
  \item{data}{Data frame in long format.}
  \item{responsesFormula}{A formula, preferably created using the \pkg{Formula} syntax, specifying the response on the left-hand side and covariates on the right-hand side.}
  \item{index}{Character vector of length two giving the names of the subject identifier and time variables in \code{data}.}
  \item{k}{Integer. Number of hidden states.}
  \item{start}{Integer controlling the initialization strategy. The default \code{0} uses deterministic initialization; \code{1} uses random initialization.}
  \item{ga1}{Numeric. Ridge-type penalty component. Default is \code{0}.}
  \item{ga2}{Numeric. Penalty strength for the dispersion of the state-specific support points. Default is \code{0}, corresponding to the unpenalized estimator.}
  \item{tol}{Numeric. Convergence tolerance for the EM algorithm.}
  \item{maxit}{Integer. Maximum number of EM iterations.}
  \item{out_se}{Logical. If \code{TRUE}, standard errors for the regression coefficients are computed.}
  \item{output}{Logical. If \code{TRUE}, posterior probabilities and predicted latent effects are returned.}
}
\details{
The model assumes a first-order hidden Markov chain and a generalized logit measurement model with covariates. When \code{ga2 > 0}, the likelihood is penalized by a quadratic term acting on the dispersion of the latent support points. Estimation is performed by an EM algorithm with a modified M-step.
}
\value{
A list containing model estimates and fit measures. Main components include:
\item{mu}{Estimated cut-point or intercept component.}
\item{al}{Estimated state-specific support points, centered according to the fitted latent distribution.}
\item{be}{Estimated regression coefficients.}
\item{la}{Estimated initial probabilities.}
\item{PI}{Estimated transition probability matrix.}
\item{par}{Internal parameter vector.}
\item{lk}{Final penalized log-likelihood value.}
\item{lkv}{Trace of the log-likelihood values across iterations.}
\item{np}{Number of free parameters used for information criteria.}
\item{aic}{Akaike information criterion.}
\item{bic}{Bayesian information criterion.}
\item{n}{Number of subjects.}
\item{k}{Number of hidden states.}
\item{TT}{Number of time occasions.}
When \code{out_se = TRUE}, additional components \code{sebe} and \code{J1} are returned. When \code{output = TRUE}, posterior probabilities and predicted quantities are also returned.
}
\examples{
\dontrun{
fit <- Estimate.HMCovManifest(
  data = dat,
  responsesFormula = Y ~ X1 + X2,
  index = c("id", "time"),
  k = 3,
  ga2 = 0.01
)
fit$al
fit$be
}
}
\keyword{models}
\keyword{multivariate}
