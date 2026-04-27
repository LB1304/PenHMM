\name{Draw.HMCovManifest}
\alias{Draw.HMCovManifest}
\title{Simulate Binary Longitudinal Data from a Hidden Markov Model with Covariates}
\description{
Simulates a longitudinal binary response from a discrete-time hidden Markov model with covariates in the measurement model. Optionally, the lagged response is included among the covariates.
}
\usage{
Draw.HMCovManifest(n, TT, al, be, piv, Pi, Lagged = FALSE)
}
\arguments{
  \item{n}{Integer. Number of subjects.}
  \item{TT}{Integer. Number of time occasions.}
  \item{al}{Numeric vector of state-specific support points. Its length determines the number of hidden states.}
  \item{be}{Numeric vector of regression coefficients for the covariates. If \code{Lagged = TRUE}, the last coefficient is used for the lagged response.}
  \item{piv}{Numeric vector of initial probabilities for the hidden states.}
  \item{Pi}{Transition probability matrix for the hidden Markov chain. Columns correspond to the previous state and rows to the next state, consistently with the simulation routine.}
  \item{Lagged}{Logical. If \code{TRUE}, the lagged response is included as an additional covariate in the measurement model.}
}
\details{
For each subject, a hidden Markov chain is generated using the initial probabilities \code{piv} and transition matrix \code{Pi}. Conditional on the latent state and covariates, the response is generated from a Bernoulli distribution with logit depending on the corresponding support point and regression coefficients.
}
\value{
A list with components:
\item{data}{A data frame in long format containing subject identifier, time occasion, simulated covariates, response variable, and, when requested, the lagged response.}
\item{U}{A matrix containing the simulated hidden state sequence for each subject and time occasion.}
}
\examples{
set.seed(1)
Pi <- matrix(c(0.90, 0.10, 0.10, 0.90), nrow = 2)
piv <- c(0.5, 0.5)
out <- Draw.HMCovManifest(
  n = 20,
  TT = 5,
  al = c(-2, 2),
  be = c(0.5, -0.5),
  piv = piv,
  Pi = Pi
)
head(out$data)
}
\keyword{models}
\keyword{datagen}
