Draw.HMCovManifest <- function(n, TT, al, be, piv, Pi, Lagged = FALSE) {
  
  k = length(al)
  nc = ifelse(Lagged, length(be) - 1, length(be))
  
  #---- Covariates ----
  X = array(rnorm(n*TT*nc), c(n, TT, nc))
  dimnames(X) = list("unit" = 1:n, "time" = 1:TT, "covariate" = 1:nc)
  
  #---- Latent variable ----
  U = matrix(0, n, TT)
  dimnames(U) = list("unit" = 1:n, "time" = 1:TT)
  for(i in 1:n){
    U[i,1] = which(rmultinom(1,1,piv)==1)
    for(t in 2:TT){
      U[i,t] = which(rmultinom(1,1,Pi[,U[i,t-1]])==1)
    }
  }
  
  #---- Response variable and Lagged response ----
  P = matrix(0, n, TT)
  Y = matrix(0, n, TT)
  if (Lagged) {
    LS = matrix(0, n, TT)
    LS[, 1] = 1*(runif(1) < 0.05)
  }
  for (t in 1:TT) for (i in 1:n) {
    if (Lagged) {
      if (t > 1) LS[i, t] = Y[i, t-1]
      p = exp(al[U[i, t]] + c(X[i, t, ], LS[i, t]) %*% be)
    } else {
      p = exp(al[U[i, t]] + c(X[i, t, ]) %*% be)
    }
    p = p/(1+p)
    P[i, t] = p
    Y[i, t] = 1*(runif(1) < p)
  }
  dimnames(Y) = list("unit" = 1:n, "time" = 1:TT)
  
  #---- Output dataset ----
  XX = matrix(aperm(X, c(2, 1, 3)), n*TT, nc)
  if (Lagged) {
    data = data.frame(id = rep(1:n, each = TT), time = rep(1:TT, n), 
                      X = XX, Y = c(t(Y)), LY = c(t(LS)))
  } else {
    data = data.frame(id = rep(1:n, each = TT), time = rep(1:TT, n), 
                      X = XX, Y = c(t(Y)))
  }
  
  return(list(data = data, U = U))
}



