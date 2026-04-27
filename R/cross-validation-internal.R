Compute.Test.LLK <- function(data, index, responsesFormula, par, la, PI, k, ga1 = 0, ga2 = 0) {
  
  ## S, X, yv, lev
  prepared_data <- PrepareData(
    data = data,
    responsesFormula = responsesFormula,
    index = index
  )
  S <- as.matrix(prepared_data$Y[, , 1])
  X <- prepared_data$X
  yv <- prepared_data$freq
  
  lev = max(S) + 1
  nt = prod(lev)
  ns = nrow(S)
  TT = ncol(S)
  n = sum(yv)
  Y0 = S + 1
  S = array(0, c(nt, ns, TT))
  for (i in 1:ns) for (t in 1:TT) {
    ind = Y0[i, t]
    S[ind, i, t] = 1
  }
  if (is.matrix(X)) {
    X = array(X, c(ns, TT, 1))
  }
  nc = dim(X)[3]
  ne = lev - 1
  XX = X
  X = array(0, c(ne, nc, ns, TT))
  for (i in 1:ns) for (t in 1:TT) {
    if (lev == 2) X[, , i, t] = XX[i, t, ]
    else X[, , i, t] = rep(1, ne) %o% XX[i, t, ]
  }
  
  ## Xd, indn
  X1 = matrix(X, ne * nc, ns * TT)
  out1 = t(unique(t(X1)))
  nd = ncol(out1)
  indn = rep(0, ns * TT)
  tmp = ne * nc
  for (jd in 1:nd) {
    ind = which(colSums(X1 == out1[, jd]) == tmp)
    indn[ind] = jd
  }
  indn = matrix(indn, ns, TT)
  Xd = array(out1, c(ne, nc, nd))
  
  ## I, one
  I = diag(ne)
  one = matrix(1, ne, 1)
  
  ## lm, Lm
  lm = c(1, rep(0, lev - 1))
  Lm = rbind(rep(0, lev - 1), diag(lev - 1)) - rbind(diag(lev - 1), rep(0, lev - 1))
  
  ## par0
  par0 = par[1:(lev - 2 + k)]
  
  ## Eta01
  Eta01 = LMest:::prod_array(Xd, par[(lev + k - 1):length(par)])
  
  ## Pio
  Pio = array(0, c(ns, k, TT))
  for (c in 1:k) {
    u = matrix(0, 1, k)
    u[c] = 1
    u = u[-1]
    D0 = cbind(I, t(as.matrix(u)) %x% one)
    agg = as.vector(D0 %*% par0)
    Eta1 = Eta01 + agg %o% rep(1, nd)
    Qv1 = LMest:::expit(Eta1)
    Qv1 = pmin(pmax(Qv1, 1e-100), 1 - 1e-100)
    Pv1 = lm %o% rep(1, nd) + Lm %*% Qv1
    Pv1 = pmin(pmax(Pv1, 1e-100), 1 - 1e-100)
    for (t in 1:TT) {
      Pio[, c, t] = colSums(S[, , t] * Pv1[, indn[, t]])
    }
  }
  
  ## Q
  Q = LMest:::rec1(Pio, la, PI)
  
  ## pim
  if (k == 1) {
    pim = Q[, , TT]
  } else {
    pim = rowSums(Q[, , TT])
  }
  pim <- pmin(pmax(pim, 10^-100), 1 - 10^-100)
  
  if (k > 1) {
    alpha <- par[(ne + 1):(ne + k - 1)]
    penalty <- ComputePenalization(alpha, ga1, ga2, k)
  } else {
    penalty <- 0
  }
  lk <- sum(yv * log(pim)) - penalty
  
  return(lk)
}






