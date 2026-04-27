Compute.LogLikSta <- function (tau, u, V, G2, outl = TRUE, ga1, ga2, alpha) {
  k = length(u)
  PI = exp(G2 %*% tau)
  PI = t(matrix(PI, k, k))
  PI = diag(1/rowSums(PI)) %*% PI
  PI1 = PI
  for (i in 1:10000) PI1 = PI1 %*% PI
  la = colMeans(PI1)
  la = la/sum(la)
  
  if (k > 1) {
    penalty <- ComputePenalization(alpha, ga1, ga2, k)
  } else {
    penalty <- 0
  }
  lk <- u %*% log(la) + sum(V * log(PI)) - penalty
  
  flk = -lk
  if (outl) 
    out = list(flk = flk, la = la, PI = PI)
  else flk
}

