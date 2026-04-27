ComputePenalization <- function(alpha, ga1, ga2, k) {
  E <- rbind(rep(0, k - 1), diag(k - 1))
  
  # Penalizzazione nello spazio completo
  one <- rep(1, k)
  I_k <- diag(k)
  M <- I_k - (1 / k) * (one %*% t(one))
  A <- ga1 * I_k + ga2 * M
  
  # Penalizzazione nello spazio ridotto
  A_red <- t(E) %*% A %*% E
  
  as.numeric(t(alpha) %*% A_red %*% alpha)
}


ComputePenalization.Gradient <- function(alpha, ga1, ga2, k) {
  E <- rbind(rep(0, k - 1), diag(k - 1))
  
  one <- rep(1, k)
  I_k <- diag(k)
  M <- I_k - (1 / k) * (one %*% t(one))
  A <- ga1 * I_k + ga2 * M
  
  A_red <- t(E) %*% A %*% E
  
  2 * A_red %*% alpha
}


ComputePenalization.Hessian <- function(ga1, ga2, k) {
  E <- rbind(rep(0, k - 1), diag(k - 1))
  
  one <- rep(1, k)
  I_k <- diag(k)
  M <- I_k - (1 / k) * (one %*% t(one))
  A <- ga1 * I_k + ga2 * M
  
  2 * t(E) %*% A %*% E
}

