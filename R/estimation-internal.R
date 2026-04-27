EstimateModel <- function(S, X, yv = rep(1, nrow(S)), k, tol = 1e-8, maxit = 1e4, 
                          mu = NULL, al = NULL, be = NULL, la = NULL, PI = NULL, 
                          start = 0, output = FALSE, out_se = FALSE, ga1 = 0, ga2 = 0) {
  
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
  if (is.matrix(X)) 
    X = array(X, c(ns, TT, 1))
  nc = dim(X)[3]
  ne = lev - 1
  XX = X
  X = array(0, c(ne, nc, ns, TT))
  for (i in 1:ns) for (t in 1:TT) {
    if (lev == 2) 
      X[, , i, t] = XX[i, t, ]
    else X[, , i, t] = rep(1, ne) %o% XX[i, t, ]
  }
  out = LMest:::marg_param(lev, "g")
  Cm = out$C
  Mm = out$M
  Gm = cbind(-rep(1, lev - 1), diag(lev - 1))
  Hm = rbind(rep(0, lev - 1), diag(lev - 1))
  GHt = t(Gm) %*% t(Hm)
  lm = c(1, rep(0, lev - 1))
  Lm = rbind(rep(0, lev - 1), diag(lev - 1)) - rbind(diag(lev - 1), rep(0, lev - 1))
  G2 = NULL
  H2 = NULL
  IPI = NULL
  if (k > 1) {
    for (c in 1:k) {
      G2c = diag(k)[, -c]
      H2c = diag(k)[-c, ]
      if (k == 2) 
        H2c[c] = -1
      else H2c[, c] = -1
      if (is.null(G2)) 
        G2 = G2c
      else if (k == 2) 
        G2 = LMest:::blkdiag(matrix(G2, ncol = 1), matrix(G2c, ncol = 1))
      else G2 = LMest:::blkdiag(G2, G2c)
      if (is.null(H2)) 
        H2 = H2c
      else if (k == 2) 
        H2 = LMest:::blkdiag(matrix(H2, nrow = 1), matrix(H2c, nrow = 1))
      else H2 = LMest:::blkdiag(H2, H2c)
      IPI = c(IPI, c + seq(0, k * (k - 1), k))
    }
  }
  mu_inp = mu
  if (is.null(mu)) {
    Pim = apply(S, c(1, 2), sum) + 0.05 * TT
    Eta = Cm %*% log(Mm %*% Pim)
    Eta = Eta %x% matrix(1, 1, TT)
    eta = as.vector(Eta)
    Z = matrix(aperm(X, c(1, 4, 3, 2)), ns * ne * TT, dim(X)[2])
    Z = cbind(matrix(1, ns * TT, 1) %x% diag(ne), Z)
    par = MASS:::ginv(t(Z) %*% Z) %*% t(Z) %*% eta
    mu = par[1:ne]
    par = par[-(1:ne)]
    be = par
    if (k == 1) 
      al = NULL
    else {
      if (start == 1) 
        al = rnorm(k) * k
      else al = seq(-k, k, 2 * k/(k - 1))
      mu = mu + al[1]
      al = al[-1] - al[1]
    }
    if (k == 1) 
      PI = 1
    else {
      PI = matrix(1, k, k) + 9 * diag(k)
      PI = diag(1/rowSums(PI)) %*% PI
    }
    if (start == 1) {
      la = matrix(runif(k), k, 1)
      la = la/sum(la)
    } else {
      la = matrix(1, k, 1)/k
    }
  }
  if (start == 2) {
    if (is.null(mu_inp)) 
      stop("initial value of the cut-points (mu) must be given in input")
    mu = mu_inp
    if (is.null(be)) 
      stop("initial value of the regression parameters (be) must be given in input")
    be = be
    if (is.null(al)) 
      stop("initial value of the support points (al) must be given in input")
    mu = mu + al[1]
    al = al[-1] - al[1]
    if (is.null(la)) 
      stop("initial value of the initial probabilities (la) must be given in input")
    la = la
    if (is.null(PI)) 
      stop("initial value of the transition probabilities (PI) must be given in input")
    PI = PI
  }
  par = c(mu, al, be)
  if (k == 1) 
    tau = NULL
  else {
    tau = H2 %*% log(PI[IPI])
  }
  la = as.vector(la)
  PI = PI %x% matrix(1, 1, 1)
  X1 = matrix(X, ne * nc, ns * TT)
  out1 = t(unique(t(X1)))
  nd = ncol(out1)
  indn = rep(0, ns * TT)
  INDN = vector("list", nd)
  tmp = ne * nc
  
  for (jd in 1:nd) {
    ind = which(colSums(X1 == out1[, jd]) == tmp)
    indn[ind] = jd
    INDN[[jd]]$ind = ind
  }
  indn = matrix(indn, ns, TT)
  Xd = array(out1, c(ne, nc, nd))
  LLm1 = array(t(Lm), c(ncol(Lm), nrow(Lm), nd))
  
  I = diag(ne)
  one = matrix(1, ne, 1)
  Pio = array(0, c(ns, k, TT))
  par0 = par[1:(lev - 2 + k)]
  Eta01 = LMest:::prod_array(Xd, par[(lev + k - 1):length(par)])
  for (c in 1:k) {
    u = matrix(0, 1, k)
    u[c] = 1
    u = u[-1]
    D0 = cbind(I, matrix(u, nrow = 1) %x% one)
    agg = D0 %*% par0
    Eta1 = Eta01 + agg %*% rep(1, nd)
    Qv1 = LMest:::expit(Eta1)
    Qv1 = pmin(pmax(Qv1, 10^-100), 1 - 10^-100)
    Pv1 = lm %o% rep(1, nd) + Lm %*% Qv1
    Pv1 = pmin(pmax(Pv1, 10^-100), 1 - 10^-100)
    for (t in 1:TT) Pio[, c, t] = colSums(S[, , t] * Pv1[, indn[, t]])
  }
  Q = LMest:::rec1(Pio, la, PI)
  if (k == 1) pim = Q[, , TT]
  else pim = rowSums(Q[, , TT])
  
  if (k > 1) {
    alpha <- par[(ne + 1):(ne + k - 1)]
    penalty <- ComputePenalization(alpha, ga1, ga2, k)
  } else {
    penalty <- 0
  }
  lk <- sum(yv * log(pim)) - penalty
  
  it = 0
  lko = lk - 10^10
  lkv = NULL
  dis = 0
  while ((lk - lko)/abs(lk) > tol & it < maxit) {
    it = it + 1
    lko = lk
    paro = par
    tauo = tau
    out = LMest:::rec3(Q, yv, PI, Pio, pim)
    U = out$U
    V = out$V
    if (k > 1) {
      alpha <- par[(ne + 1):(ne + k - 1)]
      out = optim(tau, Compute.LogLikSta, gr = NULL, rowSums(U[, , 1]), V, G2, outl = FALSE, ga1, ga2, alpha, method = "BFGS")
      tau = out$par
      out = Compute.LogLikSta(tau, rowSums(U[, , 1]), V, G2, outl = TRUE, ga1, ga2, alpha)
      
      la = out$la
      PI = out$PI
    }
    
    la = as.vector(la)
    PI = PI %x% matrix(1, 1, 1)
    U = aperm(U, c(2, 1, 3))
    s = 0
    FF = 0
    for (c in 1:k) {
      u = matrix(0, 1, k)
      u[c] = 1
      u = u[-1]
      D0 = cbind(I, t(as.matrix(u)) %x% one)
      agg = as.vector(D0 %*% par0)
      Eta1 = Eta01 + agg %o% rep(1, nd)
      Qv1 = LMest:::expit(Eta1)
      Qv1 = pmin(pmax(Qv1, 10^-100), 1 - 10^-100)
      Pit1 = lm %o% rep(1, nd) + Lm %*% Qv1
      Pit1 = pmin(pmax(Pit1, 10^-100), 1 - 10^-100)
      QQv1 = Qv1 * (1 - Qv1)
      DPv1 = 1/Pit1
      RRtc1 = array(0, c(ne, lev, nd))
      for (j1 in 1:ne) for (j2 in 1:lev) RRtc1[j1, j2, ] = QQv1[j1, ] * DPv1[j2, ]
      RRtc1 = RRtc1 * LLm1
      XXRi1 = array(0, c(dim(D0)[1], dim(D0)[2] + dim(Xd)[2], nd))
      for (h2 in 1:nd) {
        if (lev == 2) 
          XXRi1[, , h2] = c(D0, Xd[, , h2])
        else XXRi1[, , h2] = cbind(D0, Xd[, , h2])
      }
      XXRi1 = aperm(XXRi1, c(2, 1, 3))
      pc = U[, c, ]
      pc = as.vector(pc)
      nt = dim(S)[1]
      YGP = matrix(S, nt, ns * TT) - Pit1[, as.vector(indn)]
      Om = array(0, c(lev, lev, nd))
      for (r1 in 1:lev) for (r2 in 1:lev) {
        if (r2 == r1) {
          Om[r1, r2, ] = Pit1[r1, ] - Pit1[r1, ] * Pit1[r2, ]
        } else {
          Om[r1, r2, ] = -Pit1[r1, ] * Pit1[r2, ]
        }
      }
      for (jd in 1:nd) {
        ind = INDN[[jd]]$ind
        pci = pc[ind]
        if (lev == 2) 
          XRi = (XXRi1[, , jd] %o% RRtc1[, , jd]) %*% GHt
        else XRi = (XXRi1[, , jd] %*% RRtc1[, , jd]) %*% GHt
        if (length(ind) == 1) {
          s = s + XRi %*% (YGP[, ind] * pci)
        } else {
          s = s + XRi %*% (YGP[, ind] %*% pci)
        }
        FF = FF + sum(pci) * (XRi %*% Om[, , jd]) %*% t(XRi)
      }
    }
    
    if (k > 1) {
      alpha <- par[(ne + 1):(ne + k - 1)]
      penalty_s <- ComputePenalization.Gradient(alpha, ga1, ga2, k)
      penalty_FF <- ComputePenalization.Hessian(ga1, ga2, k)
    } else {
      penalty_s <- 0
      penalty_FF <- 0
    }
    s[(ne + 1):(ne + k - 1)] <- s[(ne + 1):(ne + k - 1)] - penalty_s
    FF[(ne + 1):(ne + k - 1), (ne + 1):(ne + k - 1)] <- FF[(ne + 1):(ne + k - 1), (ne + 1):(ne + k - 1)] - penalty_FF
    
    dpar = MASS:::ginv(FF) %*% s
    mdpar = max(abs(dpar))
    if (mdpar > 1) 
      dpar = dpar/mdpar
    par = par + dpar
    par0 = par[1:(lev - 2 + k)]
    Eta01 = LMest:::prod_array(Xd, par[(lev + k - 1):length(par)])
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
    Q = LMest:::rec1(Pio, la, PI)
    if (k == 1) 
      pim = Q[, , TT]
    else pim = rowSums(Q[, , TT])
    pim <- pmin(pmax(pim, 10^-100), 1 - 10^-100)
    
    if (k > 1) {
      alpha <- par[(ne + 1):(ne + k - 1)]
      penalty <- ComputePenalization(alpha, ga1, ga2, k)
    } else {
      penalty <- 0
    }
    lk <- sum(yv * log(pim)) - penalty
    
    dis = max(abs(c(par - paro, tau - tauo)))
    lkv = c(lkv, lk)
  }
  par1 = NULL
  if (k > 1) 
    par1 = tau
  par1 = c(par1, par)
  
  mu = par[1:ne]
  al = 0
  if (k > 1) {
    al = c(al, par[(ne + 1):(ne + k - 1)])
  }
  
  mu = mu + as.vector(al %*% la)
  al = al - as.vector(al %*% la)
  be = par[(ne + k):length(par)]
  np = k * (k - 1)
  np = np + (ne + (k - 1) + nc)
  aic = -2 * lk + 2 * (np)
  bic = -2 * lk + log(n) * (np)
  
  if (output) {
    out = Compute.LogLikObs(par1, S, Xd, yv, indn, lev, k,  G2, IPI, outp = TRUE)
    U = out$U
    PRED0 = array(0, c(ns, k, TT))
    PRED1 = matrix(0, ns, TT)
    for (t in 1:TT) {
      PRED0[, , t] = U[, , t]
      PRED1[, t] = U[, , t] %*% al
    }
    if (any(yv != 1)) {
      PRED0 = PRED0/yv
      PRED1 = PRED1/yv
    }
  }
  
  if (out_se) {
    out = Compute.LogLikObs(par1, S, Xd, yv, indn, lev, k, G2, IPI, outp = TRUE)
    nx = length(par1)
    d0 = out$s
    ny = length(d0)
    D = matrix(0, nx, ny)
    for (i in 1:nx) {
      o = matrix(0, nx, 1)
      o[i] = 10^-6
      out = Compute.LogLikObs(par1 + o, S, Xd, yv, indn, lev, k, G2, IPI, outp = TRUE)
      d1 = out$s
      d = (d1 - d0)/10^-6
      D[i, ] = t(d)
    }
    J1 = D
    J1 = -(J1 + t(J1))/2
    
    if (rcond(J1) < 10^-15) 
      print(c("rcond of information = ", rcond(J1)))
    se1 = diag(MASS:::ginv(J1))
    if (k > 1) {
      se1 = se1[-(1:(k * (k - 1)))]
    }
    sebe = sqrt(se1[(ne + k):length(se1)])
  }
  out = list(mu = mu, al = al, be = be, la = la, PI = PI, par = par, 
             lk = lk, lkv = lkv, np = np, aic = aic, bic = bic, 
             n = ns, k = k, TT = TT)
  if (out_se) {
    out$sebe = sebe
    out$J1 = J1
  }
  if (output) {
    if (k > 1) {
      Pmarg <- as.matrix(la)
      for (t in 2:TT) Pmarg = cbind(Pmarg, t(PI) %*% Pmarg[, t - 1])
    }
    else Pmarg = NULL
    out = c(out, list(V = PRED0, PRED1 = PRED1, S = S, yv = yv, Pmarg = Pmarg))
  }
  return(out)
}









Compute.LogLikObs <- function (par, Y, Xd, yv, indn, lev, k, G2, IPI, outp = FALSE) {
  
  n = dim(Y)[2]
  TT = dim(Y)[3]
  nd = dim(Xd)[3]
  np = lev - 1
  Gm = cbind(-rep(1, lev - 1, 1), diag(lev - 1))
  Hm = rbind(rep(0, lev - 1), diag(lev - 1))
  GHt = t(Gm) %*% t(Hm)
  lm = c(1, rep(0, lev - 1))
  Lm = rbind(rep(0, lev - 1), diag(lev - 1)) - rbind(diag(lev - 1), rep(0, lev - 1))
  INDN = vector("list", nd)
  for (jd in 1:nd) INDN[[jd]]$ind = which(as.vector(indn) == jd)
  out = LMest:::trans_par(par, lev, k, 0, G2, IPI, 0)
  la = out$la
  PI = out$PI
  par = out$par
  tau = out$tau
  las = as.vector(la %x% 1)
  PIs = PI
  
  LLm1 = array(0, c(dim(Lm)[2], dim(Lm)[1], nd))
  for (re in 1:nd) LLm1[, , re] = t(Lm)
  I = diag(np)
  one = matrix(1, np, 1)
  Pio = array(0, c(n, k, TT))
  
  par0 = par[1:(lev - 2 + k)]
  Eta01 = LMest:::prod_array(Xd, par[(lev + k - 1):length(par)])
  
  for (c in 1:k) {
    u = matrix(0, 1, k)
    u[c] = 1
    u = u[-1]
    D0 = cbind(I, matrix(u, nrow = 1) %x% one)
    agg = D0 %*% par0
    Eta1 = Eta01 + agg %*% matrix(1, 1, nd)
    Qv1 = LMest:::expit(Eta1)
    Qv1 = pmin(pmax(Qv1, 10^-100), 1 - 10^-100)
    Pv1 = lm %*% matrix(1, 1, nd) + Lm %*% Qv1
    Pv1 = pmin(pmax(Pv1, 10^-100), 1 - 10^-100)
    for (t in 1:TT) Pio[, c, t] = colSums(Y[, , t] * Pv1[, indn[, t]])
  }
  if (k == 1) 
    PIs = as.matrix(PIs)
  Q = LMest:::rec1(Pio, las, PIs)
  if (k == 1) 
    pim = Q[, , TT]
  else pim = rowSums(Q[, , TT])
  
  out = LMest:::rec3(Q, yv, PIs, Pio, pim)
  U = out$U
  V = out$V
  s1 = NULL
  if (k > 1) {
    out = LMest:::stationary(tau, k, G2, IPI)
    d21 = out$d0
    d22 = out$d1
    Mar = diag(k)
    u1 = Mar %*% rowSums(U[, , 1])
    V1 = Mar %*% V %*% t(Mar)
    s2 = d22 %*% (u1/la) + d21 %*% (V1[IPI]/PI[IPI])
  } else {
    s2 = NULL
  }
  s3 = NULL
  U = aperm(U, c(2, 1, 3))
  s4 = 0
  for (c in 1:k) {
    u = matrix(0, 1, k)
    u[c] = 1
    u = u[-1]
    D0 = cbind(I, matrix(u, nrow = 1) %x% one)
    agg = D0 %*% par0
    Eta1 = Eta01 + agg %*% matrix(1, 1, nd)
    Qv1 = LMest:::expit(Eta1)
    Qv1 = pmin(pmax(Qv1, 10^-100), 1 - 10^-100)
    Pit1 = lm %*% matrix(1, 1, nd) + Lm %*% Qv1
    Pit1 = pmin(pmax(Pit1, 10^-100), 1 - 10^-100)
    QQv1 = Qv1 * (1 - Qv1)
    DPv1 = 1/Pit1
    RRtc1 = array(0, c(np, lev, nd))
    for (j1 in 1:np) for (j2 in 1:lev) RRtc1[j1, j2, ] = QQv1[j1, ] * DPv1[j2, ]
    RRtc1 = RRtc1 * LLm1
    XXRi1 = array(0, c(dim(D0)[1], dim(D0)[2] + dim(Xd)[2], nd))
    for (h2 in 1:nd) {
      if (lev == 2) 
        XXRi1[, , h2] = c(D0, Xd[, , h2])
      else XXRi1[, , h2] = cbind(D0, Xd[, , h2])
    }
    XXRi1 = aperm(XXRi1, c(2, 1, 3))
    pc = U[, c, ]
    pc = as.vector(pc)
    nt = dim(Y)[1]
    YGP = matrix(Y, nt, n * TT) - Pit1[, as.vector(indn)]
    for (jd in 1:nd) {
      ind = INDN[[jd]]$ind
      pci = pc[ind]
      if (lev == 2) 
        XRi = (XXRi1[, , jd] %o% RRtc1[, , jd]) %*% GHt
      else XRi = (XXRi1[, , jd] %*% RRtc1[, , jd]) %*% GHt
      if (length(ind) == 1) {
        s4 = s4 + XRi %*% (YGP[, ind] * pci)
      } else {
        s4 = s4 + XRi %*% (YGP[, ind] %*% pci)
      }
    }
  }
  s = c(s1, s2, s3, s4)
  
  out = list(U = U, s = s)
  out
}


