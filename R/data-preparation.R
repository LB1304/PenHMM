PrepareData <- function(data, responsesFormula, index) {
  
  id.which <- which(names(data) == index[1])
  tv.which <- which(names(data) == index[2])
  
  id <- data[, id.which]
  tv <- data[, tv.which]
  data.new <- data[, -c(id.which, tv.which), drop = FALSE]
  
  formula <- Formula:::Formula(responsesFormula)
  formula.length <- length(formula)
  
  Y <- Formula::model.part(
    formula,
    data = model.frame(formula, data = data.new, na.action = NULL),
    lhs = 1
  )
  Y <- data.matrix(Y)
  Y_names <- colnames(Y)
  
  X <- NULL
  if (formula.length[2] != 0) {
    X <- model.matrix(
      formula,
      model.frame(formula = formula, data.new, na.action = NULL)
    )
    X <- data.matrix(X)[, -1]
  }
  
  idu <- unique(id)
  n <- length(idu)
  TT <- max(tv)
  
  X <- as.matrix(X)
  nxMan <- ncol(X)
  XX <- array(NA, c(n, TT, nxMan))
  
  Y <- as.matrix(Y)
  ny <- ncol(Y)
  YY <- array(NA, c(n, TT, ny))
  
  for (i in seq_len(n)) {
    ind <- which(id == idu[i])
    tmp <- 0
    
    for (t in tv[ind]) {
      tmp <- tmp + 1
      XX[i, t, ] <- X[ind[tmp], ]
      YY[i, t, ] <- Y[ind[tmp], ]
    }
  }
  
  freq <- rep(1, nrow(YY))
  
  Y <- YY
  X <- XX
  
  if (min(Y, na.rm = TRUE) > 0) {
    for (j in seq_len(dim(Y)[3])) {
      Y[, , j] <- Y[, , j] - min(Y[, , j], na.rm = TRUE)
    }
  }
  
  dimnames(Y)[[3]] <- list(Y_names)
  
  return(list(Y = Y, X = X, freq = freq))
}


