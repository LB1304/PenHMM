Estimate.HMCovManifest <- function(data, responsesFormula, index, k, start = 0, ga1 = 0, ga2 = 0, 
                                   tol = 1e-8, maxit = 1e4, out_se = FALSE, output = FALSE) {
  
  prepared_data <- PrepareData(
    data = data,
    responsesFormula = responsesFormula,
    index = index
  )
  
  out <- EstimateModel(
    S = as.matrix(prepared_data$Y[, , 1]),
    X = prepared_data$X,
    yv = prepared_data$freq,
    k = k,
    tol = tol,
    maxit = maxit,
    start = start,
    output = output,
    out_se = out_se,
    ga1 = ga1,
    ga2 = ga2
  )
  
  out <- c(out, call = match.call())
  attributes(out)$responsesFormula <- responsesFormula
  attributes(out)$whichid <- prepared_data$whichid
  attributes(out)$whichtv <- prepared_data$whichtv
  attributes(out)$id <- prepared_data$id
  attributes(out)$time <- prepared_data$time
  
  return(out)
}


