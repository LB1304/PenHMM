Select.CV <- function(data, index, responsesFormula, k_values, gamma_values, n_folds = 5, 
                      start = 0, tol = 1e-4, maxit = 1e4, output = FALSE, out_se = FALSE, 
                      parallel = FALSE, workers = 1, future_seed = TRUE, se_multiplier = 0.1,
                      selection_rule = c("smallest_k_largest_gamma", "largest_gamma_smallest_k", "max_mean"), 
                      shuffle_folds = FALSE, seed = NULL, keep_fits = TRUE) {
  
  selection_rule <- match.arg(selection_rule)
  
  # ---------------------------------------------------------------------------
  # Basic checks
  # ---------------------------------------------------------------------------
  id_col <- which(names(data) == index[1])
  
  if (length(id_col) != 1) {
    stop("The first element of 'index' must identify exactly one column in 'data'.")
  }
  
  if (length(k_values) < 1) {
    stop("'k_values' must contain at least one value.")
  }
  
  if (length(gamma_values) < 1) {
    stop("'gamma_values' must contain at least one value.")
  }
  
  id <- data[, id_col]
  subject_ids <- unique(id)
  n_subjects <- length(subject_ids)
  
  if (n_folds > n_subjects) {
    stop("'n_folds' cannot be larger than the number of subjects.")
  }
  
  # ---------------------------------------------------------------------------
  # Fold construction at subject level
  # ---------------------------------------------------------------------------
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  ids_for_folds <- subject_ids
  
  if (shuffle_folds) {
    ids_for_folds <- sample(ids_for_folds)
  }
  
  fold_id <- as.integer(cut(seq_along(ids_for_folds),
                            breaks = n_folds,
                            labels = FALSE))
  
  fold_table <- data.frame(
    id = ids_for_folds,
    Fold = fold_id
  )
  
  # ---------------------------------------------------------------------------
  # Grid of models to estimate
  # ---------------------------------------------------------------------------
  grid_CV <- expand.grid(
    k = k_values,
    Gamma = gamma_values,
    Fold = seq_len(n_folds),
    KEEP.OUT.ATTRS = FALSE
  )
  
  grid_CV <- grid_CV[order(grid_CV$k, grid_CV$Gamma, grid_CV$Fold), ]
  rownames(grid_CV) <- NULL
  
  # ---------------------------------------------------------------------------
  # Internal functions
  # ---------------------------------------------------------------------------
  run_one_est_CV <- function(j) {
    
    k_j <- grid_CV$k[j]
    gamma_j <- grid_CV$Gamma[j]
    fold_j <- grid_CV$Fold[j]
    
    test_ids <- fold_table$id[fold_table$Fold == fold_j]
    data_train <- data[!(id %in% test_ids), , drop = FALSE]
    
    Estimate.HMCovManifest(
      data = data_train,
      index = index,
      responsesFormula = responsesFormula,
      k = k_j,
      start = start,
      tol = tol,
      maxit = maxit,
      output = output,
      out_se = out_se,
      ga2 = gamma_j
    )
  }
  
  compute_one_test_llk <- function(j, estimates) {
    
    k_j <- grid_CV$k[j]
    gamma_j <- grid_CV$Gamma[j]
    fold_j <- grid_CV$Fold[j]
    
    test_ids <- fold_table$id[fold_table$Fold == fold_j]
    data_test <- data[id %in% test_ids, , drop = FALSE]
    
    est_j <- estimates[[j]]
    
    Compute.Test.LLK(
      data = data_test,
      index = index,
      responsesFormula = responsesFormula,
      k = k_j,
      par = est_j$par,
      la = est_j$la,
      PI = est_j$PI,
      ga2 = gamma_j
    )
  }
  
  # ---------------------------------------------------------------------------
  # Apply function: sequential or parallel
  # ---------------------------------------------------------------------------
  if (parallel) {
    if (!requireNamespace("future", quietly = TRUE)) {
      stop("Package 'future' is required when parallel = TRUE.")
    }
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Package 'future.apply' is required when parallel = TRUE.")
    }
    
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    
    future::plan(future::multisession, workers = workers)
    
    apply_fun <- function(X, FUN) {
      future.apply::future_lapply(X, FUN, future.seed = future_seed)
    }
  } else {
    apply_fun <- function(X, FUN) {
      lapply(X, FUN)
    }
  }
  
  # ---------------------------------------------------------------------------
  # Model estimation
  # ---------------------------------------------------------------------------
  estimates <- apply_fun(seq_len(nrow(grid_CV)), run_one_est_CV)
  
  # ---------------------------------------------------------------------------
  # Test log-likelihood computation
  # ---------------------------------------------------------------------------
  test_llk <- apply_fun(
    seq_len(nrow(grid_CV)),
    function(j) compute_one_test_llk(j, estimates)
  )
  
  grid_results <- cbind(
    grid_CV,
    Test_LLK = unlist(test_llk)
  )
  
  # ---------------------------------------------------------------------------
  # Cross-validation summary
  # ---------------------------------------------------------------------------
  mean_table <- aggregate(
    Test_LLK ~ k + Gamma,
    data = grid_results,
    FUN = mean
  )
  
  se_table <- aggregate(
    Test_LLK ~ k + Gamma,
    data = grid_results,
    FUN = function(x) stats::sd(x) / sqrt(length(x))
  )
  
  cv_table <- merge(
    mean_table,
    se_table,
    by = c("k", "Gamma"),
    suffixes = c("_Mean", "_SE")
  )
  
  names(cv_table)[names(cv_table) == "Test_LLK_Mean"] <- "Mean_LLK"
  names(cv_table)[names(cv_table) == "Test_LLK_SE"] <- "SE_LLK"
  
  cv_table <- cv_table[order(cv_table$k, cv_table$Gamma), ]
  rownames(cv_table) <- NULL
  
  # ---------------------------------------------------------------------------
  # Selection rule
  # ---------------------------------------------------------------------------
  best_idx <- which.max(cv_table$Mean_LLK)
  best_mean <- cv_table$Mean_LLK[best_idx]
  se_at_best <- cv_table$SE_LLK[best_idx]
  threshold <- best_mean - se_multiplier * se_at_best
  
  cv_table$Best_Abs_LLK <- best_mean
  cv_table$SE_at_Best <- se_at_best
  cv_table$Threshold <- threshold
  cv_table$Is_Candidate <- cv_table$Mean_LLK >= threshold
  
  candidates <- cv_table[cv_table$Is_Candidate, , drop = FALSE]
  
  if (selection_rule == "smallest_k_largest_gamma") {
    candidates <- candidates[order(candidates$k, -candidates$Gamma), ]
    selected <- candidates[1, , drop = FALSE]
  }
  
  if (selection_rule == "largest_gamma_smallest_k") {
    candidates <- candidates[order(-candidates$Gamma, candidates$k), ]
    selected <- candidates[1, , drop = FALSE]
  }
  
  if (selection_rule == "max_mean") {
    selected <- cv_table[best_idx, , drop = FALSE]
  }
  
  # ---------------------------------------------------------------------------
  # Output
  # ---------------------------------------------------------------------------
  out <- list(
    selected_k = selected$k,
    selected_gamma = selected$Gamma,
    selected_row = selected,
    cv_table = cv_table,
    fold_table = fold_table,
    grid_results = grid_results,
    call = match.call()
  )
  
  if (keep_fits) {
    out$fits <- estimates
  }
  
  return(out)
}