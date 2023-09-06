
start_logger <- function(log_file, debug, verbose = TRUE){
  config_level <- if (debug) "debug" else "default"
  Sys.setenv(R_CONFIG_ACTIVE = config_level)
  if (is.null(log_file)){
    invisible(flog.appender(appender.console()))
  }else{
    invisible(flog.appender(appender.tee(log_file)))
  }
  invisible(flog.threshold(if(debug) "DEBUG" else "INFO"))
  if (!verbose) (flog.threshold("ERROR"))
}


lsqlincon_tmp <- function (C, d, A = NULL, b = NULL, Aeq = NULL, beq = NULL, lb = NULL, ub = NULL) {
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("quadprog needed for this function to work. Please install it.", 
         call. = FALSE)
  }
  stopifnot(is.numeric(C), is.numeric(d))
  if (is.null(A) && !is.null(b) || !is.null(A) && is.null(b)) 
    stop("If any, both 'A' and 'b' must be NULL.")
  if (is.null(Aeq) && !is.null(beq) || !is.null(Aeq) && is.null(beq)) 
    stop("If any, both 'Aeq' and 'beq' must be NULL.")
  if (!is.matrix(C)) 
    C <- matrix(C, 1)
  mc <- nrow(C)
  nc <- ncol(C)
  n <- nc
  if (length(d) != mc) 
    stop("Dimensions of 'C' and 'd' do not fit.")
  if (is.null(A) && is.null(Aeq) && is.null(lb) && is.null(ub)) 
    return(qr.solve(C, d))
  if (!is.null(A)) {
    if (!is.matrix(A)) 
      A <- matrix(A, 1)
    ma <- nrow(A)
    na <- ncol(A)
    if (na != n) 
      stop("Number of columns of 'A' does not fit with 'C'.")
    A <- -A
    b <- -b
  }
  if (is.null(Aeq)) {
    meq <- 0
  }
  else {
    if (!is.matrix(Aeq)) 
      Aeq <- matrix(Aeq, 1)
    meq <- nrow(Aeq)
    neq <- ncol(Aeq)
    if (neq != n) 
      stop("Number of columns of 'Aeq' does not fit with 'C'.")
  }
  if (is.null(lb)) {
    diag_lb <- NULL
  }
  else {
    if (length(lb) == 1) {
      lb <- rep(lb, n)
    }
    
    else if (length(lb) != n) {
      stop("Length of 'lb' and dimensions of C do not fit.")
    }
    diag_lb <- diag(n)
  }
  if (is.null(ub)) {
    diag_ub <- NULL
  }
  else {
    if (length(ub) == 1) {
      ub <- rep(ub, n)
    }
    else if (length(ub) != n) {
      stop("Length of 'ub' and dimensions of C do not fit.")
    }
    diag_ub <- -diag(n)
    ub <- -ub
  }
  # Dmat <- t(C) %*% C
  # dvec <- t(C) %*% d
  Amat <- rbind(Aeq, A, diag_lb, diag_ub)
  bvec <- c(beq, b, lb, ub)
  browser()
  rslt <- quadprog::solve.QP(Dmat, dvec, t(Amat), bvec, meq = meq)
  rslt$solution
}

init_means_vars <- function(C1_names, C2_names, feature_ids, source_ids, tau){
  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)
  k <- length(source_ids)
  m <- length(feature_ids)
  p1 <- length(C1_names)
  p2 <- length(C2_names)
  C1_names_ <- matrix("0",k*p1,1)
  if (p1){
    for (j in 1:k){
      C1_names_[((j-1)*p1+1):(p1*j)] <- unlist(lapply(C1_names,function(i) paste(source_ids[j],".",i,sep="")))
    }
  }
  # init
  mus_hat <- matrix(0, nrow=m, ncol=k)
  sigmas_hat <- matrix(0, nrow=m, ncol=k)
  deltas_hat <- matrix(0, nrow=m, ncol=p2)
  gammas_hat <- matrix(0, nrow=m, ncol=p1*k)
  deltas_hat_pvals <- matrix(1, nrow=m, ncol=p2)
  gammas_hat_pvals <- matrix(1, nrow=m, ncol=p1*k)
  gammas_hat_pvals.joint <- matrix(1, nrow=m, ncol=p1)
  if (is.null(tau)){
    tau_hat <- config[["tau_hat_init"]]
  }else{
    tau_hat <- tau
  }
  # set row and column names
  rownames(mus_hat) <- feature_ids
  colnames(mus_hat) <- source_ids
  rownames(sigmas_hat) <- feature_ids
  colnames(sigmas_hat) <- source_ids
  rownames(deltas_hat) <- feature_ids
  colnames(deltas_hat) <- C2_names
  rownames(gammas_hat) <- feature_ids
  colnames(gammas_hat) <- C1_names_
  rownames(deltas_hat_pvals) <- feature_ids
  colnames(deltas_hat_pvals) <- C2_names
  rownames(gammas_hat_pvals) <- feature_ids
  colnames(gammas_hat_pvals) <- C1_names_
  rownames(gammas_hat_pvals.joint) <- feature_ids
  colnames(gammas_hat_pvals.joint) <- C1_names
  
  return( list("mus_hat" = mus_hat, "sigmas_hat" = sigmas_hat, "gammas_hat" = gammas_hat, "deltas_hat"= deltas_hat, "deltas_hat_pvals" = deltas_hat_pvals, "gammas_hat_pvals" = gammas_hat_pvals, "gammas_hat_pvals.joint" = gammas_hat_pvals.joint, "tau_hat" = tau_hat) )
}

tca.validate_input <- function(X, W, C1, C1.map, C2, refit_W, refit_W.features, refit_W.sparsity, refit_W.sd_threshold, tau, constrain_mu, parallel, num_cores, max_iters, log_file, debug){
  
  flog.debug("Validating input types...")
  
  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)
  
  assert(is.matrix(X), "X must be of class 'matrix'")
  assert(is.matrix(W), "W must be of class 'matrix'")
  assert(is.null(C1) | is.matrix(C1), "C1 must be of class 'matrix' or NULL")
  assert(is.null(C1.map) | is.matrix(C1.map), "C1.map must be of class 'matrix' or NULL")
  assert(is.null(C2) | is.matrix(C2), "C2 must be of class 'matrix' or NULL")
  assert(is.null(refit_W.features) | is.character(refit_W.features), "refit_W.features must be of class 'character' or NULL")
  assert(is.null(tau) | is.numeric(tau), "tau must be of class 'numeric' or NULL")
  
  assert(is.numeric(refit_W.sparsity), "refit_W.sparsity must be of class 'numeric'")
  assert(is.numeric(refit_W.sd_threshold), "refit_W.sparsity must be of class 'numeric'")
  assert(is.numeric(max_iters), "max_iters must be of class 'numeric'")
  assert(is.logical(refit_W), "refit_W must be of class 'logical'")
  assert(is.logical(parallel), "parallel must be of class 'logical'")
  assert(is.logical(debug), "debug must be of class 'logical'")
  assert(is.character(log_file) | is.null(log_file), "log_file must be of class 'character' or NULL")
  
  flog.debug("Validating input stucture and values...")
  assert(!is.null(rownames(X)) & !is.null(colnames(X)), "X must have row names and column names")
  assert(!is.null(rownames(W)) & !is.null(colnames(W)), "W must have row names and column names")
  if (!is.null(C1)) assert(!is.null(rownames(C1)) & !is.null(colnames(C1)), "C1 must have row names and column names")
  if (!is.null(C2)) assert(!is.null(rownames(C2)) & !is.null(colnames(C2)), "C2 must have row names and column names")
  
  flog.debug("Validating input conditions...")
  if (refit_W) assert(refit_W.sparsity <= nrow(X) , "argument refit_W.sparsity must satisfy refit_W.sparsity <= nrow(X)")
  if (refit_W) assert(refit_W.sd_threshold >= 0 , "argument refit_W.sd_threshold must satisfy refit_W.sd_threshold >= 0")
  if (!(is.null(C1.map))) assert(constrain_mu , "argument C1.map cannot be useed with constrain_mu set to FALSE")
  
  flog.debug("Validating matrix dimensions...")
  assert(dim(X)[2] == dim(W)[1] , "The number of columns in X is inconsistent with the number of rows in W")
  if (!is.null(C1)) assert(dim(X)[2] == dim(C1)[1] , "The number of columns in X is inconsistent with the number of rows in C1")
  if (!is.null(C2)) assert(dim(X)[2] == dim(C2)[1] , "The number of columns in X is inconsistent with the number of rows in C2")
  
  flog.debug("Validating the order of observations across matrices...")
  assert(all(colnames(X) == rownames(W)), "The order of observations in W (in the rows) must match the order of the observations in X (in the columns).")
  if (!is.null(C1)) assert(all(colnames(X) == rownames(C1)), "The order of observations in C1 (in the rows) must match the order of the observations in X (in the columns).")
  if (!is.null(C2)) assert(all(colnames(X) == rownames(C2)), "The order of observations in C2 (in the rows) must match the order of the observations in X (in the columns).")
  
  flog.debug("Validating that W is non-negative and each row sums up to 1...")
  assert(all(W >= 0), "The entries of W must be non-negative.")
  assert(all(abs(rowSums(W) - 1) < 0.0001), "Each row in W must sum up to 1.")
  
  if (!is.null(C1.map)) assert(sum(C1.map == 1 | C1.map == 0) == ncol(C1.map)*nrow(C1.map), "The entries of C1.map must all be 0 or 1")
  
  if (!is.null(tau)) assert(tau >= 0, "tau must be non-negative")
  
  th <- config[["min_sd"]]**2
  assert(sum(rowVars(X) < th) == 0, paste("X must not include features with variance less than ",as.character(th),sep=""))
  
}