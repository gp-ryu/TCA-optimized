source('~/.scripts/zz_header.r')
library(TCA)
library(Rcpp)
library(RcppEigen)
library(matrixcalc)
library(pracma)
library(pbmcapply)
library(futile.logger)
library(nloptr)
library(matrixStats)

sourceCpp('anRpackage/src/matrixMultiplications.cpp')
sourceCpp('anRpackage/src/estimate_Z_j.cpp')
cores = 16

#load('../hannum.chr22.RData')
args <- commandArgs(trailingOnly = T)

#tca.mdl <- tca(hannum$X, hannum$W, C1 = hannum$cov[,c('age','gender')], C2 = hannum$cov[,3:ncol(hannum$cov)])

#### tca_cpp ####

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

minus_log_likelihood_tau <- function(U,W_squared,sigmas,const,tau){
  m <- ncol(U)
  res <- matrix(0,m,2)
  tmp <- lapply(1:m, function(j)
  {V <- tcrossprod(W_squared,t(sigmas[j,]**2))+tau**2;
  V_squared <- V**2;
  return (c(-0.5*(const-sum(log(V))-sum(U[,j]/V)), -(tau*(sum(U[,j]/V_squared)-sum(1./V))))) } )
  for (j in 1:m){
    res[j,] = tmp[[j]]
  }
  res <- colSums(res)
  return(list( "objective" = res[1], "gradient" = res[2] ))
}

# Returns the (minus) log likelihood of the model for a particular feature j in a list together with the gradient with respect to sigmas of feature j
# Input:
#  sigmas_hat (of a particular feature j)
#  U_j is the j-th columns of U, where U = (tcrossprod(W,mus_hat) + tcrossprod(C2,deltas_hat) + tcrossprod(C1_,gammas_hat) - X)**2
#  const <- -n*log(2*pi)
#  W_squared <- W**2
#  tau
minus_log_likelihood_sigmas <- function(sigmas,U_j,W_squared,const,tau){
  k <- length(sigmas)
  n <- nrow(W_squared)
  V <- tcrossprod(W_squared,t(sigmas**2))+tau**2
  V_squared <- V**2
  return(list("objective"= -0.5*(const-sum(log(V))-sum(U_j/V)),
              "gradient" = -(colSums(W_squared*repmat(sigmas,n,1)*t(repmat(U_j,k,1))/repmat(V_squared,1,k)) - colSums(W_squared*repmat(sigmas,n,1)/repmat(V,1,k)))))
}


# Returns the (minus) log likelihood of the model for a particular observation i in a list together with the gradient with respect to w_i
# Input:
# C_tilde - if p1=0 then C_tilde = matrix(0,m,k), otherwise it is calculated as follows:
# for (h in 1:k){
#    C_tilde[,h] = tcrossprod(gammas[,(1+(h-1)*p1):(h*p1)],c1_i)
# }
# sigmas_squared = sigmas**2
# crossprod_deltas_c2_i = tcrossprod(deltas,t(C2[i,]))
minus_log_likelihood_w <- function(w_i, x_i, c1_i, mus, tau, gammas, const, C_tilde, sigmas_squared, crossprod_deltas_c2_i){
  k <- length(w_i)
  m <- dim(mus)[1]
  c1_i_ <- calc_C1_W_interactions(w_i,c1_i)
  V <- tcrossprod(sigmas_squared,w_i**2)+tau**2
  V_rep <- repmat(V,1,k)
  U_i <- tcrossprod(mus,w_i) + crossprod_deltas_c2_i + tcrossprod(gammas,c1_i_) - t(x_i)
  U_i_squared <- U_i**2
  w_i_rep <- repmat(w_i,m,1)
  fval <- -0.5*(const-sum(log(V))-sum(U_i_squared/V))
  gval <- colSums(w_i_rep*sigmas_squared/V_rep) + colSums(( (mus+C_tilde)*repmat(U_i,1,k)*V_rep - w_i_rep*sigmas_squared*repmat(U_i_squared,1,k) ) / repmat(V**2,1,k))
  return(list("objective"= fval, "gradient" = gval))
}

#' @importFrom matrixcalc hadamard.prod
calc_C1_W_interactions <- function(W,C1){
  n <- nrow(W)
  k <- ncol(W)
  p1 <- ncol(C1)
  if (p1){
    return( hadamard.prod(Reshape(Reshape(apply(W, 2, function(v) repmat(v,1,p1)), n*p1*k,1), n,p1*k), repmat(C1, 1, k)) )
  }else{
    return(matrix(0,n,0))
  }
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
#' @importFrom pracma repmat
#' @importFrom pracma lsqlincon
#' @importFrom nloptr nloptr
#' @importFrom matrixStats colVars
tca.fit_means_vars <- function(X, W, mus_hat, sigmas_hat, tau_hat, C2, deltas_hat, C1, gammas_hat, C1.map, tau, vars.mle, constrain_mu, max_iters, parallel, num_cores){
  
  flog.debug("Starting function 'tca.fit_means_vars'...")
  
  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)
  nloptr_opts = list("algorithm"=config[["nloptr_opts_algorithm"]], "xtol_rel"=config[["nloptr_opts_xtol_rel"]], "print_level" = config[["nloptr_opts_print_level"]], "check_derivatives" = config[["nloptr_opts_check_derivatives"]])
  
  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(W)
  p1 <- ncol(C1)
  p2 <- ncol(C2)
  
  C1_ <- calc_C1_W_interactions(W,C1)
  
  cl <- if (parallel) init_cluster(num_cores) else NULL
  
  const <- -n*log(2*pi)
  W_squared <- W**2
  W_squared_ <- cbind(W_squared,numeric(n) + 1)
  
  ll_prev <- -Inf
  # Perform an alternative optimization of the means (mus, deltas, gammas) and variances (sigmas and tau)
  for (iter in 1:max_iters){
    
    flog.info("Iteration %s out of %s internal iterations...", iter, max_iters)
    
    # (1) Estimate the means (mus, deltas, gammas)
    
    ub <- numeric(k+p2+p1*k)+config[["lsqlincon_inf"]]
    lb <- numeric(k+p2+p1*k)-config[["lsqlincon_inf"]]
    if (constrain_mu){
      ub[1:k] <- max(X)
      ub[1:k] <- ub[1:k] - config[["mu_epsilon"]]
      lb[1:k] <- min(X)
      lb[1:k] <- lb[1:k] + config[["mu_epsilon"]]
    }
    
    if (p1){
      # All zeros will get ub and lb set to 0 (i.e. effect sizes will not be estimated); use a small value (config[["nloptr_opts_xtol_rel"]]) for stability; otherwise nloptr may return an error in some cases.
      ub[seq(k+p2+1,k+p2+p1*k)] <- Reshape(C1.map,p1*k,1)*config[["lsqlincon_inf"]] + config[["nloptr_opts_xtol_rel"]]
      lb[seq(k+p2+1,k+p2+p1*k)] <- Reshape(C1.map,p1*k,1)*(-config[["lsqlincon_inf"]]) - config[["nloptr_opts_xtol_rel"]]
    }
    
    flog.debug("Calculate W_norms and related quantities")
    if (sum(colSums(mus_hat)) == 0){
      # Use the following for getting initial estimates of mus, deltas, gammas; under the assumptions that tau=0 and sigmas_{1j},...,sigmas_{kj} for each j.
      W_norms <- rowSums(W**2)**0.5
      # Since W_norms is the same for all features in this case can already calculate the following quantities
      W_tilde <- W/t(repmat(W_norms,k,1))
      C1_tilde <- if (p1>0) C1_/t(repmat(W_norms,k*p1,1)) else C1_
      C2_tilde <- if (p2>0) C2/t(repmat(W_norms,p2,1)) else C2
    }else{
      flog.debug("Calculate W_norms")
      W_norms <- (tcrossprod(W**2,sigmas_hat**2)+tau_hat**2 )**0.5
    }
    X_tilde <- X/W_norms
    
    flog.debug("Estimate mus, deltas, gammas.")
    
    if (sum(colSums(mus_hat)) == 0){
      #if (parallel) clusterExport(cl, varlist = c("W_tilde","C2_tilde","C1_tilde","X_tilde","lb","ub","k","p1","p2","lsqlincon"), envir=environment())
      res <- pbmclapply(ignore.interactive = T, 1:m, mc.cores = cores, function(j) lsqlincon(cbind(W_tilde,C2_tilde,C1_tilde), X_tilde[,j], lb = lb,  ub = ub))
    }else{
      #if (parallel) clusterExport(cl, c("W_norms","W","C2","C1_","X_tilde","lb","ub","k","p1","p2","lsqlincon"), envir=environment())
      res <- pbmclapply(ignore.interactive = T, 1:m,mc.cores = cores, function(j) lsqlincon(cbind(W/t(repmat(W_norms[,j],k,1)),
                                                      if (p2>0) C2/t(repmat(W_norms[,j],p2,1)) else C2,
                                                      if (p1>0) C1_/t(repmat(W_norms[,j],k*p1,1)) else C1_),
                                                X_tilde[,j], lb = lb,  ub = ub))
    }
    
    # Update estimtes
    tic(); flog.info('update estimates1')
    
    for (j in 1:m){
      mus_hat[j,] = res[[j]][1:k]
      deltas_hat[j,seq(1,p2,length=p2)] = res[[j]][seq(k+1,k+p2,length=p2)]
      gammas_hat[j,seq(1,k*p1,length=k*p1)] = res[[j]][seq(k+p2+1,k+p2+p1*k,length=p1*k)]
    }
    toc()
    
    # (2) Estimate the variances (sigmas, tau)
    
    # Calculate some quantities that will be repeatedly used throughout the optimization in this step
    U <- (tcrossprod(W,mus_hat) + tcrossprod(C2,deltas_hat) + tcrossprod(C1_,gammas_hat) - X)**2
    
    if (vars.mle){
      
      if (sum(colSums(sigmas_hat)) == 0){
        # Set a starting point for the optimization
        flog.debug("Get initial estimates of sigmas")
        row_names <- rownames(sigmas_hat)
        col_names <- colnames(sigmas_hat)
        sigmas_hat <- t(repmat((colVars(X)/k)**0.5,k,1))
        rownames(sigmas_hat) <- row_names
        colnames(sigmas_hat) <- col_names
      }
      
      # (2.2) Estimate sigmas
      flog.debug("Estimate sigmas.")
      lb <- numeric(k)+config[["min_sd"]]
      ub <- numeric(k)+Inf
      if (parallel) clusterExport(cl, c("lb","ub","n","k","U","const","W_squared","sigmas_hat","tau_hat","nloptr_opts","minus_log_likelihood_sigmas"), envir=environment())
      res <- pbmclapply(ignore.interactive = T, 1:m,mc.cores = cores, function(j) nloptr( x0=sigmas_hat[j,],
                                              eval_f = function(x,U_j,W_squared,const,tau_hat) minus_log_likelihood_sigmas(x,U_j,W_squared,const,tau_hat),
                                              lb = lb,
                                              ub = ub,
                                              opts = nloptr_opts,
                                              U_j = U[,j],
                                              W_squared = W_squared,
                                              const = const,
                                              tau_hat = tau_hat)$solution)
      
      tic(); flog.info('update estimates')
      for (j in 1:m){
        sigmas_hat[j,] = res[[j]]
      }
      toc()
      
      # (2.3) Estimate tau
      tic(); flog.info('estimate tau')
      if (is.null(tau)){
        flog.debug("Estimate tau.")
        lb <- config[["min_sd"]]
        ub <- Inf
        tau_hat = nloptr(x0=tau_hat,
                         eval_f = function(x,U,W_squared,sigmas_hat,const) minus_log_likelihood_tau(U,W_squared,sigmas_hat,const,x),
                         lb = lb,
                         ub = ub,
                         opts = nloptr_opts,
                         U = U,
                         W_squared = W_squared,
                         sigmas_hat = sigmas_hat,
                         const = const)$solution
      }
      
      
    }else{
      # Learn the variances using gmm
      flog.debug("Estimate sigmas, tau")
      lb <- numeric(k+1)+config[["min_sd"]]
      # calculate weights matrix V
      if (sum(colSums(sigmas_hat)) == 0){
        V <- matrix(1,n,m)
      }else{
        V <- abs(U - tcrossprod(W_squared_, cbind(sigmas_hat**2, matrix(tau_hat**2,m,1))))
        V[V < config[["V_weight_limit"]]] <- config[["V_weight_limit"]]
      }
      if (parallel) clusterExport(cl, c("lb","U","W_squared_","lsqlincon","V","n"), envir=environment())
      res <- pbmclapply(ignore.interactive = T, 1:m,mc.cores = cores, function(j) {
        x <- W_squared_/t(repmat(V[,j],ncol(W_squared_),1));
        # For numeric stability, normalize the design matrix and adjust the final estimats accordingly
        norms <- (colSums(x**2))**0.5;
        x <- x/repmat(norms,n,1);
        lsqlincon(x, U[,j]/V[,j], lb = lb*norms)/norms
      })
      tau_squared_hat <- 0
      for (j in 1:m){
        sigmas_hat[j,] <- sqrt(res[[j]][1:k])
        tau_squared_hat <- tau_squared_hat + res[[j]][k+1]
      }
      tau_hat <- sqrt(tau_squared_hat/m)
    }
    toc()
    
    flog.debug("Test for convergence.")
    ll_new <- -minus_log_likelihood_tau(U,W_squared,sigmas_hat,const,tau_hat)[[1]]
    # Test for convergence
    ll_diff = ll_new-ll_prev
    flog.debug("ll_new=%s, ll_prev=%s, ll_diff=%s, diff_threshold=%s",ll_new,ll_prev,ll_diff,config[["epsilon"]]*abs(ll_new))
    if (ll_diff < config[["epsilon"]]*abs(ll_new)){
      flog.info("Internal loop converged.")
      flog.debug("break")
      break
    }
    ll_prev <- ll_new
    
  }
  
  if (parallel) stop_cluster(cl)
  
  return(list("mus_hat"=mus_hat, "sigmas_hat"=sigmas_hat, "tau_hat"=tau_hat, "deltas_hat"=deltas_hat, "gammas_hat"=gammas_hat, "tau_hat"=tau_hat))
  
}


tca.fit_W <- function(X, W, mus_hat, sigmas_hat, tau_hat, C2, deltas_hat, C1, gammas_hat, max_iters, parallel, num_cores){
  
  flog.debug("Starting function 'tca.fit_W'...")
  
  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)
  nloptr_opts = list("algorithm"=config[["nloptr_opts_fit_W_algorithm"]], "xtol_rel"=config[["nloptr_opts_xtol_rel"]], "print_level" = config[["nloptr_opts_print_level"]], "check_derivatives" = config[["nloptr_opts_check_derivatives"]])
  
  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(W)
  p1 <- ncol(C1)
  lb <- numeric(k)
  ub <- numeric(k)+1
  ones <- numeric(k)+1
  sigmas_squared <- sigmas_hat**2
  W_hat <- matrix(0, n, k)
  colnames(W_hat) <- colnames(W)
  rownames(W_hat) <- rownames(W)
  const <- -m*log(2*pi)
  
  cl <- if (parallel) init_cluster(num_cores) else NULL
  if (parallel) clusterExport(cl, c("lb","ub","ones","nloptr_opts","X","W","C1","C2","n","k","p1","m","const","tau_hat","mus_hat","gammas_hat","deltas_hat","sigmas_squared","minus_log_likelihood_w","calc_C1_W_interactions"), envir=environment())
  res <- pbmclapply(ignore.interactive = T, 1:n,mc.cores = cores, function(i) nloptr( x0=W[i,],
                                          eval_f = function(x, x_i, c1_i, mus, tau, gammas, const, C_tilde, sigmas_squared, crossprod_deltas_c2_i) minus_log_likelihood_w(t(x), x_i, c1_i, mus, tau, gammas, const, C_tilde, sigmas_squared, crossprod_deltas_c2_i),
                                          eval_g_eq = function(x, x_i, c1_i, mus, tau, gammas, const, C_tilde, sigmas_squared, crossprod_deltas_c2_i) list("constraints" = crossprod(x,ones)-1, "jacobian" = ones),
                                          lb = lb,
                                          ub = ub,
                                          opts = nloptr_opts,
                                          x_i = t(X[i,]),
                                          c1_i = t(C1[i,]),
                                          mus = mus_hat,
                                          tau = tau_hat,
                                          gammas = gammas_hat,
                                          const = const,
                                          C_tilde = if (p1>0) apply(as.matrix(1:k),1,function(h) tcrossprod(gammas_hat[,(1+(h-1)*p1):(h*p1)],t(C1[i,]))) else matrix(0,m,k),
                                          sigmas_squared = sigmas_squared,
                                          crossprod_deltas_c2_i = tcrossprod(deltas_hat,t(C2[i,])) )$solution )
  if (parallel) stop_cluster(cl)
  for (i in 1:n){
    W_hat[i,] = res[[i]]
  }
  return(W_hat)
}


tca.fit <- function(X, W, C1, C1.map, C2, refit_W, tau, vars.mle, constrain_mu, parallel, num_cores, max_iters){
  
  flog.debug("Starting function 'tca.fit'...")
  
  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)
  
  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(W)
  p1 <- ncol(C1)
  p2 <- ncol(C2)
  
  # init estimates
  init <- init_means_vars(colnames(C1), colnames(C2), colnames(X), colnames(W), tau)
  mus_hat <- init[["mus_hat"]]
  sigmas_hat <- init[["sigmas_hat"]]
  deltas_hat <- init[["deltas_hat"]]
  gammas_hat <- init[["gammas_hat"]]
  tau_hat <- init[["tau_hat"]]
  deltas_hat_pvals <- if(constrain_mu) NULL else init[["deltas_hat_pvals"]]
  gammas_hat_pvals <- if(constrain_mu) NULL else init[["gammas_hat_pvals"]]
  gammas_hat_pvals.joint <- if(constrain_mu) NULL else init[["gammas_hat_pvals.joint"]]
  
  const <- -n*log(2*pi)
  ll_prev <- -Inf
  for (iter in 1:max_iters){
    if (refit_W) flog.info("Iteration %s out of %s external iterations (fitting all parameters including W)...", iter, max_iters)
    
    flog.info("Fitting means and variances...")
    res <- tca.fit_means_vars(X, W, mus_hat, sigmas_hat, tau_hat, C2, deltas_hat, C1, gammas_hat, C1.map, tau, vars.mle, constrain_mu, max_iters, parallel, num_cores)
    mus_hat <- res[["mus_hat"]]
    sigmas_hat <- res[["sigmas_hat"]]
    deltas_hat <- res[["deltas_hat"]]
    gammas_hat <- res[["gammas_hat"]]
    tau_hat <- res[["tau_hat"]]
    
    if (!refit_W) break
    flog.info("Fitting W...")
    W <- tca.fit_W(X, W, mus_hat, sigmas_hat, tau_hat, C2, deltas_hat, C1, gammas_hat, max_iters, parallel, num_cores)
    
    C1_ <- calc_C1_W_interactions(W,C1)
    flog.debug("Test for convergence.")
    U <- (tcrossprod(W,mus_hat) + tcrossprod(C2,deltas_hat) + tcrossprod(C1_,gammas_hat) - X)**2
    W_squared <- W**2
    ll_new <- -minus_log_likelihood_tau(U,W_squared,sigmas_hat,const,tau_hat)[[1]]
    flog.debug("~~Main loop ll=%s",ll_new)
    # Test for convergence
    ll_diff = ll_new-ll_prev
    flog.debug("ll_new=%s, ll_prev=%s, ll_diff=%s, threshold=%s",ll_new,ll_prev,ll_diff,config[["epsilon"]]*abs(ll_new))
    if (ll_diff < config[["epsilon"]]*abs(ll_new)){
      flog.info("External loop converged.")
      flog.debug("break")
      break
    }
    ll_prev <- ll_new
    
  }
  
  if(!constrain_mu){
    flog.info("Calculate p-values for deltas and gammas.")
    C1_ <- calc_C1_W_interactions(W,C1)
    W_norms <- (tcrossprod(W**2,sigmas_hat**2)+tau_hat**2 )**0.5
    X_tilde <- X/W_norms
    cl <- NULL
    if (parallel){
      cl <- init_cluster(num_cores)
      clusterExport(cl, varlist = c("C1_","W_norms","X_tilde","W","k","p2","C2","p1","p1","lm"), envir=environment())
    }
    
    #### opt needed ####
    
    # res <- pbmclapply(ignore.interactive = T, 1:m, mc.cores = cores, function(j) {
    res <- ({j = 1
    
      colnames(C1_) <- c(seq(1:ncol(C1_)) %>% paste0('V',.))
      df = data.frame(y = X_tilde[,j],
                      cbind(W/t(repmat(W_norms[,j],k,1)),
                            if (p2>0) C2/t(repmat(W_norms[,j],p2,1)) else C2,
                            if (p1>0) C1_/t(repmat(W_norms[,j],k*p1,1)) else C1_)
                      );
      df_y = X_tilde[,j] %>% as.matrix
      df_x = matrix(c(W/t(repmat(W_norms[,j],k,1)),
                      if (p2>0) C2/t(repmat(W_norms[,j],p2,1)) else C2,
                      if (p1>0) C1_/t(repmat(W_norms[,j],k*p1,1)) else C1_),
                    nrow = nrow(df_y))
      colnames(df_x) <- c(colnames(W), colnames(C2), colnames(C1_))
      rownames(df_x) <- rownames(df_y)
      
      
      mdl1.fit <- lm(y~., data = df)
      mdl1.coef <- summary(mdl1.fit)$coefficients
      mdl1.cov.names <- colnames(df)[colnames(df) != 'y']

      
      mdl2.fit <- .lm.fit(df_x, df_y)
      mdl2.coef <- summary(mdl1.fit$coefficients)
      mdl2.cov.names <- mdl2.coef %>% rownames %>% .[-1]
      
      browser()
      deltas_gammas_hat_pvals <- sapply(mdl1.cov.names, function(x){
        if (x %in% rownames(mdl1.coef)){
          return(mdl1.coef[x,'Pr(>|t|)'])
        }
        else {
          return(NA)
        }
      })
      
      #deltas_gammas_hat_pvals <- summary(mdl1.fit)$coefficients[2:(1+k+p1*k+p2),4];
      gammas_hat_pvals.joint <- numeric(p1)+1
      if (p1){
        C1_alt <- C1_/t(repmat(W_norms[,j],k*p1,1))
        for (d in 1:p1){
          C1_null <- C1_alt[,setdiff(1:(p1*k),seq(d,k*p1,p1))]
          df = data.frame(y = X_tilde[,j], cbind(W/t(repmat(W_norms[,j],k,1)),if (p2>0) C2/t(repmat(W_norms[,j],p2,1)) else C2,C1_null));
          mdl0.fit <- lm(y~., data = df)
          anova.fit <- anova(mdl0.fit,mdl1.fit);
          gammas_hat_pvals.joint[d] <- anova.fit$`Pr(>F)`[2];
        }
      }
      return(c(deltas_gammas_hat_pvals,gammas_hat_pvals.joint));
      
      } )
    ####  
    
    if (parallel) stop_cluster(cl)
    for (j in 1:m){
      deltas_hat_pvals[j,] <- res[[j]][(k+1):(k+p2)]
      gammas_hat_pvals[j,] <- res[[j]][(k+p2+1):(k+p2+p1*k)]
      gammas_hat_pvals.joint[j,] <- res[[j]][(k+p2+p1*k+1):(k+p2+p1*k+p1)]
    }
    # res <- pblapply(1:m,function(j) {
    #   df = data.frame(y = X_tilde[,j], cbind(W/t(repmat(W_norms[,j],k,1)),if (p2>0) C2/t(repmat(W_norms[,j],p2,1)) else C2,if (p1>0) C1_/t(repmat(W_norms[,j],k*p1,1)) else C1_));
    #   return(summary(lm(y~., data = df))$coefficients[2:(1+k+p1*k+p2),4]);
    #   }, cl = cl )
    # if (parallel) stop_cluster(cl)
    # for (j in 1:m){
    #   deltas_hat_pvals[j,] <- res[[j]][(k+1):(k+p2)]
    #   gammas_hat_pvals[j,] <- res[[j]][(k+p2+1):(k+p2+p1*k)]
    # }
  }
  return (list("W" = W, "mus_hat" = mus_hat, "sigmas_hat" = sigmas_hat, "tau_hat" = tau_hat, "deltas_hat" = deltas_hat, "gammas_hat" = gammas_hat, "deltas_hat_pvals" = deltas_hat_pvals, "gammas_hat_pvals"= gammas_hat_pvals, "gammas_hat_pvals.joint" = gammas_hat_pvals.joint))
}

tca_pbmc <- function(X, W, C1 = NULL, C1.map = NULL, C2 = NULL, refit_W = FALSE, refit_W.features = NULL, refit_W.sparsity = 500, refit_W.sd_threshold = 0.02, tau = NULL, vars.mle = FALSE, constrain_mu = FALSE, parallel = FALSE, num_cores = NULL, max_iters = 10, log_file = "TCA.log", debug = FALSE, verbose = TRUE){
  
  start_logger(log_file, debug, verbose)
  
  flog.info("Starting tca...")
  
  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)
  
  # progress bar options
  #op <- if (verbose) pboptions(nout = config$nout, type = config[["type"]]) else pboptions(nout = config$nout, type = "none")
  
  X <- if (is.matrix(X)) X else as.matrix(X)
  W <- if (is.matrix(W)) W else as.matrix(W)
  C1 <- if (is.matrix(C1) | is.null(C1)) C1 else as.matrix(C1)
  C2 <- if (is.matrix(C2) | is.null(C2)) C2 else as.matrix(C2)
  
  flog.info("Validating input...")
  tca.validate_input(X %>% head, W, C1, C1.map, C2, refit_W, refit_W.features, refit_W.sparsity, refit_W.sd_threshold, tau, constrain_mu, parallel, num_cores, max_iters, log_file, debug)
  if (is.null(C1)) C1 <- matrix(0, nrow=ncol(X), ncol=0)
  if (is.null(C2)) C2 <- matrix(0, nrow=ncol(X), ncol=0)
  
  C1.map <- if (is.null(C1.map)) matrix(1,ncol(C1),ncol(W)) else C1.map
  
  msg <- "Fitting the TCA model..."
  if (refit_W){
    flog.info("Starting re-estimation of W...")
    if (is.null(refit_W.features)){
      flog.info("Performing feature selection using refactor...")
      ref <- refactor.run(X = X, k = ncol(W), sparsity = refit_W.sparsity, C = if (ncol(cbind(C1,C2))) cbind(C1,C2) else NULL, C.remove = FALSE, sd_threshold = refit_W.sd_threshold, num_comp = NULL, rand_svd = config$rand_svd, log_file = NULL, logger_on = FALSE, verbose = verbose)
      refit_W.features <- ref$ranked_list[1:refit_W.sparsity]
    }
    X_sub <- subset(X, subset = rownames(X) %in% refit_W.features)
    flog.info("Fitting the TCA model using the selected features for re-estimating W...")
    mdl0 <- tca.fit(X = t(X_sub), W = W, C1 = C1, C1.map = C1.map, C2 = C2, refit_W = TRUE, tau = tau, vars.mle = vars.mle, constrain_mu = constrain_mu, parallel = parallel, num_cores = num_cores, max_iters = max_iters)
    W <- mdl0[["W"]]
    msg <- "Fitting the TCA model given the updated W..."
  }
  X <- t(X)
  flog.info(msg)
  mdl <- tca.fit(X = X, W = W, C1 = C1, C1.map = C1.map, C2 = C2, refit_W = FALSE, tau = tau, vars.mle = vars.mle, constrain_mu = constrain_mu, parallel = parallel, num_cores = num_cores, max_iters = max_iters)
  W <- mdl[["W"]]
  mus_hat <- mdl[["mus_hat"]]
  sigmas_hat <- mdl[["sigmas_hat"]]
  tau_hat <- mdl[["tau_hat"]]
  deltas_hat <- mdl[["deltas_hat"]]
  gammas_hat <- mdl[["gammas_hat"]]
  deltas_hat_pvals <- mdl[["deltas_hat_pvals"]]
  gammas_hat_pvals <- mdl[["gammas_hat_pvals"]]
  gammas_hat_pvals.joint <- mdl[["gammas_hat_pvals.joint"]]
  flog.info("Finished tca.")
  return (list("W" = W, "mus_hat" = mus_hat, "sigmas_hat" = sigmas_hat, "tau_hat" = tau_hat, "deltas_hat" = deltas_hat, "gammas_hat" = gammas_hat, "deltas_hat_pvals" = deltas_hat_pvals, "gammas_hat_pvals" = gammas_hat_pvals, "gammas_hat_pvals.joint" = gammas_hat_pvals.joint, "C1" = C1, "C2" = C2) )
}



#### Misc ####
#RcppEigen.package.skeleton()
# 
# X = hannum$X %>% t
# W = hannum$W
# C1= tca.mdl$C1
# n = dim(X)[1]
# m <- dim(X)[2]
# k = ncol(hannum$W)
# p1= ncol(C1)
# tau_hat <- tca.mdl$tau_hat
# 
# mus_hat_j = tca.mdl$mus_hat[1,]
# sigmas_hat_j = tca.mdl$sigmas_hat[1,]
# gammas_hat_j = tca.mdl$gammas_hat[1,]
# 
# W_prime <- replicate(dim(X)[2], matrix(0,k,k), simplify=F)
# for (i in 1:n){
#   W_prime[[i]] = tcrossprod(W[i,],W[i,])/(tau_hat**2)
# }
# C2_prime <- (X - tcrossprod(tca.mdl$C2, tca.mdl$deltas_hat))/(tca.mdl$tau_hat**2)

#### CPP estimate_Z_j_cpp ####
estimate_Z_cpp <- function(X, W, mus_hat, sigmas_hat, tau_hat, C2, deltas_hat, C1, gammas_hat, scale, parallel, num_cores){
  # init
  n <- dim(X)[1]
  m <- dim(X)[2]
  k <- ncol(W)
  Z_hat <- list()
  for (h in 1:k){
    Z_hat[[h]] <- matrix(0,n,m)
    colnames(Z_hat[[h]]) <- colnames(X)
    rownames(Z_hat[[h]]) <- rownames(X)
  }
  #browser()
  # Calculate quantities that can be calculated only once
  W_prime <- replicate(n, matrix(0,k,k), simplify=F)
  for (i in 1:n){
    W_prime[[i]] = tcrossprod(W[i,],W[i,])/(tau_hat**2)
  }
  #if (m==1) deltas_hat <- t(as.matrix(deltas_hat))
  C2_prime <- (X - tcrossprod(C2,deltas_hat))/(tau_hat**2)
  cl <- if (parallel) init_cluster(num_cores) else NULL
  if (parallel) clusterExport(cl, c("W","mus_hat","sigmas_hat","tau_hat","C1","gammas_hat","W_prime","C2_prime","estimate_Z_j"), envir=environment())
  # Estimate Z
  
  res <- pbmclapply(ignore.interactive = T, 1:m, mc.cores = cores, function(j){estimate_Z_j_cpp(W, mus_hat[j,], sigmas_hat[j,], tau_hat, C1, gammas_hat[j,], W_prime, C2_prime[,j], scale)} )
  
  if (parallel) stop_cluster(cl)
  for (j in 1:m){
    for (h in 1:k){
      Z_hat[[h]][,j] = res[[j]][,h]
    }
  }

  # add rownames and colnames and transpose matrices

  for (h in 1:k){
    rownames(Z_hat[[h]]) <- rownames(X)
    colnames(Z_hat[[h]]) <- colnames(X)
    Z_hat[[h]] <- t(Z_hat[[h]])

  }
  
  return(Z_hat)
}

#### tensor_cpp ####
assert <- function (expr, error) {
  if (! expr) stop(error, call. = FALSE)
}

tensor.validate_input <- function(X, scale, parallel, num_cores, log_file, debug){
  assert(is.matrix(X), "X must be of class 'matrix'")
  assert(!is.null(rownames(X)) & !is.null(colnames(X)), "X must have row names and column names")
  assert(is.logical(scale), "scale must be of class 'logical'")
  assert(is.logical(debug), "debug must be of class 'logical'")
  assert(is.logical(parallel), "parallel must be of class 'logical'")
  assert(is.null(num_cores) | is.numeric(num_cores), "argument num_cores must take a numric value or NULL")
  assert(is.character(log_file) | is.null(log_file), "log_file must be of class 'character' or NULL")
  
}


tensor_cpp <- function(X, tca.mdl, scale = FALSE, parallel = FALSE, num_cores = NULL, log_file = "TCA.log", debug = FALSE, verbose = TRUE){
  
  config <- config::get(file = system.file("extdata", "config.yml", package = "TCA"), use_parent = FALSE)
  
  tensor.validate_input(X, scale, parallel, num_cores, log_file, debug)
  
  # progress bar options
#  op <- if (verbose) pboptions(nout = config$nout, type = config[["type"]]) else pboptions(nout = config$nout, type = "none")
  
  flog.info("Starting tensor for estimating Z...")
  
  X <- if (is.matrix(X)) X else as.matrix(X)
  
  W <- tca.mdl[["W"]]
  mus_hat <- tca.mdl[["mus_hat"]]
  sigmas_hat <- tca.mdl[["sigmas_hat"]]
  tau_hat <- tca.mdl[["tau_hat"]]
  deltas_hat <- tca.mdl[["deltas_hat"]]
  gammas_hat <- tca.mdl[["gammas_hat"]]
  C1 <- tca.mdl[["C1"]]
  C2 <- tca.mdl[["C2"]]

  return( estimate_Z_cpp(t(X), W, mus_hat, sigmas_hat, tau_hat, C2, deltas_hat, C1, gammas_hat, scale, parallel, num_cores) )
}

#### run tensor ####

#input_X <- fread(args[[1]], nThread = 4, verbose = T)
input_X <- vroom::vroom(args[[1]], delim = ',', altrep = T, num_threads = 2, progress = T ) %>% 
    column_to_rownames('probe') %>% 
    as.matrix

# if (args[[3]])
# tmp_out <- epidish(tmp, ref.m = centDHSbloodDMC.m, method = 'RPC')

epi_W <- vroom::vroom(args[[2]]) %>% 
    column_to_rownames('probe') %>% 
    as.matrix

gc()

tca.mdl <- tca_pbmc(input_X, epi_W, constrain_mu = T, debug = T)
write_rds(tca.mdl, file = paste0(args[[2]],'.rds'))
tensor_res <- tensor_cpp(input_X, tca.mdl)
names(tensor_res) <- colnames(tca.mdl$W)
for(cells in colnames(tca.mdl$W)){
  print (cells)
  tensor_res[[cells]] %>% as.data.table(keep.rownames = T) %>% setnames('rn','probe') %>% fwrite(., paste0(args[[3]], '_', cells,'.txt'), sep = '\t')
}













