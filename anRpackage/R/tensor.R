# CPP estimate_Z_j_cpp ####
estimate_Z_cpp <- function(X, W, mus_hat, sigmas_hat, tau_hat, C2, deltas_hat, C1, gammas_hat, scale, cohort = "TCA"){
  # init
  n <- dim(X)[1]
  m <- dim(X)[2]
  k <- ncol(W)
  # Z_hat <- list()
  # for (h in 1:k){
  #   Z_hat[[h]] <- matrix(0,m,n)
  #   colnames(Z_hat[[h]]) <- rownames(X)
  #   rownames(Z_hat[[h]]) <- colnames(X)
  # }
  # Calculate quantities that can be calculated only once
  flog.debug('W_prime calc')
  W_prime <- W/tau_hat
  
  # W_prime vectorize
  # W_prime <- replicate(n, matrix(0,k,k), simplify=F)
  # for (i in 1:n){
  #   W_prime[[i]] = tcrossprod(W[i,])/(tau_hat**2)
  # }
  
  
  #if (m==1) deltas_hat <- t(as.matrix(deltas_hat))
  if(ncol(C2) != 0){
    C2_prime <- (X - tcrossprod(C2,deltas_hat))/(tau_hat**2)
  }else{
    C2_prime = X/tau_hat^2
    flog.info('Starting estimate_Z_j_cpp ...')
    res <- pbmclapply(ignore.interactive = T, 1:m, mc.cores = cores, function(j){
      estimate_Z_j_cpp(X , W , mus_hat[j,], sigmas_hat[j,]^2, tau_hat, C1, gammas_hat[j,], W_prime, C2_prime[,j], scale) %>% as.data.frame() 
    }) %>% rbindlist 
  } 
  flog.info('Write tensor_res to rds')
  gc()
  write_rds(res, file = paste0(cohort,'_tensor_res.rds'))
  for(h in 1:k){
    flog.info(paste0('Writing ',colnames(W)[h]))
    res[[h]] %>% unlist %>% 
      round(7) %>% 
      matrix(., byrow = T, nrow = m, ncol = n,
             dimnames = list(colnames(X), rownames(X))) %>% as.data.frame %>% rownames_to_column(., 'probe') %>%
      fwrite(., paste0(cohort,'_beta_inner_decompose_', colnames(W)[h], '.csv'), sep = ',')
  }
  
  
  # flog.info("file writing")
  # write_tensor_header(paste0(cohort,'_beta_inner_decompose_',colnames(W),".txt"), 
  #                           rownames(W))
  # write_tensor_list(paste0(cohort, '_beta_inner_decompose_',colnames(W),".txt"), 
  #                       rownames(W), colnames(X), res)
  # write_tensor_list(paste0(cohort, '_beta_inner_decompose_',colnames(W),".txt"), 
                        # rownames(W), colnames(X)[1], matrix(res[[1]], ncol = 7) %>% list)


  # Estimate Z
  # rm(C2_prime, X)
  # gc()
  # add rownames and colnames and transpose matrices
  # for (h in 1:k){
  #   rownames(Z_hat[[h]]) <- rownames(X)
  #   colnames(Z_hat[[h]]) <- colnames(X)
  #   Z_hat[[h]] <- t(Z_hat[[h]])
  # }
  # return(Z_hat)
}





estimate_Z_j_tmp <- function(X, W, mus_hat, sigmas_hat, tau_hat, C2, deltas_hat, C1, gammas_hat, scale){
  # init
  # n <- dim(X)[1]
  # m <- dim(X)[2]
  k <- ncol(W)
  # Z_hat <- list()
  # for (h in 1:k){
  #   Z_hat[[h]] <- matrix(0,n,m)
  #   colnames(Z_hat[[h]]) <- colnames(X)
  #   rownames(Z_hat[[h]]) <- rownames(X)
  # }
  
  flog.info('Starting estimate_Z_j_standalone ...')
  Z_hat <- estimate_Z_j_standalone(X, W, mus_hat, sigmas_hat, tau_hat, scale)
  
  # flog.debug('W_prime calc')
  # W_prime <- replicate(n, matrix(0,k,k), simplify=F)
  # W_prime <- (W/tau_hat)^2
  # for (i in 1:n){
  #   W_prime[[i]] = tcrossprod(W[i,])/(tau_hat**2)
  # }
  #if (m==1) deltas_hat <- t(as.matrix(deltas_hat))
  # if(ncol(C2) != 0){
  #   C2_prime <- (X - tcrossprod(C2,deltas_hat))/(tau_hat**2)
  # }else{
  #   C2_prime <- X
  # }
  
  # 
  # 
  # for (i in 1:n) {
  #   BA_inv <- W_prime[[i]] %*% Sig_j_orig
  #   g <- sum(diag(BA_inv))
  #   S_ij <- Sig_j_orig - (1 / (1 + g)) * (Sig_j_orig %*% BA_inv)
  #   C1_prime_j <- C1 %*% diag(gammas_hat[, j])
  #   
  #   S_ij_tmptmp <- crossprod(W[i, ], C2_prime[, j])
  #   S <- t(S_ij) %*% (solve(Sig_j_orig) %*% (mus_hat + C1_prime_j) + S_ij_tmptmp)
  #   
  #   Z_j_hat[i, ] <- S
  # }
  # 
  
  
  # res <- pbmclapply(ignore.interactive = T, 1:m, mc.cores = cores, function(j){
  #   # estimate_Z_j_cpp(W, mus_hat[j,], sigmas_hat[j,]^2, tau_hat, C1, gammas_hat[j,], W_prime, C2_prime[,j], scale)} )
  #   estimate_Z_j_tmp(W, mus_hat[j,], sigmas_hat[j,]^2, tau_hat, C1, gammas_hat[j,], W_prime, C2_prime[,j], scale)} )
  
  gc()
  # if (parallel) stop_cluster(cl)
  # tmp <- lapply(res, function(matrix) matrix[, h]) 
  # flog.info('list dim change')
  # for (j in 1:m){
  #   for (h in 1:k){
  #     Z_hat[[h]][,j] = res[[j]][,h]
  #   }
  # }
  # add rownames and colnames and transpose matrices
  
  for (h in 1:k){
    rownames(Z_hat[[h]]) <- rownames(X)
    colnames(Z_hat[[h]]) <- colnames(X)
    Z_hat[[h]] <- t(Z_hat[[h]])
  }
  return(Z_hat)
}





# tensor_cpp ####
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


tensor_cpp <- function(X, tca.mdl, scale = FALSE, parallel = FALSE, num_cores = NULL, log_file = "TCA.log", debug = FALSE, verbose = TRUE, cohort = args[[3]]){
  
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
  
  return( estimate_Z_cpp(t(X), W, mus_hat, sigmas_hat, tau_hat, C2, deltas_hat, C1, gammas_hat, scale, cohort) )
}

#
tensor_orgn <- function(X, tca.mdl, scale = FALSE, parallel = FALSE, num_cores = NULL, log_file = "TCA.log", debug = FALSE, verbose = TRUE){

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

  return( estimate_Z_orgn(t(X), W, mus_hat, sigmas_hat, tau_hat, scale) )
}



tensor_tmp <- function(X, tca.mdl, scale = FALSE, parallel = FALSE, num_cores = NULL, log_file = "TCA.log", debug = FALSE, verbose = TRUE){

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

  return( estimate_Z_j_standalone(t(X), W, mus_hat, sigmas_hat, tau_hat, scale) )
}

estimate_Z_orgn <- function(X, W, mus_hat, sigmas_hat, tau_hat, C2, deltas_hat, C1, gammas_hat, scale, parallel, num_cores){
  flog.info("Estimate tensor...")
  # init
  n <- nrow(X)
  m <- ncol(X)
  k <- ncol(W)
  scale = F
  Z_hat <- list()
  for (h in 1:k){
    Z_hat[[h]] <- matrix(0,n,m)
    colnames(Z_hat[[h]]) <- colnames(X)
    rownames(Z_hat[[h]]) <- rownames(X)
  }
  # Calculate quantities that can be calculated only once
  W_prime <- replicate(n, matrix(0,k,k), simplify=F)
  for (i in 1:n){
    W_prime[[i]] = tcrossprod(W[i,],W[i,])/(tau_hat**2)
  }
  #if (m==1) deltas_hat <- t(as.matrix(deltas_hat))

  C2_prime <- (X )/(tau_hat**2)
  # cl <- if (parallel) init_cluster(num_cores) else NULL
  # if (parallel) clusterExport(cl, c("W","mus_hat","sigmas_hat","tau_hat","C1","gammas_hat","W_prime","C2_prime","estimate_Z_j"), envir=environment())
  # Estimate Z
  browser()
  res <- pbmclapply(1:1, function(j)
    estimate_Z_j_orgn(W, mus_hat[j,], sigmas_hat[j,], tau_hat, gammas_hat[j,], W_prime, C2_prime[,j], scale)
    )
  # if (parallel) stop_cluster(cl)
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
  flog.info("Finished estimating tensor.")
  return(Z_hat)
}

estimate_Z_j_orgn <- function(W, mus_hat_j, sigmas_hat_j, tau_hat,gammas_hat_j, W_prime, C2_prime_j, scale){
  n <- nrow(W)
  k <- ncol(W)
  # p1 <- ncol(C1)
  Z_j_hat <- matrix(0,n,k)
  Sig_j_orig <- diag(sigmas_hat_j**2)
  Sig_j <- matrix.inverse(Sig_j_orig)
  # C1_prime <- tcrossprod(C1, t(Reshape(gammas_hat_j,p1,k)))

  browser()
  # for (i in 1:1){
  i = 1
    #S_ij_inv <- W_prime[[i]]+Sig_j
    #S_ij <- matrix.inverse(S_ij_inv)
    ## the above two lines are the straightforward (slower) way to calculate S_ij; calculate 'matrix.inverse(W_prime[[i]]+Sig_j)' using the lemma from Miller 1981:
    BA_inv <- W_prime[[i]] %*% Sig_j_orig
    g <- sum(diag(BA_inv))
    S_ij <- Sig_j_orig - ((1/(1+g))*(Sig_j_orig %*% BA_inv))
    Z_j_hat[i,] = crossprod(S_ij, ( tcrossprod(Sig_j,t(mus_hat_j)) + W[i,]*C2_prime_j[i] ) )
    if (scale) Z_j_hat[i,] <- Z_j_hat[i,] / diag(S_ij)**0.5
    print(BA_inv); S_ij; Sig_j_orig; W[i,]*C2_prime_j[i]
  # }
  return(Z_j_hat)
}