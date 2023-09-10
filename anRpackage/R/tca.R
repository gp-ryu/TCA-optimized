
# tca_cpp ####

#' @importFrom pracma repmat
#' @importFrom pracma lsqlincon
#' @importFrom nloptr nloptr
#' @importFrom matrixStats colVars

minus_log_likelihood_tau <- function(U,W_squared,sigmas,const,tau, isll = FALSE){
  m <- ncol(U)
  # res <- matrix(0,m,2)
  V <- innerprodMM(W_squared, sigmas^2) + tau^2
  # sumU <- sum(U)
  # sumV <- sum(V)
  # UbyV <- U/V
  if(isll) return((-.5*(const*m - sumlog_p(U) - sum(U/V))))
  return(list( "objective" = sum((-.5*(const*m - sumlog_p(V) - sum(U/V)))), "gradient" = sum((-tau * (sum(U/V^2) - sum(1/V)))) ))
 
  { 
  # t2 <- (-tau * (colSums(UbyVbyV) - colSums(onebyV)))
  # tmp <- pbmclapply(ignore.interactive = T, 1:m, mc.cores = cores, function(j){
  #   V_tmp <- V[,j]
  #   #V_squared <- V**2;
  #   return (c(-0.5*(const-sum(log(V_tmp))-sum(U[,j]/V_tmp)), -(tau*(sum(U[,j]/V_tmp^2)-sum(1./V_tmp))))) 
  #   })
  # res <- tmp %>% unlist %>% matrix(., nrow = 2) %>% t %>% colSums()
  # for (j in 1:m){
  #   res[j,] = tmp[[j]]
  # }
  # res <- colSums(res)}
  }
}
  

minus_log_likelihood_sigmas <- function(sigmas,U_j,W_squared,const,tau){
  k <- length(sigmas)
  n <- nrow(W_squared)
  #V <- tcrossprod(W_squared,t(sigmas**2))+tau**2
  V <- (W_squared %*% sigmas**2) + tau**2
  V_squared <- V**2
  
  t1 <- (-.5*(const - sum(log(V) - sum(U_j/V))))
  t2_1 <- crossprod(W_squared, ((U_j / V_squared) %*% sigmas)) %>% Diag
  t2_2 <- crossprod(W_squared, ((1./V) %*% sigmas)) %>% Diag
  t2 <- -(t2_1 - t2_2)
  return(list("objective" = t1, "gradient" = t2))
  
  # return(list("objective"= -0.5*(const-sum(log(V))-sum(U_j/V)),
  #             "gradient" = -(colSums(W_squared*repmat(sigmas,n,1)*t(repmat(U_j,k,1))/repmat(V_squared,1,k)) - colSums(W_squared*repmat(sigmas,n,1)/repmat(V,1,k)))))
}


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
  W_squared <- W^2
  W_squared_ <- cbind(W_squared,numeric(n) + 1)
  
  ll_prev <- -Inf
  # Perform an alternative optimization of the means (mus, deltas, gammas) and variances (sigmas and tau)
  for (iter in 1:max_iters){
    
    flog.info("Iteration %s out of %s internal iterations...", iter, max_iters)
    
    # (1) Estimate the µ, ∂, γ--------------------
    
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
    
    if (sum(mus_hat) == 0){
      # Use the following for getting initial estimates of mus, deltas, gammas; under the assumptions that tau=0 and sigmas_{1j},...,sigmas_{kj} for each j.
      flog.debug("Calculate W_norms")
      W_norms <- sqrt(rowSums(W**2))
      # Since W_norms is the same for all features in this case can already calculate the following quantities
      W_tilde <- W/W_norms
      C1_tilde <- if (p1>0) C1_/W_norms else C1_
      C2_tilde <- if (p2>0) C2/W_norms else C2
      
      flog.debug("Estimate mus, deltas, gammas.")
      Dmat_tmp <- AtA(cbind(W_tilde,C2_tilde,C1_tilde))
      dvec_tmp <- crossprod(cbind(W_tilde,C2_tilde,C1_tilde), X/W_norms)
      Amat_tmp <- rbind( diag(length(lb)), -diag(length(ub)))
      bvec_tmp <- c(lb, -ub)
    
      res <- pbmclapply(ignore.interactive = T, 1:m, mc.cores = cores, function(j){
        quadprog::solve.QP(Dmat = Dmat_tmp, 
                           dvec = dvec_tmp[,j] %>% matrix(ncol = 1),
                           Amat = t(Amat_tmp),
                           bvec = bvec_tmp,
                           meq = 0
        )$solution
      })
      rm(Dmat_tmp, dvec_tmp, Amat_tmp, bvec_tmp, W_norms)
      gc()
    }else{
      flog.debug("Calculate W_norms")
      #W_norms <- (tcrossprod(W**2,sigmas_hat**2)+tau_hat**2 )**0.5
      W_norms <- eigenW_norms(W^2 ,t(sigmas_hat^2), tau_hat, 1)
      X_tilde <- X/W_norms
      flog.debug("Estimate mus, deltas, gammas.")
      if(p1 == 0 || p2 == 0) conds = T
      res <- pbmclapply(ignore.interactive = T, 1:m,mc.cores = cores, function(j){
        if(conds){
          lsqlincon_opt(W/W_norms[,j], X_tilde[,j], lb = lb, ub = ub)
        }else{
          pracma::lsqlincon_opt(cbind(W/W_norms[,j],
                                  if (p2>0) C2/t(repmat(W_norms[,j],p2,1)) else C2,
                                  if (p1>0) C1_/t(repmat(W_norms[,j],k*p1,1)) else C1_),
                            X_tilde[,j], lb = lb,  ub = ub)
        }
      })
      rm(W_norms, X_tilde)
      gc()
    }
    # X_tilde <- X/W_norms
    
    
    if (sum(colSums(mus_hat)) == 0){
      #if (parallel) clusterExport(cl, varlist = c("W_tilde","C2_tilde","C1_tilde","X_tilde","lb","ub","k","p1","p2","lsqlincon"), envir=environment())
      # res <- pbmclapply(ignore.interactive = T, 1:m, mc.cores = cores, function(j) {
      #   lsqlincon(cbind(W_tilde,C2_tilde,C1_tilde), 
      #             X_tilde[,j], 
      #             lb = lb,  
      #             ub = ub)})
    }else{
      #if (parallel) clusterExport(cl, c("W_norms","W","C2","C1_","X_tilde","lb","ub","k","p1","p2","lsqlincon"), envir=environment())
      # if(p1==0 || p2==0){
      #   x_tmp <- crossprod(W,W)
      #   
      # }else if(p2>0){
      #   x_tmp <- crossprod()
      # }
      
    }
    
    # (1.1) Update estimtes --------------------
    flog.debug('Update estimates mus, deltas, gammas')
    
    if(p1==0 || p2==0){
      mus_hat <- res %>% unlist %>% matrix(., nrow = k) %>% t
      rm(res)
    }else{
      tmp_res <- res %>% unlist %>% matrix(., nrow = k + p2 + p1*k) %>% t
      mus_hat <- tmp_res[,1:k]
      deltas_hat <- tmp_res[,k+1:k+p2]
      gammas_hat <- tmp_res[,k+p2+1,k+p2+p1*k]
      rm(res)
    }
    
    # pbmclapply(1:m, mc.cores = cores, function(j){
    #   mus_hat[j,] <<- res[[j]][1:k]
    #   deltas_hat[j,seq(1,p2,length=p2)] = res[[j]][seq(k+1,k+p2,length=p2)]
    #   gammas_hat[j,seq(1,k*p1,length=k*p1)] = res[[j]][seq(k+p2+1,k+p2+p1*k,length=p1*k)]
    # })
    # for (j in 1:m){
    #   mus_hat[j,] = res[[j]][1:k]
    #   deltas_hat[j,seq(1,p2,length=p2)] = res[[j]][seq(k+1,k+p2,length=p2)]
    #   gammas_hat[j,seq(1,k*p1,length=k*p1)] = res[[j]][seq(k+p2+1,k+p2+p1*k,length=p1*k)]
    # }
    
    
    # (2) Estimate σ, τ --------------------
    
    # Calculate some quantities that will be repeatedly used throughout the optimization in this step
    if(p1==0 || p2==0) {#U <- ((eigenProduct(W,t(mus_hat)) - X)^2) 
    }else{
      # U <- (eigenProduct(W,mus_hat %>% t) + tcrossprod(C2,deltas_hat) + tcrossprod(C1_,gammas_hat) - X)**2
    }
    if (vars.mle){
      
      U <- (eigenProduct(W,mus_hat %>% t) + tcrossprod(C2,deltas_hat) + tcrossprod(C1_,gammas_hat) - X)**2
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
      res <- pbmclapply(ignore.interactive = T, 1:m,mc.cores = cores, function(j){
        nloptr( x0=sigmas_hat[j,],
                eval_f = function(x,U_j,W_squared,const,tau_hat) minus_log_likelihood_sigmas(x,U_j,W_squared,const,tau_hat),
                lb = lb,
                ub = ub,
                opts = nloptr_opts,
                U_j = U[,j],
                W_squared = W_squared,
                const = const,
                tau_hat = tau_hat)$solution
      })
      
      flog.info('update estimates')
      sigmas_hat <- res %>% unlist %>% matrix(., nrow = k) %>% t 
      
      # for (j in 1:m){
      #   sigmas_hat[j,] = res[[j]]
      # }
      
      # (2.3) Estimate tau
      flog.info('estimate tau')
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
      browser()
      # calculate weights matrix V
      conds   <- (sum(sigmas_hat) == 0)
      if (conds){
        V <- matrix(1,n,m)
        # norms_tmp <- W_squared_^2 %>% colSums %>% repmat(., m,1) %>% t %>% eigenSqrt()
        norms_tmp <- colNorm(W_squared_)
        # norms_tmp <- sqrt(colSums(W_squared_^2))
        x_tmp <- W_squared_
        Dmat_tmp <- AtA(W_squared_/norms_tmp)
        dvec_tmp <- crossprod(W_squared_/(norms_tmp), (tcrossprod(W,mus_hat ) + tcrossprod(C2,deltas_hat) + tcrossprod(C1_,gammas_hat) - X)**2)
      }else{
        V <- abs((tcrossprod(W,mus_hat ) + tcrossprod(C2,deltas_hat) + tcrossprod(C1_,gammas_hat) - X)**2 - tcrossprod(W_squared_, cbind(sigmas_hat**2, matrix(tau_hat**2,m,1))))
        V[V < config[["V_weight_limit"]]] <- config[["V_weight_limit"]]
        
        d_tmp <- ((tcrossprod(W,mus_hat ) + tcrossprod(C2,deltas_hat) + tcrossprod(C1_,gammas_hat) - X)**2)/V
        x_tmp <- W_squared_
      }
      # if (parallel) clusterExport(cl, c("lb","U","W_squared_","lsqlincon","V","n"), envir=environment())
      
      Aeq     = NULL
      A       = NULL
      diag_ub = NULL
      ub      = NULL
      beq     = NULL
      b       = NULL
      res <- pbmclapply(ignore.interactive = T, 1:m, mc.cores = cores, function(j) {
        if (conds){
          dvec_tmp <- dvec_tmp[,j]
          norms_tmptmp <- norms_tmp
        }else{
          x_tmptmp <-  x_tmp/V[,j]
          # norms_tmptmp <-colNorm(x_tmptmp); names(norms_tmptmp) <- colnames(x_tmptmp)
          norms_tmptmp <-sqrt(colSums(x_tmptmp^2))
          x_tmptmp <- x_tmptmp/norms_tmptmp 
          
          Dmat_tmp <- crossprod(x_tmptmp)
          dvec_tmp <- crossprod(x_tmptmp, d_tmp[,j]) %>% .[,1]
        }
        rslt <- quadprog::solve.QP(Dmat_tmp,
                                   dvec_tmp,
                                   rbind(Aeq, A, diag(length(lb)), diag_ub),
                                   c(beq, b, lb*norms_tmptmp, ub),
                                   meq = 0)$solution
         rslt / norms_tmptmp
      }
        )
      # res <- pbmclapply(ignore.interactive = T, 1:m,mc.cores = cores, function(j) {
      #   if (!conds){
      #     dvec_tmp <- dvec_tmp[,j]
      #     bvec <- c(beq, b, lb*norms_tmp, ub)
      #   }else{
      #     # x_tmp <- ( x_tmp/matrix(unlist(V[, j])%>%rep(., 7),ncol=ncol(W_squared_)) )
      #     x_tmp <- ( x_tmp%*%diag(1/V[,j]) )
      #     norms_tmp <- sqrt(colSums(x_tmp**2))
      #     x_tmp <- ( x_tmp/t(matrix(norms%>%rep(.,nrow(x_tmp)),ncol=nrow(x_tmp))) )
      #     
      #     Dmat_tmp <- crossprod(x_tmp)
      #     dvec_tmp <- crossprod(x_tmp, d_tmp[,j])
      #     bvec <- c(beq, b, lb*norms, ub)
      #   }
      #   
      #   diag_lb <- diag(length(lb))
      #   Amat <- rbind(Aeq, A, diag_lb, diag_ub)
      #   rslt <- quadprogpp::QP.Solve(Dmat_tmp, 
      #                              -dvec_tmp, 
      #                              Amat,  
      #                              -bvec)
      #   rslt / norms_tmp
      # })
      rm(norms_tmp, d_tmp, x_tmp, V) %>% suppressWarnings
      gc()
      
      # res <- pbmclapply(ignore.interactive = T, 1:m,mc.cores = cores, function(j) {
      #   x <- W_squared_/V[,j]
      #   # For numeric stability, normalize the design matrix and adjust the final estimats accordingly
      #   norms <- (colSums(x**2))**0.5;
      #   x <- x%*%diag(1/norms);
      #   lsqlincon(x, U[,j]/V[,j], lb = lb*norms)/norms
      # })
      
      
      ## (2.1) Update estimates --------------------
      browser()
      tmp_res <- res %>% unlist %>% matrix(., nrow = k+1) %>% t 
      sigmas_hat <- sqrt(tmp_res[,1:k]) 
      tau_hat <- sqrt(mean(tmp_res[,k+1]))
      rm(res, tmp_res)
      gc()
      
    }
    
    # (3) Test for convergence --------------------
    flog.debug("Test for convergence.")
    ll_new <- -minus_log_likelihood_tau((tcrossprod(W,mus_hat ) + tcrossprod(C2,deltas_hat) + tcrossprod(C1_,gammas_hat) - X)**2,W_squared,sigmas_hat,const,tau_hat, isll = T) 
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
    U <- (eigenProduct(W,mus_hat %>% t) + tcrossprod(C2,deltas_hat) + tcrossprod(C1_,gammas_hat) - X)**2
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

lsqlincon_opt <- function (C, d, A = NULL, b = NULL, Aeq = NULL, beq = NULL, lb = NULL, 
                           ub = NULL) 
{
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
  Dmat <- AtA(C) # t(C) %*% C
  dvec <- outerprodMV_p(C, d) # t(C) %*% d
  Amat <- rbind(Aeq, A, diag_lb, diag_ub)
  bvec <- c(beq, b, lb, ub)
  rslt <- quadprog::solve.QP(Dmat, dvec, t(Amat), bvec, meq = meq)
  rslt$solution
}
