##
covariance_estimation_k <- function(
    data_Y,
    alpha_r,   # length k
    mu_r,      # length k
    eta_r,     # length k
    Sigma_r,   # k x k, used by the W-sampler; updated each EM iter
    burn_in = 10,
    num_iterations = 100,
    M = 10,          # total draws per E-step (including burn-in)
    sampler = sample_w_function_k,  # function(Y, alpha_val, Wdata, eta_val, mu_val, Sigma_val)
    psd_jitter = 1e-10              # tiny PSD guard
){
  Y <- as.matrix(data_Y)
  n <- nrow(Y); k <- ncol(Y)
  
  stopifnot(M > burn_in,
            length(alpha_r) == k,
            length(mu_r)    == k,
            length(eta_r)   == k,
            is.matrix(Sigma_r),
            all(dim(Sigma_r) == c(k,k)))
  
  eps <- 1e-12
  pos <- function(x) pmax(as.numeric(x), eps)
  
  alpha_eta <- alpha_r * eta_r  # length k, reused
  
  # initialize W seeds
  Wdata0 <- matrix(1, n, k)
  
  for (iter in seq_len(num_iterations)) {
    
    # accumulate crossproducts of (A - V)
    S_sum <- matrix(0, k, k)
    Meff  <- 0L
    
    for (m in seq_len(M)) {
      Wm <- sampler(
        Y        = Y,
        alpha_val = alpha_r,
        Wdata     = Wdata0,
        eta_val   = eta_r,
        mu_val    = mu_r,
        Sigma_val = Sigma_r
      )
      Wm <- as.matrix(Wm)
      # guard
      if (!all(dim(Wm) == c(n,k))) stop("sampler must return n x k weights.")
      Wm <- pmax(Wm, eps)
      
      if (m > burn_in) {
        sw  <- sqrt(Wm)           # n x k
        isw <- 1 / sw             # n x k
        
        # A = (Y - mu) ∘ (1/√W), V = (alpha*eta) ∘ √W
        Yc  <- sweep(Y, 2, mu_r, "-")       # n x k
        A   <- Yc * isw                     # n x k
        V   <- sweep(sw, 2, alpha_eta, "*") # n x k
        
        Z   <- A - V                        # n x k
        # BLAS-3: k x n  %*% n x k
        S_sum <- S_sum + crossprod(Z)       # (A-V)^T (A-V)
        Meff  <- Meff + 1L
        
        # carry forward W for next proposal state (mixes chain faster)
        Wdata0 <- Wm
      } else {
        Wdata0 <- Wm
      }
    }
    
    # average over effective draws and normalize by n
    Sigma_new <- (S_sum / max(Meff,1L)) / n
    
    # symmetrize + tiny PSD guard
    Sigma_new <- 0.5 * (Sigma_new + t(Sigma_new))
    eig <- eigen(Sigma_new, symmetric = TRUE)
    if (min(eig$values) < psd_jitter) {
      eig$values <- pmax(eig$values, psd_jitter)
      Sigma_new <- eig$vectors %*% diag(eig$values, k) %*% t(eig$vectors)
    }
    
    Sigma_r <- Sigma_new
  }
  
  list(Sigma_mle = Sigma_r, iter = num_iterations)
}
