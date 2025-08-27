# rMGAL generator for arbitrary dimension k
# Y = mu + diag(W) diag(alpha) eta + diag(W)^{1/2} X,  X ~ N_k(0, Sigma), W_j ~ Gamma(1/alpha_j, 1/alpha_j)
# Returns an n x k matrix; optionally returns the sampled weights W as well.

rMGAL <- function(n,
                  mu_vec,
                  gamma_parameters,   # length k (your "alpha" vector); all > 0
                  eta_vec,            # length k
                  sigma_matrix,       # k x k, symmetric PD
                  seed = NULL,
                  return_weights = FALSE,
                  psd_jitter = 1e-10) {

  if (!is.null(seed)) set.seed(seed)

  # ---- Dimension checks ----
  k <- length(mu_vec)
  if (length(eta_vec) != k)
    stop("eta_vec must have the same length as mu_vec.")
  if (length(gamma_parameters) == 1L) {
    gamma_parameters <- rep(gamma_parameters, k)
  }
  if (length(gamma_parameters) != k)
    stop("gamma_parameters must be length k (or a scalar to recycle).")
  if (any(gamma_parameters <= 0))
    stop("All entries in gamma_parameters must be positive.")
  if (!is.matrix(sigma_matrix) || any(dim(sigma_matrix) != c(k, k)))
    stop("sigma_matrix must be a k x k matrix.")

  # ---- Symmetrize + tiny PSD guard for Sigma ----
  Sigma <- 0.5 * (sigma_matrix + t(sigma_matrix))
  ev <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < psd_jitter) {
    eig <- eigen(Sigma, symmetric = TRUE)
    vals <- pmax(eig$values, psd_jitter)
    Sigma <- eig$vectors %*% diag(vals, k) %*% t(eig$vectors)
  }

  # ---- Draw latent Gaussian X ~ N_k(0, Sigma) ----
  X <- MASS::mvrnorm(n = n, mu = rep(0, k), Sigma = Sigma)  # n x k

  # ---- Draw weights W_j ~ Gamma(shape=1/alpha_j, rate=1/alpha_j) ----
  # mean = 1, var = alpha_j
  W <- sapply(seq_len(k), function(j) rgamma(n, shape = 1 / gamma_parameters[j],
                                             rate  = 1 / gamma_parameters[j]))
  # ensure matrix shape (n x k)
  W <- matrix(W, nrow = n, ncol = k, dimnames = list(NULL, NULL))

  # ---- Build components (all vectorized) ----
  alpha_eta <- gamma_parameters * eta_vec       # length k
  mu_mat    <- matrix(rep(mu_vec, each = n), nrow = n)  # n x k
  term2     <- sweep(W, 2, alpha_eta, `*`)      # n x k : W .* (alpha .* eta)
  term3     <- sqrt(W) * X                      # n x k : diag(W)^{1/2} X (rowwise)

  Y <- mu_mat + term2 + term3                   # n x k

  if (return_weights) {
    return(list(Y = Y, W = W))
  } else {
    return(Y)
  }
}
