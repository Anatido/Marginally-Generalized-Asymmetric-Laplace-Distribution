# k-dimensional Gibbs/MH sampler for latent weights W (n x k)
# - Preserves your 2D logic: GIG proposals per coordinate; acceptance uses full quadratic form.
# - Sigma_val is k x k covariance; alpha_val, eta_val, mu_val are length-k.
# - Returns an n x k matrix of positive weights.
sample_w_function_k <- function(
    Y,
    alpha_val,
    Wdata,
    eta_val,
    mu_val,
    Sigma_val
){
  Y <- as.matrix(Y)
  n <- nrow(Y); k <- ncol(Y)
  
  stopifnot(length(alpha_val) == k,
            length(eta_val)   == k,
            length(mu_val)    == k,
            is.matrix(Sigma_val),
            all(dim(Sigma_val) == c(k,k)))
  
  # --- numeric guards/helpers
  eps <- 1e-12
  posify <- function(x) pmax(as.numeric(x), eps)
  finfix <- function(x, alt = 1) { x[!is.finite(x)] <- alt; x }
  
  # --- standardize by Sigma's diagonal (to mirror your gamma constructs)
  sig_diag <- diag(Sigma_val)
  sd_vec   <- sqrt(posify(sig_diag))       # length-k, sqrt(σ_jj)
  # standardized "gammas"
  g_eta <- (alpha_val * eta_val) / sd_vec  # α_j η_j / √σ_jj   (your gamma_1, gamma_3, …)
  g_y   <- (Y - matrix(mu_val, n, k, byrow = TRUE)) /
    matrix(sd_vec, n, k, byrow = TRUE)       # (y_j - μ_j)/√σ_jj (your gamma_2, gamma_4, …)
  
  # correlation + precision in standardized basis
  Dm12   <- diag(1/sd_vec, k)
  R      <- Dm12 %*% Sigma_val %*% Dm12       # correlation matrix
  Rinvt  <- tryCatch(solve(R), error = function(e) NULL)
  if (is.null(Rinvt)) {
    # tiny jitter if needed
    eig <- eigen(R, symmetric = TRUE)
    eig$values <- pmax(eig$values, eps)
    R <- eig$vectors %*% diag(eig$values, k) %*% t(eig$vectors)
    Rinvt <- solve(R)
  }
  Rinvt <- (Rinvt + t(Rinvt))/2  # symmetrize
  
  # per-coordinate GIG "shape" lambda = 1/alpha - 1/2  (your p1, p2)
  p_vec <- 1/alpha_val - 1/2
  
  # convenience closures for proposal params (match your 2D a1,b1,a2,b2 forms):
  # chi_j = Rinvt[j,j] * (g_y_ij)^2
  # psi_j = Rinvt[j,j] * (g_eta_j)^2 + 2/alpha_j
  proposal_params <- function(j, gi_y_row) {
    chi <- Rinvt[j,j] * (gi_y_row[j]^2)
    psi <- Rinvt[j,j] * (g_eta[j]^2) + 2/alpha_val[j]
    list(chi = chi, psi = psi, lambda = p_vec[j])
  }
  
  # log-target for a single coordinate w_j given others (row-wise)
  # This encodes: prior ∝ w_j^(1/α_j - 3/2), and the quadratic form
  #   t = A√w - B/√w,  with A = diag(g_eta), B = diag(g_y_row)
  #   loglik ∝ -0.5 * t^T Rinvt t  - 0.5 * (2/α_j) w_j   (the latter sits in psi above; here it’s in full t^T term)
  log_target_wj <- function(wj, j, w_row, gi_y_row) {
    # build √w and 1/√w with current row; replace the jth entry
    sw  <- sqrt(posify(w_row))
    isw <- 1 / sw
    sw[j]  <- sqrt(posify(wj))
    isw[j] <- 1 / sw[j]
    
    # A√w - B/√w
    tvec <- g_eta * sw - gi_y_row * isw
    
    # quadratic term
    qf <- as.numeric(t(tvec) %*% Rinvt %*% tvec)  # t^T Rinvt t
    
    # prior-ish log term consistent with your acceptance algebra:
    # (1/α_j - 1/2 - 1)*log w_j = (1/α_j - 3/2) * log w_j
    prior_log <- (1/alpha_val[j] - 1.5) * log(posify(wj))
    
    # total (up to additive constants that cancel in MH ratio)
    -(0.5)*qf + prior_log
  }
  
  # initialize outputs
  out <- matrix(NA_real_, n, k)
  colnames(out) <- paste0("w_", seq_len(k))
  
  # --- main loop over observations
  for (i in seq_len(n)) {
    # starting point (row i)
    w_row <- as.numeric(Wdata[i, ])
    if (length(w_row) != k) {
      stop("Wdata must be n x k with initial weights per observation.")
    }
    w_row <- finfix(posify(w_row), 1)
    
    gi_y_row <- as.numeric(g_y[i, ])
    
    # coordinate-wise MH over j = 1..k
    for (j in seq_len(k)) {
      # current value
      wj_curr <- w_row[j]
      
      # proposal ~ GIG(lambda=p_j, chi=Rinv[j,j]*(g_y_ij)^2, psi=Rinv[j,j]*(g_eta_j)^2 + 2/α_j)
      pp <- proposal_params(j, gi_y_row)
      wj_prop <- GIGrvg::rgig(1, lambda = pp$lambda, chi = posify(pp$chi), psi = posify(pp$psi))
      
      # log target difference
      lt_prop <- log_target_wj(wj_prop, j, w_row, gi_y_row)
      lt_curr <- log_target_wj(wj_curr, j, w_row, gi_y_row)
      
      # proposal density correction (independence MH)
      # add log q(curr) - log q(prop)
      lq_curr <- GIGrvg::dgig(wj_curr, lambda = pp$lambda, chi = posify(pp$chi), psi = posify(pp$psi), log = TRUE)
      lq_prop <- GIGrvg::dgig(wj_prop, lambda = pp$lambda, chi = posify(pp$chi), psi = posify(pp$psi), log = TRUE)
      
      log_alpha <- (lt_prop - lt_curr) + (lq_curr - lq_prop)
      if (is.finite(log_alpha) && log(runif(1)) < min(log_alpha, 0)) {
        w_row[j] <- wj_prop  # accept
      } # else reject and keep current
    }
    
    out[i, ] <- w_row
  }
  
  out
}
