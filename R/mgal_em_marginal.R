# mgal_em_marginal.R
# Component-wise (marginal) EM updates for alpha, mu, eta, diag(Sigma)

## Expections for E-steps in Marginal EM's

Ew <- function(chi, psi, lambda) {
  sqrt(chi / psi) * BesselK(sqrt(chi * psi), nu = lambda + 1, expon.scaled = TRUE) /
    BesselK(sqrt(chi * psi), nu = lambda, expon.scaled = TRUE)
}

##
Ew_inv <- function(chi, psi, lambda) {
  sqrt(psi / chi) * BesselK(sqrt(chi * psi), nu = lambda + 1, expon.scaled = TRUE) /
    BesselK(sqrt(chi * psi), nu = lambda, expon.scaled = TRUE) - 2 * lambda /
    chi
}

##
Ew_log <- function(chi, psi, lambda) {
  #if(chi == 0 && lambda > 0){
  #digamma(lambda) - log(0.5*psi)
  #}else{
  
  f <- function(nu, x) {
    BesselK(x, nu, expon.scaled = TRUE)
  }
  nu <- lambda
  h <- 1e-6
  x <- sqrt(chi * psi)
  (log(f(nu + h, x)) - log(f(nu, x))) / h  + log(sqrt(chi / psi))
}


##
parameter_estimate_non_sym_k <- function(
    data_Y,          # n x k matrix/data.frame
    alpha_r,         # length-k vector
    mu_r,            # length-k vector
    eta_r,           # length-k vector
    Sigma_r,         # length-k vector (variances; diagonal Σ)
    num_iterations = 100,
    tolerance = 1e-3 # kept for API compatibility; not used
){
  data_Y   <- as.matrix(data_Y)
  n        <- nrow(data_Y)
  k        <- ncol(data_Y)
  
  stopifnot(length(alpha_r) == k,
            length(mu_r)    == k,
            length(eta_r)   == k,
            length(Sigma_r) == k)
  
  # -- small numerical guardrails
  eps        <- 1e-12
  posify     <- function(x, minv = eps) pmax(as.numeric(x), minv)
  finite_or  <- function(x, alt) ifelse(is.finite(x), x, alt)
  
  # robust alpha solver per component: solve digamma(1/a) + log(a) = target
  find_alpha <- function(target, lower = 1e-6, upper = 1e6) {
    # Reparam: a = exp(t), t in [log(lower), log(upper)] for better stability
    f_t <- function(t) {
      a <- exp(t)
      # digamma(1/a) + log(a) - target
      digamma(1 / a) + t - target
    }
    lo <- log(lower); hi <- log(upper)
    
    # Ensure bracket: expand if needed (limited times)
    f_lo <- f_t(lo); f_hi <- f_t(hi)
    expand_ct <- 0L
    while (sign(f_lo) == sign(f_hi) && expand_ct < 6L) {
      lo <- lo - 2
      hi <- hi + 2
      f_lo <- f_t(lo); f_hi <- f_t(hi)
      expand_ct <- expand_ct + 1L
    }
    
    out <- tryCatch(
      {
        uniroot(f_t, c(lo, hi))$root |> exp()
      },
      error = function(e) {
        # Fallback: clamp to a reasonable value
        # Use exp(target) bounded to [lower, upper] as a crude safe guess
        a_guess <- exp(target)
        a_guess <- min(max(a_guess, lower), upper)
        a_guess
      }
    )
    as.numeric(out)
  }
  
  # vectorized alpha finder
  find_alpha_vec <- function(target_vec) {
    vapply(target_vec, find_alpha, numeric(1))
  }
  
  # main EM loop with fixed iterations (no convergence stopping)
  for (iteration in seq_len(num_iterations)) {
    
    # ---- E STEP ----
    # Per-dim quantities (length k)
    p_vec   <- 1 / alpha_r - 1/2
    sig_sd  <- sqrt(posify(Sigma_r))
    psi_vec <- (alpha_r * eta_r / sig_sd)^2 + 2 / alpha_r
    
    # Pre-allocate expectation matrices (n x k)
    a_mat <- matrix(0, n, k) # E[w]
    b_mat <- matrix(0, n, k) # E[1/w]
    c_mat <- matrix(0, n, k) # E[log w]
    
    # For each observation i, compute chi and expectations for all k dims
    # (kept structurally close to your original inner loop)
    for (i in seq_len(n)) {
      y_i   <- data_Y[i, ]
      chi_i <- ((y_i - mu_r) / sig_sd)^2
      
      # Ew(lambda = p, chi, psi), Ew(lambda = -p, chi = psi, psi = chi), Ew_log(lambda = p, chi, psi)
      a_row <- mapply(function(p, chi, psi) Ew(lambda = p,   chi = chi, psi = psi),
                      p_vec, chi_i, psi_vec)
      b_row <- mapply(function(p, chi, psi) Ew(lambda = -p,  chi = psi, psi = chi),
                      p_vec, chi_i, psi_vec)
      c_row <- mapply(function(p, chi, psi) Ew_log(lambda = p, chi = chi, psi = psi),
                      p_vec, chi_i, psi_vec)
      
      # guard against non-finite values
      a_mat[i, ] <- finite_or(a_row, 1)         # neutral-ish fallback
      b_mat[i, ] <- finite_or(b_row, 1)         # "
      c_mat[i, ] <- finite_or(c_row, log(1))    # "
    }
    
    # ---- M STEP ----
    
    # 1) update alpha (per-dim)
    # mean_vector <- colMeans(c_matrix - a_matrix + 1)
    mean_vec <- colMeans(c_mat - a_mat + 1, na.rm = TRUE)
    alpha_r1 <- find_alpha_vec(mean_vec)
    alpha_r1 <- posify(alpha_r1)  # enforce positivity
    
    # 2) update mu (per-dim)
    # num <- colMeans(a)*colMeans(Y*b) - colMeans(Y)
    # den <- colMeans(a)*colMeans(b) - 1
    ca <- colMeans(a_mat, na.rm = TRUE)
    cb <- colMeans(b_mat, na.rm = TRUE)
    cY <- colMeans(data_Y, na.rm = TRUE)
    cYb <- colMeans(data_Y * b_mat, na.rm = TRUE)
    
    num <- ca * cYb - cY
    den <- ca * cb - 1
    den <- ifelse(abs(den) < eps, sign(den) * eps, den)  # prevent zero denom
    
    mu_r1 <- num / den
    
    # 3) update eta (per-dim) using previous mu (your original logic)
    centered_prev_mu <- sweep(data_Y, 2, mu_r, "-")
    eta_num <- colMeans(centered_prev_mu, na.rm = TRUE)
    eta_den <- ca * alpha_r1
    eta_den <- ifelse(abs(eta_den) < eps, sign(eta_den) * eps, eta_den)
    eta_r1  <- eta_num / eta_den
    
    # 4) update Sigma (diagonal, per-dim) using previous mu (your original logic)
    # sigma_part1 <- colSums(b * centered^2)
    # sigma_part2 <- (colSums(centered)^2) / colSums(a)
    # sigma <- (1/n) * (part1 - part2)
    sum_b_center2 <- colSums(b_mat * centered_prev_mu^2, na.rm = TRUE)
    sum_center    <- colSums(centered_prev_mu, na.rm = TRUE)
    sum_a         <- colSums(a_mat, na.rm = TRUE)
    
    sum_a <- ifelse(abs(sum_a) < eps, sign(sum_a) * eps, sum_a)
    
    sigma_part1 <- sum_b_center2
    sigma_part2 <- (sum_center^2) / sum_a
    
    sigma_r1 <- (sigma_part1 - sigma_part2) / n
    sigma_r1 <- posify(sigma_r1)  # enforce non-negativity
    
    # commit updates for next iteration
    alpha_r <- alpha_r1
    mu_r    <- mu_r1
    eta_r   <- eta_r1
    Sigma_r <- sigma_r1
  }
  
  list(
    alpha_mle = as.numeric(alpha_r),
    mu_mle    = as.numeric(mu_r),
    eta_mle   = as.numeric(eta_r),
    Sigma_mle = as.numeric(Sigma_r),
    iter      = num_iterations
  )
}
