## Example
data_ <- rMGAL(n = 500,
                        mu_vec = c(3,2),
                        eta_vec =  c(1.4, 1),
                        gamma_parameters = c(1.5,1),
                        sigma_matrix = matrix(data = c(10, 5,5, 10), 
                                              nrow = 2, 
                                              ncol = 2))


## STAGE 1: marginal EM
param_estimate1 <- parameter_estimate_non_sym_k(data_Y = data_,
                                              alpha_r = c(1,1),
                                              eta_r = c(0, 0),
                                              mu_r = c(1,1),
                                              Sigma_r = c(1,1),
                                              num_iterations = 100) 
## STAGE 2: Covariance matrix
cov_estimate_init <- covariance_estimation_k(data_Y = data_,
                            alpha_r = round(param_estimate1$alpha_mle, 4),
                            eta_r = round(param_estimate1$eta_mle,4),
                            mu_r = round(param_estimate1$mu_mle, 4),
                            Sigma_r = diag(param_estimate1$Sigma_mle),
                            num_iterations = 50,
                            burn_in = 10,
                            M = 40)

covariance_matrix_est <- matrix(c(0,cov_estimate_init$Sigma_mle[1,2],
                                cov_estimate_init$Sigma_mle[1,2], 0), 
                                nrow = 2, ncol = 2) + diag(param_estimate1$Sigma_mle)
