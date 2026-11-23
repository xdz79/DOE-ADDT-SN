# In this procedure, we aim to find the optimal x_m under the robust design criteria
# discussed in Section 5. The polymer material example is considered.
# Table 9

rm(list = ls())

# Parameter setting

# Define the 5x5 matrix A,
test_t <- c(8, 25, 75, 130, 180)
theta_sn <- c(-4.965450e+00,  1.370726e+03,  1.454751e+01, -5.956617e+03,
              2.584227e-02,  1.249420e+00)

N <- 75

A_matrix <- function(x, theta_nstar, theta_sn, w) {
  # A
  m <- 3
  x_l <- 25
  x_h <- 80
  x_m <- x
  test_s <- c(x_l, x_m, x_h)
  A_m <- matrix(0, nrow = 5, ncol = 5)
  
  for (i in 1:3) {
    for (j in 1:5) {
      xi <- 1/(273.15 + test_s[i])
      
      alpha_sni <- exp(theta_sn[1] + theta_sn[2] * xi)
      beta_sni <- exp(theta_sn[3] + theta_sn[4] * xi)
      eta_snij <- exp(-beta_sni * test_t[j]^alpha_sni)
      
      alpha_nstari <- exp(theta_nstar[1] + theta_nstar[2] * xi)
      beta_nstari <- exp(theta_nstar[3] + theta_nstar[4] * xi)
      eta_nstarij <- exp(-beta_nstari * test_t[j]^alpha_nstari)
      
      etan11 <- -alpha_nstari * beta_nstari * eta_nstarij * test_t[j]^alpha_nstari *
        log(test_t[j]) * (1 + alpha_nstari * log(test_t[j]) * (1 - beta_nstari * test_t[j]^alpha_nstari))
      etan12 <- etan11 * xi
      etan13 <- -alpha_nstari * beta_nstari * eta_nstarij * test_t[j]^alpha_nstari *
        log(test_t[j]) * (1 - beta_nstari * test_t[j]^alpha_nstari)
      etan14 <- etan13 * xi
      etan22 <- etan11 * xi^2
      etan23 <- etan14
      etan24 <- etan23 * xi
      etan33 <- -beta_nstari * eta_nstarij * test_t[j]^alpha_nstari *
        (1 - beta_nstari * test_t[j]^alpha_nstari)
      etan34 <- etan33 * xi
      etan44 <- etan33 * xi^2
      
      eta_matrix <- matrix(c(etan11, etan12, etan13, etan14,
                             etan12, etan22, etan23, etan24,
                             etan13, etan23, etan33, etan34,
                             etan14, etan24, etan34, etan44), nrow = 4, ncol = 4)
      etan1 <- -alpha_nstari * beta_nstari * eta_nstarij * test_t[j]^alpha_nstari * log(test_t[j])
      etan2 <- etan1 * xi
      etan3 <- -beta_nstari * eta_nstarij * test_t[j]^alpha_nstari
      etan4 <- etan3 * xi
      etan_v <- matrix(c(etan1, etan2, etan3, etan4), ncol = 1)
      
      Ni <- N/5 * w[i]
      #### n_i = N/5 * w[i]
      A_m[1:4, 1:4] <- A_m[1:4, 1:4] + Ni/theta_nstar[5]^2 * ((eta_snij - eta_nstarij) * eta_matrix -
                                                                etan_v %*% t(etan_v))
      A_m[5, 5] <- A_m[5, 5] + Ni * (1/theta_nstar[5]^2 - 3/theta_nstar[5]^4
                                     * ((1 - 2/pi * theta_sn[6]^2/(1 + theta_sn[6]^2)) * theta_sn[5]^2 +
                                          (eta_snij - eta_nstarij)^2))
      
      A_m[1:4, 5] <- A_m[1:4, 5] - 2 * Ni/theta_nstar[5]^3 * (eta_snij - eta_nstarij) * etan_v
      A_m[5, 1:4] <- t(A_m[1:4, 5])
    }
  }
  return(A_m)
}


B_matrix <- function(x, theta_nstar, theta_sn, w) {
  # B
  m <- 3
  x_l <- 25
  x_h <- 80
  x_m <- x
  test_s <- c(x_l, x_m, x_h)
  B_m <- matrix(0, nrow = 5, ncol = 5)
  B_m11 <- matrix(0, nrow = 4, ncol = 4)
  B_m12 <- matrix(0, nrow = 4, ncol = 1)
  
  B_m21 <- matrix(0, nrow = 4, ncol = 1)
  B_m221 <- matrix(0, nrow = 4, ncol = 1)
  B_m222 <- 0
  
  B_m31 <- 0
  B_m32 <- 0
  
  for (i in 1:3) {
    for (j in 1:5) {
      xi <- 1/(273.15 + test_s[i])
      Ni <- N/5 * w[i]
      
      alpha_sni <- exp(theta_sn[1] + theta_sn[2] * xi)
      beta_sni <- exp(theta_sn[3] + theta_sn[4] * xi)
      eta_snij <- exp(-beta_sni * test_t[j]^alpha_sni)
      
      alpha_nstari <- exp(theta_nstar[1] + theta_nstar[2] * xi)
      beta_nstari <- exp(theta_nstar[3] + theta_nstar[4] * xi)
      eta_nstarij <- exp(-beta_nstari * test_t[j]^alpha_nstari)
      
      
      etan1 <- -alpha_nstari * beta_nstari * eta_nstarij * test_t[j]^alpha_nstari * log(test_t[j])
      etan2 <- etan1 * xi
      etan3 <- -beta_nstari * eta_nstarij * test_t[j]^alpha_nstari
      etan4 <- etan3 * xi
      etan_v <- matrix(c(etan1, etan2, etan3, etan4), ncol = 1)
      
      delta <- theta_sn[6] / sqrt(1 + theta_sn[6]^2)
      b <- sqrt(2/pi)
      
      # Common constants
      const1 <- (1 - b^2 * delta^2) * theta_sn[5]^2
      const2 <- 1/2 * (4 - pi) * theta_sn[5]^3 * b^3 * delta^3
      
      B_m11 <- B_m11 + Ni/(theta_nstar[5]^4) * const1 * etan_v %*% t(etan_v)
      
      B_m12 <- B_m12 + Ni/(theta_nstar[5]^2) * (eta_snij - eta_nstarij) * etan_v
      
      B_m21 <- B_m21 + Ni/(theta_nstar[5]^5) *
        (const2 + 2 * (eta_snij - eta_nstarij) * const1) * etan_v
      
      B_m221 <- B_m221 + Ni/theta_nstar[5]^2 * (eta_snij - eta_nstarij) * etan_v
      
      B_m222 <- B_m222 + (-Ni/theta_nstar[5] + Ni/theta_nstar[5]^3 *
                            (const1 + (eta_snij - eta_nstarij)^2))
      
      B_m31 <- B_m31 + Ni/theta_nstar[5]^6 * (4 * const1 * (eta_snij - eta_nstarij)^2 +
                                                4 * const2 * (eta_snij - eta_nstarij) + 2 * theta_sn[5]^4 *
                                                (1 - 2 * b^2 * delta^2 + (pi - 2) * b^4 * delta^4))
      B_m32 <- B_m32 + (-Ni/theta_nstar[5] + Ni/theta_nstar[5]^3 * (const1 + (eta_snij - eta_nstarij)^2))
      
    } # end j
  } # end i
  B_m[1:4, 1:4] <- B_m11 + B_m12 %*% t(B_m12)
  B_m[1:4, 5] <- B_m21 + B_m221 * B_m222
  B_m[5, 5] <- B_m31 + B_m32^2
  B_m[5, 1:4] <- t(B_m[1:4, 5])
  return(B_m)
}

C_matrix <- function(x, theta_nstar, theta_sn, w) {
  A <- A_matrix(x, theta_nstar, theta_sn, w)
  B <- B_matrix(x, theta_nstar, theta_sn, w)
  # Regularization
  # A <- A + diag(1e-8, nrow(A))
  A_inv <- solve(A)
  C <- A_inv %*% B %*% A_inv
  return(C)
}

# To estimate the quasi MLE under the design xi
b <- sqrt(2/pi)
Loglifun <- function(y, x, theta_sn, w) {
  # y[1] = gamma1, y[2] = gamma2, y[3] = gamma3, y[4] = gamma4, y[5] = sigma
  # x is the lower level
  TEM <- rep(c(25, x, 80), each = 5)
  TIME <- rep(c(8, 25, 75, 130, 180), 3)
  xsi <- 1/(273.15 + TEM)
  alphai <- exp(y[1] + y[2] * xsi)
  betai <- exp(y[3] + y[4] * xsi)
  etai <- exp(-betai * (TIME)^alphai)
  
  alphaisn <- exp(theta_sn[1] + theta_sn[2] * xsi)
  betaisn <- exp(theta_sn[3] + theta_sn[4] * xsi)
  Ni <- rep(N/5 * w, each = 5)
  
  RATIO1 <- exp(-betaisn * (TIME)^alphaisn)
  # 
  # ((1 - b^2 * delta^2) * sig0^2
  
  LLfun <- sum(Ni * (RATIO1 - etai)^2)
  
  return(LLfun)
}


# test 
x_ini <- c(-4.965450e+00,  1.370726e+03,  1.454751e+01, -5.956617e+03)

Estimate_thetanstar <- function(x, theta_sn, x_ini, w) {
  # 1. Estimate gamma (minimize the unweighted sum of squares using nlm)
  fit <- nlm(f = Loglifun, p = x_ini, x = x, theta_sn = theta_sn, w)
  gamma_est <- fit$estimate
  
  # 2. Compute the weighted sum of squares and the constant term to obtain the estimate of sigma
  TEM <- rep(c(25, x, 80), each = 5)
  TIME <- rep(c(8, 25, 75, 130, 180), 3)
  xsi <- 1/(273.15 + TEM)
  
  # eta under the working model
  alphai <- exp(gamma_est[1] + gamma_est[2] * xsi)
  betai  <- exp(gamma_est[3] + gamma_est[4] * xsi)
  eta_N <- exp(-betai * TIME^alphai)
  
  # eta under the true model
  alphaisn <- exp(theta_sn[1] + theta_sn[2] * xsi)
  betaisn  <- exp(theta_sn[3] + theta_sn[4] * xsi)
  eta_SN <- exp(-betaisn * TIME^alphaisn)
  
  # Constant term
  delta <- theta_sn[6] / sqrt(1 + theta_sn[6]^2)
  b <- sqrt(2/pi)
  const_term <- (1 - b^2 * delta^2) * theta_sn[5]^2
  
  Ni <- rep(N/5 * w, each = 5)
  
  # Total sum of squares (5 replicates for each (i, j))
  
  weighted_sq <- sum(Ni * (eta_SN - eta_N)^2)
  total_const <- N * const_term   # 15 (i, j) combinations
  S <- weighted_sq + total_const
  
  sigma_est <- sqrt(S / N)
  
  return(c(gamma_est, sigma_est))
}


# DM-optimal design criterion
# DM_cri <- function(x){
#   theta_nstar <- Estimate_thetanstar(x)
#   res <- -log(det(C_matrix(x, theta_nstar)))
#   return(res)
# }

DM_cri <- function(x, theta_sn, x_ini, w) {
  theta_nstar <- Estimate_thetanstar(x = x, theta_sn = theta_sn, x_ini = x_ini, w)
  C <- C_matrix(x, theta_nstar = theta_nstar, theta_sn = theta_sn, w = w)
  
  det_val <- tryCatch(det(C), error = function(e) 0)
  if (det_val <= 0) return(1e10)
  log(det_val)
}


# p-quantile under the normal model
t_p_normal <- function(theta_n, y_f, p, x0) {
  # theta_n: c(gamma1, gamma2, gamma3, gamma4, sigma)
  gamma1 <- theta_n[1]; gamma2 <- theta_n[2]
  gamma3 <- theta_n[3]; gamma4 <- theta_n[4]; sigma <- theta_n[5]
  alpha0 <- exp(gamma1 + gamma2 * x0)
  beta0  <- exp(gamma3 + gamma4 * x0)
  zp <- qnorm(p)
  A <- y_f - zp * sigma
  if (A <= 0) return(NA)
  (log(A) / (-beta0))^(1/alpha0)
}

grad_tp_normal <- function(theta, y_f, p, x0) {
  # theta: parameter vector c(gamma1, gamma2, gamma3, gamma4, sigma)
  gamma1 <- theta[1]
  gamma2 <- theta[2]
  gamma3 <- theta[3]
  gamma4 <- theta[4]
  sigma  <- theta[5]
  
  alpha0 <- exp(gamma1 + gamma2 * x0)   # alpha under the use condition
  beta0  <- exp(gamma3 + gamma4 * x0)   # beta under the use condition
  
  zp <- qnorm(p)                        # standard normal p-quantile
  A <- y_f - zp * sigma                 # intermediate quantity in the quantile formula
  
  # Numerical safeguard: A must be > 0, otherwise the quantile is undefined
  if (A <= 0) return(rep(NA, 5))
  
  t_p <- (log(A) / (-beta0))^(1 / alpha0)
  log_t <- log(t_p)
  
  # Gradient components
  d_gamma1 <- -t_p * log_t
  d_gamma2 <- -t_p * log_t * x0
  d_gamma3 <- -t_p / alpha0
  d_gamma4 <- -t_p / alpha0 * x0
  d_sigma  <- t_p * zp / (alpha0 * A * (-log(A)))
  
  c(d_gamma1, d_gamma2, d_gamma3, d_gamma4, d_sigma)
}

library(sn)
t_p_SN <- function(theta_sn, y_f, p, x0) {
  # theta_sn: c(gamma1, gamma2, gamma3, gamma4, sigma, lambda)
  gamma1 <- theta_sn[1]; gamma2 <- theta_sn[2]
  gamma3 <- theta_sn[3]; gamma4 <- theta_sn[4]
  sigma  <- theta_sn[5]; lambda <- theta_sn[6]
  alpha0 <- exp(gamma1 + gamma2 * x0)
  beta0  <- exp(gamma3 + gamma4 * x0)
  b <- sqrt(2/pi)
  delta <- lambda / sqrt(1 + lambda^2)
  mu_z <- b * delta
  z_lambda_p <- qsn(p, alpha = lambda)
  A <- y_f + sigma * mu_z - sigma * z_lambda_p
  if (A <= 0) return(NA)
  (log(A) / (-beta0))^(1/alpha0)
}


AMSE_criterion <- function(x,  y_f = 0.6, p = 0.1, 
                           x0 = 1/(273.15 + 25), 
                           theta_sn, x_ini, w = w) {
  # x  
  # Other arguments: number of runs, time points, critical value, quantile,
  # transformed value of the use condition, true parameters, working parameters
  
  theta_nstar <- Estimate_thetanstar(x = x, theta_sn = theta_sn, x_ini = x_ini, w)
  
  # 1. Compute the matrix C
  C <- C_matrix(x = x, theta_nstar = theta_nstar, theta_sn = theta_sn, w = w)    # should use theta_sn, theta_nstar, l, T_set internally
  
  
  # 2. Compute the gradient c (5-dimensional)
  c_vec <- grad_tp_normal(theta_nstar, y_f, p, x0)
  if (any(is.na(c_vec))) return(Inf)
  
  # 3. Asymptotic variance
  AVar <- as.numeric(t(c_vec) %*% C %*% c_vec)
  
  # 4. Asymptotic bias
  t_true <- t_p_SN(theta_sn, y_f, p, x0)
  t_work <- t_p_normal(theta_nstar, y_f, p, x0)
  if (is.na(t_true) || is.na(t_work)) return(Inf)
  ABias <- t_work - t_true
  
  # 5. AMSE
  AMSE <- AVar + ABias^2 
  return(AMSE)
}


# theta_sn <- c(-4.965450e+00,  1.370726e+03,  1.454751e+01, -5.956617e+03,
#                       2.584227e-02,  1.249420e+00)


# w <- c(1/3, 1/3, 1/3)
w <- c(4/7, 2/7, 1/7)
# w <- c(1/7, 2/7, 4/7)
theta_sn <- c(-4.965450e+00,  1.370726e+03,  1.454751e+01, -5.956617e+03,
              2.584227e-02,  1.249420e+00)
# theta_nstar <- c(-4.826, 1269.850, 13.171, -5403.990, 0.002)  

gam1_set <- -4.965450e+00
gam2_set <- 1.370726e+03
gam3_set <- c(14, 15)
# gam3_set <- 1.454751e+01
gam4_set <- c(-5900, -6000)
# gam4_set <- -5.956617e+03
gam5_set <- 2.584227e-02
lam_set <- c(-150, -15, -1.5, -0.15, 0.15, 1.5, 15, 150)
theta_sn_set <- expand.grid(gam1_set, gam2_set, gam3_set, gam4_set, gam5_set, lam_set)
# 

# theta_sn_set <- matrix(theta_sn, nrow = 1)
### grid method
# Robust_set <- matrix(0, nrow = nrow(theta_sn_set), ncol = ncol(theta_sn_set) + 2)
# for (i in 1:nrow(theta_sn_set)) {
# 
#   theta_sn <- unlist(theta_sn_set[i, ])
#   x_set <- seq(25, 79, 0.1)
#   yM_set <- sapply(x_set, AMSE_criterion, theta_sn = theta_sn, x_ini = x_ini, w = w)
#   plot(x_set, log(yM_set), type = "l")
#   print(x_set[which.min(yM_set)])
#   #
#   # x_set <- seq(25, 79, 1)
#   yD_set <- sapply(x_set, DM_cri, theta_sn = theta_sn, x_ini = x_ini, w = w)
#   plot(x_set, yD_set, type = "l")
#   print(x_set[which.min(yD_set)])
# 
#   Robust_set[i, ] <- c(theta_sn, x_set[which.min(yD_set)], x_set[which.min(yM_set)])
# }
# # write.csv(Robust_set, file = "Robust_set_grid.csv")
# Robust_set <- as.data.frame(Robust_set)
# 
# Robust_set <- Robust_set[order(Robust_set$V3, Robust_set$V4), ]
# Robust_set_s <- Robust_set[, c("V3", "V4", "V6", "V7", "V8")]

# write.csv(result_robust_s, file = "Robust_set_opt_s1.csv")

# library(xtable)
# xtable(Robust_set_s, digits = c(0, 0, 0, 2, 4, 4))


####### Generic one-dimensional optimization function
optimize_design <- function(criterion_fun, interval = c(25, 80), 
                            n_grid = 20, ...) {
  # criterion_fun: objective function (accepts x and the ... arguments)
  # interval: search interval
  # n_grid: number of grid points
  # ...: other arguments passed to criterion_fun (e.g., theta_sn, x_ini, y_f, p, x0)
  
  # 1. Coarse grid search
  grid_x <- seq(interval[1], interval[2], length.out = n_grid)
  grid_obj <- sapply(grid_x, function(x) {
    tryCatch(criterion_fun(x, ...), error = function(e) 1e10)
  })
  plot(grid_x, grid_obj)
  best_grid_idx <- which.min(grid_obj)
  best_grid_x <- grid_x[best_grid_idx]
  best_grid_val <- grid_obj[best_grid_idx]
  
  # 2. Local optimization (using optimize)
  # Set a narrower starting interval to ensure stable optimization
  # Center on the best grid point, take half the interval width on each side,
  # but do not exceed the original interval
  half_width <- diff(interval) / 4
  lower <- max(interval[1], best_grid_x - half_width)
  upper <- min(interval[2], best_grid_x + half_width)
  # If the interval degenerates to a point, expand it slightly
  if (lower == upper) {
    lower <- max(interval[1], lower - 0.1)
    upper <- min(interval[2], upper + 0.1)
  }
  
  opt <- tryCatch(
    optimize(f = criterion_fun, interval = c(lower, upper), 
             maximum = FALSE, ...),
    error = function(e) list(minimum = best_grid_x, objective = best_grid_val)
  )
  
  # 3. Output the result
  result <- list(
    minimum = opt$minimum,
    objective = opt$objective,
    grid_search = list(
      x = grid_x,
      obj = grid_obj,
      best_x = best_grid_x,
      best_obj = best_grid_val
    )
  )
  return(result)
}

# Example: optimize AMSE_criterion (requires theta_sn, x_ini, y_f, p, x0
# to be defined beforehand)
# Assume these variables are already defined in the global environment:
# theta_sn <- c(...)
# x_ini <- c(...)
# y_f <- 0.6; p <- 0.1; x0 <- 1/(273.15 + 25)

# Run the optimization

# w <- c(1/3, 1/3, 1/3) 

result_robust <- matrix(0, nrow = nrow(theta_sn_set), ncol = ncol(theta_sn_set) + 2)
for (i in 1:nrow(theta_sn_set)) {
  
  theta_sn <- unlist(theta_sn_set[i, ])
  
  amse_opt <- optimize_design(AMSE_criterion, interval = c(25, 80),
                              n_grid = 20,
                              y_f = 0.6, p = 0.1,
                              x0 = 1/(273.15 + 25),
                              theta_sn = theta_sn,
                              x_ini = x_ini, w = w)
  
  dm_opt <- optimize_design(DM_cri, interval = c(25, 80),
                            n_grid = 20,
                            theta_sn = theta_sn,
                            x_ini = x_ini, w = w)
  result_robust[i, ] <- c(theta_sn, dm_opt$minimum, amse_opt$minimum)
}

print(result_robust)
result_robust <- as.data.frame(result_robust)

result_robust <- result_robust[order(result_robust$V3, result_robust$V4), ]
result_robust_s <- result_robust[, c("V3", "V4", "V6", "V7", "V8")]

# write.csv(result_robust_s, file = "Robust_set_opt_s.csv")

library(xtable)
xtable(result_robust_s, digits = c(0, 0, 0, 2, 4, 4))
# View the results
# amse_opt$minimum
# amse_opt$objective