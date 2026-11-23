# Clear the environment and load packages
rm(list = ls())
library(sn)
library(ICAOD)
library(xtable)

source("merge_function.R")   # requires the merge_design_points function
# ----------------------------- Global constants -----------------------------
# Model parameter MLEs
theta_hat <- c(gamma0 = 3.9700, gamma1 = -0.2078, gamma2 = 0.6259,
               sigma = 0.1938, lambda = 1.0529)

# Design variable range (temperature, °C)
Lowerx <- 50
Upperx <- 70

# Transformation constants (from the literature)
x_bar <- -34.833      # center of the transformed temperature
tau_bar <- 2.455      # center of the transformed time

# Transformed value corresponding to the use condition (25°C)
s0 <- 25
x0 <- -11605 / (s0 + 273.15)   # or use the formula directly; note consistency with x_bar

# Critical value y_f (e.g., log of 40 Newtons? Here it is assumed to be already log-transformed)
yf <- log(40)         # since the response is log(strength), the critical value is also log-transformed
p <- 0.01             # lifetime quantile of interest (0.01)

# # Number of observations and time points in the experiment (example: l = 5, d = 6; adjust as needed)
# l <- 6                # number of observations at each stress level
# d <- 5                # time interval (weeks)
# T_set <- seq(d, l * d, d)   # observation time points (weeks)
# tau_set <- sqrt(T_set)       # transformed time

# ----------------------------- Mean function -----------------------------
# Input: parameter theta = (gamma0, gamma1, gamma2), temperature s, time t (weeks)
# Output: mean of the log strength, mu
mu_linear <- function(theta, s, t) {
  gamma0 <- theta[1]; gamma1 <- theta[2]; gamma2 <- theta[3]
  x <- -11605 / (s + 273.15)
  tau <- sqrt(t)
  gamma0 + gamma1 * ( exp((x - x_bar) * gamma2) * tau - tau_bar )
}

# ----------------------------- Fisher information matrix -----------------------------
# Note: this function needs to be derived for the new model. The code below is a placeholder;
# replace it with the correct analytical expression when using it.
# Input: vector of design points x (temperature), weight vector w, parameter vector
#        param = (gamma0, gamma1, gamma2, sigma, lambda)
# Output: 5x5 Fisher information matrix (based on skew-normal errors)
FIM_SNNA <- function(x, w, param) {
  # Extract parameters
  gamma0 <- param[1]; gamma1 <- param[2]; gamma2 <- param[3]
  sigma <- param[4]; lambda <- param[5]
  
  # Define constants
  b <- sqrt(2/pi)
  
  # Compute a0, a1, a2 (note the use of lambda)
  if (abs(lambda) < 200) {
    f.a <- function(y, lambda, k) y^k * dsn(y, 0, 1, lambda) * (zeta(1, lambda * y))^2
    err <- .Machine$double.eps^0.5
    a0 <- integrate(f.a, -Inf, Inf, lambda = lambda, k = 0, 
                    rel.tol = err)$value
    a1 <- integrate(f.a, -Inf, Inf, lambda = lambda, k = 1, 
                    rel.tol = err)$value
    a2 <- integrate(f.a, -Inf, Inf, lambda = lambda, k = 2, 
                    rel.tol = err)$value
  } else {
    u <- 1 + 8 * (lambda/pi)^2
    b_const <- sqrt(2/pi)
    a0 <- 1.019149098 * b_const^2 / sqrt(u)
    a1 <- 1.020466516 * (-lambda * b_const^3 / sqrt(u^3 * (1 + lambda^2/u)))
    a2 <- 1.009258704 * b_const^2 / sqrt(u)^3
  }
  
  # Initialize the accumulation matrices
  SXXT <- matrix(0, ncol = 3, nrow = 3)
  SX   <- matrix(0, ncol = 1, nrow = 3)
  
  for (i in 1:length(x)) {
    FXi <- matrix(0, ncol = 3, nrow = 3)
    fxi <- matrix(0, ncol = 1, nrow = 3)
    
    xsi <- -11605 / (x[i] + 273.15)   # transformed temperature
    
    for (j in 1:l) {
      tau <- sqrt(T_set[j])            # transformed time
      # Mean function (assuming mu_linear is already defined)
      mu_val <- mu_linear(c(gamma0, gamma1, gamma2), x[i], T_set[j])  # note: pass the original temperature and time
      # Gradient vector
      fxij <- c(1,
                (mu_val - gamma0) / gamma1,
                gamma1 * exp(gamma2 * (xsi - x_bar)) * tau * (xsi - x_bar))
      FXi <- FXi + outer(fxij, fxij)
      fxi <- fxi + fxij
    }
    SXXT <- SXXT + w[i] * FXi
    SX   <- SX   + w[i] * fxi
  }
  
  # Construct the blocks of the information matrix
  I_theta2      <- ((1 + lambda^2 * a0) / (l * sigma^2)) * SXXT
  I_sig_theta   <- (1/(l * sigma^2) * (b * lambda^3 * 
                                         (1 - a0*(1+lambda^2)) / (1+lambda^2)^1.5 + lambda^2 * a1)) * SX
  I_lamb_theta  <- -1/(l * sigma) * (b * a0 * lambda^2 / (1+lambda^2)^1.5 + lambda * a1) * SX
  
  I_sig2 <- (2 + lambda^2 * a2 + (a0*b^2*lambda^4 - b^2*lambda^2)/(1+lambda^2) -
               2*b^2*lambda^4/(1+lambda^2)^2 - 2*a1*b*lambda^3/(1+lambda^2)^0.5) / (sigma^2)
  I_sig_lamb <- -(lambda * a2 - a1*b*lambda^4/(1+lambda^2)^1.5 - a0*b^2*lambda^3/(1+lambda^2)^2 +
                    b^2*lambda*(1+2*lambda^2)/(1+lambda^2)^3) / sigma
  I_lamb2 <- a2 + 2*a1*b*lambda/(1+lambda^2)^1.5 + b^2*(a0*lambda^2 - 1)/(1+lambda^2)^3
  
  # Assemble the 5x5 information matrix
  FIM <- rbind(
    cbind(I_theta2, I_sig_theta, I_lamb_theta),
    cbind(t(I_sig_theta), I_sig2, I_sig_lamb),
    cbind(t(I_lamb_theta), I_sig_lamb, I_lamb2)
  )
  return(FIM)
}

# ----------------------------- c-optimal criterion -----------------------------
# Expression for the quantile t_p (in weeks)
quantile_tp <- function(param, yf, p, x0, x_bar, tau_bar) {
  gamma0 <- param[1]; gamma1 <- param[2]; gamma2 <- param[3]
  sigma <- param[4]; lambda <- param[5]
  
  b <- sqrt(2/pi)
  delta <- lambda / sqrt(1 + lambda^2)
  mu_z <- b * delta
  z_lambda_p <- qsn(p, alpha = lambda, solver = "RFB")
  
  A <- yf - gamma0 + gamma1 * tau_bar + sigma * (mu_z - z_lambda_p)
  denom <- gamma1 * exp((x0 - x_bar) * gamma2)
  t_p <- (A / denom)^2
  return(t_p)
}

num_deriv_quantile <- function(lambda, p, h = 1e-5) {
  (qsn(p, alpha = lambda + h, solver = "RFB") - qsn(p, alpha = lambda - h, solver = "RFB")) / (2 * h)
}
# Gradient vector (5 components) - based on the formulas provided
grad_tp <- function(param, yf, p, x0, x_bar, tau_bar) {
  gamma0 <- param[1]; gamma1 <- param[2]; gamma2 <- param[3]
  sigma <- param[4]; lambda <- param[5]
  
  b <- sqrt(2/pi)
  delta <- lambda / sqrt(1 + lambda^2)
  mu_z <- b * delta
  z_lambda_p <- qsn(p, alpha = lambda, solver = "RFB")
  dmuz_dlambda <- b * (1 + lambda^2)^(-1.5)
  
  A <- yf - gamma0 + gamma1 * tau_bar + sigma * (mu_z - z_lambda_p)
  denom <- gamma1 * exp((x0 - x_bar) * gamma2)
  t_p <- (A / denom)^2
  h_t <- sqrt(t_p)   # h_t(t_p) = sqrt(t_p)
  
  d_gamma0 <- -2 * h_t / denom
  d_gamma1 <- -2 * h_t * (h_t/gamma1 - tau_bar / denom)
  d_gamma2 <- -2 * h_t^2 * (x0 - x_bar)
  d_sigma  <-  2 * h_t * (mu_z - z_lambda_p) / denom
  d_lambda <- 2 * h_t * sigma*(dmuz_dlambda-num_deriv_quantile(lambda, p, h = 1e-5))/ denom
  
  return(c(d_gamma0, d_gamma1, d_gamma2, d_sigma, d_lambda))
}

# c-optimal criterion value (asymptotic variance)
Copt_criterion <- function(x, w, param, fimfunc) {
  fim <- fimfunc(x, w, param)
  # Regularized inverse
  eps <- 1e-8
  fim_reg <- fim + diag(eps, nrow(fim))
  fim_inv <- solve(fim_reg)
  c_vec <- grad_tp(param, yf, p, x0, x_bar, tau_bar)
  avar <- t(c_vec) %*% fim_inv %*% c_vec
  return(avar)
}

# c-optimal sensitivity function (equivalence theorem)
Copt_sensitivity <- function(xi_x, x, w, param, fimfunc) {
  fim <- fimfunc(x, w, param)
  eps <- 1e-8
  fim_reg <- fim + diag(eps, nrow(fim))
  M_inv <- solve(fim_reg)
  M_x <- fimfunc(xi_x, w = 1, param)
  c_vec <- grad_tp(param, yf, p, x0, x_bar, tau_bar)
  Cv_norm <- sqrt(sum(c_vec^2))
  c_vec <- c_vec / Cv_norm 
  sens <- t(c_vec) %*% M_inv %*% M_x %*% M_inv %*% c_vec - t(c_vec) %*% M_inv %*% c_vec
  return(sens)
}

# ----------------------------- D-optimal criterion -----------------------------
Dopt_criterion <- function(x, w, param, fimfunc) {
  fim <- fimfunc(x, w, param)
  return(-log(det(fim + diag(1e-8, nrow(fim)))))
}

Dopt_sensitivity <- function(xi_x, x, w, param, fimfunc) {
  fim <- fimfunc(x, w, param)
  fim_reg <- fim + diag(1e-8, nrow(fim))
  M_inv <- solve(fim_reg)
  M_x <- fimfunc(xi_x, w = 1, param)
  sens <- sum(diag(M_inv %*% M_x)) - length(param)
  return(sens)
}

# ----------------------------- Optimization settings -----------------------------

# case 1

# Table 5
# Initial parameters (MLEs)
init_par <- theta_hat

# Number of observations and time points (to be used inside the FIM function)
# Note: FIM_SNNA needs to know l and T_set, so they are set as global variables here
# (or passed as arguments)
l <-5   # the number of observation times: 3, 4, 5
d <- 6  # the distance between inspections: 3, 4, 5, 6
# equal inspection interval
T_set <- seq(d, l * d, d)

# T_set <- seq(6, 30, 6)   # 6, 12, 18, 24, 30 weeks (example)

# Optimization with ICAOD
Kv <- 6   # maximum number of design points (points will be merged afterwards)

# c-optimal design
res_c <- locally(
  fimfunc = FIM_SNNA,
  lx = Lowerx, ux = Upperx,
  inipars = init_par,
  iter = 1000, k = Kv,
  crtfunc = Copt_criterion,
  sensfunc = Copt_sensitivity,
  ICA.control = list(
    #   ncount = 50,
    checkfreq = 300,
    stop_rule = "equivalence",
    plot_cost = TRUE
  )
)

# Output the c-optimal design (merging nearby points)

merge_design_points(res_c, Kv, threshold = 0.1)
res_c$design$min_cost


# D-optimal design
res_d <- locally(
  fimfunc = FIM_SNNA,
  lx = Lowerx, ux = Upperx,
  inipars = init_par,
  iter = 1000, k = Kv,
  crtfunc = Dopt_criterion,
  sensfunc = Dopt_sensitivity,
  ICA.control = list(
    checkfreq = 200,
    stop_rule = "equivalence",
    plot_cost = FALSE,
    rseed = 1
  )
)

merge_design_points(res_d, Kv, threshold = 0.1)
res_d$design$min_cost


#Table 7 (design part)
# Number of observations and time points (to be used inside the FIM function)
# Note: FIM_SNNA needs to know l and T_set, so they are set as global variables here
# (or passed as arguments)
l <- 5  # the number of observation times: 3, 4, 5
d <- 6  # the distance between inspections: 3, 4, 5, 6
## case 2
# Initial parameters (MLEs)

lambda <- 10    # -100, -10, -1, -0.1, 0.1, 1, 10, 100

init_par <- c(theta_hat[1:4], lambda)


# equal inspection interval
T_set <- seq(d, l * d, d)

# T_set <- seq(6, 30, 6)   # 6, 12, 18, 24, 30 weeks (example)

# Optimization with ICAOD
Kv <- 6   # maximum number of design points (points will be merged afterwards)

# c-optimal design
res_c <- locally(
  fimfunc = FIM_SNNA,
  lx = Lowerx, ux = Upperx,
  inipars = init_par,
  iter = 1000, k = Kv,
  crtfunc = Copt_criterion,
  sensfunc = Copt_sensitivity,
  ICA.control = list(
    ncount = 100,
    checkfreq = 500,
    stop_rule = "equivalence",
    plot_cost = FALSE
  )
)

# Output the c-optimal design (merging nearby points)
merge_design_points(res_c, Kv, threshold = 0.1)
res_c$design$min_cost


# D-optimal design
res_d <- locally(
  fimfunc = FIM_SNNA,
  lx = Lowerx, ux = Upperx,
  inipars = init_par,
  iter = 1000, k = Kv,
  crtfunc = Dopt_criterion,
  sensfunc = Dopt_sensitivity,
  ICA.control = list(
    checkfreq = 200,
    stop_rule = "equivalence",
    plot_cost = FALSE,
    rseed = 1
  )
)

merge_design_points(res_d, Kv, threshold = 0.1)
res_d$design$min_cost



## Case 3 
## Table 7(left part, efficiencies)

l <- 5  # the number of observation times: 3, 4, 5
d <- 6  # the distance between inspections: 3, 4, 5, 6
## case 2
# Initial parameters (MLEs)

lambda <- -10  # -100, -10, -1, -0.1, 0.1, 1, 10, 100

init_par <- c(theta_hat[1:4], lambda)


# equal inspection interval
T_set <- seq(d, l * d, d)

xsn <- c(57.4524,
         70
) # supports of the optimal designs for ADDTM-SN
wsn <- c(0.7489,
         0.2511) # weights of the optimal designs for ADDTM-SN

xn <- c(55.8078,
        70
)      # supports of the optimal designs for ADDTM-N
wn <- c(0.7773,
        0.2227) # weights of the optimal designs for ADDTM-N
EFF <- Copt_criterion(x = xsn, w = wsn, param = init_par, fimfunc = FIM_SNNA)/
  Copt_criterion(x = xn, w = wn, param = init_par, fimfunc = FIM_SNNA)
print(EFF)