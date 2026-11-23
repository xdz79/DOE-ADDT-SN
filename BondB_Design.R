# In this procedure we aim to find the optimal design for predicting the p-th quantile
#
# Using ICA to find the optimal design for SNNAM
# We assume that v = (x, d), where x and d are both design variables
rm(list = ls())
library(sn)
library(ICAOD)
library(xtable)

source("merge_function.R")   # requires the merge_design_points function
# ----------------------------- Global constants -----------------------------
# Model parameter MLEs
theta_hat <- c(gamma0 = 3.9700, gamma1 = -0.2078, gamma2 = 0.6259,
               sigma = 0.1938, lambda = 1.0529)

Lowerx <- c(50, 3); Upperx <- c(70, 6) # the design region



# Transformation constants (from the literature)
x_bar <- -34.833      # center of the transformed temperature
tau_bar <- 2.455      # center of the transformed time

# Transformed value corresponding to the use condition (25 degrees C)
s0 <- 25
x0 <- -11605 / (s0 + 273.15)   # or use the formula directly; ensure consistency with the definition of x_bar

# Critical value y_f (e.g., log of 40 Newtons? Here it is assumed to be already log-transformed)
yf <- log(40)         # since the response is log(strength), the critical value is also log-transformed
p <- 0.01             # lifetime quantile of interest (0.01)
# ----------------------------- Mean function -----------------------------
# Input: parameter theta = (gamma0, gamma1, gamma2), temperature s, time t (weeks)
# Output: mean of the log strength, mu
mu_linear <- function(theta, s, t) {
  gamma0 <- theta[1]; gamma1 <- theta[2]; gamma2 <- theta[3]
  x <- -11605 / (s + 273.15)
  tau <- sqrt(t)
  gamma0 + gamma1 * ( exp((x - x_bar) * gamma2) * tau - tau_bar )
}

# Fisher information matrix of the skew-normal ADDT model for Bond B
# Design variables: temperature s and inspection interval d (both positive)
# The design vector x has length 2K: the first K entries are temperature values,
# the last K entries are the corresponding inspection intervals
# The weight vector w has length K and gives the proportion of runs for each (s, d) combination
FIM_SNNA_2D <- function(x, w, param) {
  # Extract parameters (theta = gamma0, gamma1, gamma2, sigma, lambda)
  gamma0 <- param[1]; gamma1 <- param[2]; gamma2 <- param[3]
  sigma <- param[4]; lambda <- param[5]
  
  npoint <- length(x) / 2
  s_vals <- x[1:npoint]               # temperature values
  d_vals <- x[(npoint + 1):(2 * npoint)]  # inspection intervals (weeks)
  
  # Constants
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
  # Initialize the accumulation matrices (3x3 corresponding to gamma0, gamma1, gamma2)
  SXXT <- matrix(0, ncol = 3, nrow = 3)
  SX   <- matrix(0, ncol = 1, nrow = 3)
  
  # Global variables: the number of observations l and the time points T_set (weeks)
  # must be defined outside this function. It is assumed that l and T_set exist,
  # and that T_set is a vector of length l (weeks). The constants x_bar, tau_bar, etc.
  # are also required.
  for (i in 1:npoint) {
    s_i <- s_vals[i]
    d_i <- d_vals[i]
    # Transformed temperature
    x_i <- -11605 / (s_i + 273.15)
    # Construct all observation time points (equally spaced) for this design point
    times <- seq(d_i, l * d_i, by = d_i)   # weeks
    tau_vals <- sqrt(times)                # transformed time
    
    FXi <- matrix(0, 3, 3)
    fxi <- matrix(0, 3, 1)
    
    for (j in 1:l) {
      tau <- tau_vals[j]
      # Mean function mu (on the original scale, i.e., log strength)
      mu_val <- gamma0 + gamma1 * ( exp((x_i - x_bar) * gamma2) * tau - tau_bar )
      # Gradient vector (d mu / d gamma0, d mu / d gamma1, d mu / d gamma2)
      dmu_dg0 <- 1
      dmu_dg1 <- exp((x_i - x_bar) * gamma2) * tau - tau_bar
      dmu_dg2 <- gamma1 * exp((x_i - x_bar) * gamma2) * tau * (x_i - x_bar)
      grad_mu <- c(dmu_dg0, dmu_dg1, dmu_dg2)
      FXi <- FXi + outer(grad_mu, grad_mu)
      fxi <- fxi + grad_mu
    }
    SXXT <- SXXT + w[i] * FXi
    SX   <- SX   + w[i] * fxi
  }
  
  # Construct the blocks of the information matrix
  I_theta2      <- ((1 + lambda^2 * a0) / (l * sigma^2)) * SXXT
  I_sig_theta   <- (1/(l * sigma^2) * (b * lambda^3 * (1 - a0*(1+lambda^2)) / (1+lambda^2)^1.5 + lambda^2 * a1)) * SX
  I_lamb_theta  <- -1/(l * sigma) * (b * a0 * lambda^2 / (1+lambda^2)^1.5 + lambda * a1) * SX
  
  I_sig2 <- (2 + lambda^2 * a2 + (a0*b^2*lambda^4 - b^2*lambda^2)/(1+lambda^2) -
               2*b^2*lambda^4/(1+lambda^2)^2 - 2*a1*b*lambda^3/(1+lambda^2)^0.5) / (sigma^2)
  I_sig_lamb <- -(lambda * a2 - a1*b*lambda^4/(1+lambda^2)^1.5 - a0*b^2*lambda^3/(1+lambda^2)^2 +
                    b^2*lambda*(1+2*lambda^2)/(1+lambda^2)^3) / sigma
  I_lamb2 <- a2 + 2*a1*b*lambda/(1+lambda^2)^1.5 + b^2*(a0*lambda^2 - 1)/(1+lambda^2)^3
  
  # Assemble the 5x5 Fisher information matrix
  FIM <- rbind(
    cbind(I_theta2, I_sig_theta, I_lamb_theta),
    cbind(t(I_sig_theta), I_sig2, I_sig_lamb),
    cbind(t(I_lamb_theta), I_sig_lamb, I_lamb2)
  )
  return(FIM)
}


num_deriv_quantile <- function(lambda, p, h = 1e-5) {
  (qsn(p, alpha = lambda + h, solver = "RFB") - 
     qsn(p, alpha = lambda - h, solver = "RFB")) / (2 * h)
}
# C-optimal design --------------------------------------------------------
#
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


# Part 2 D-opt -------------------------------------------------------------



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
# Table 6 (left part)

# Initial parameters (MLEs)
init_par <- theta_hat

# Number of observations and time points (to be used inside the FIM function)
# Note: FIM_SNNA_2D needs to know l and T_set, so they are set as global variables here
# (or passed as arguments)
l <- 3  # the number of observation times: 3, 4, 5


# T_set <- seq(6, 30, 6)   # 6, 12, 18, 24, 30 weeks (example)

# Optimization with ICAOD
Kv <- 6   # maximum number of design points (points will be merged afterwards)

# c-optimal design
res_c <- locally(
  fimfunc = FIM_SNNA_2D,
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
print(l)
merge_2d_design_points(res = res_c, Kv = Kv, threshold = 1, digits = 4)
print(res_c$design$min_cost)


# Case 1 # Table 6 (right part)
# D-optimal design
res_d <- locally(
  fimfunc = FIM_SNNA_2D,
  lx = Lowerx, ux = Upperx,
  inipars = init_par,
  iter = 1000, k = Kv,
  crtfunc = Dopt_criterion,
  sensfunc = Dopt_sensitivity,
  ICA.control = list(
    checkfreq = 200,
    stop_rule = "equivalence",
    plot_cost = TRUE,
    rseed = 1
  )
)

print(l)
merge_2d_design_points(res = res_d, Kv = Kv, threshold = 1, digits = 4)
print(res_d$design$min_cost)


# case 2 Table 8 (left part)
# Number of observations and time points (to be used inside the FIM function)
# Note: FIM_SNNA_2D needs to know l and T_set, so they are set as global variables here
# (or passed as arguments)
l <- 5  # the number of observation times: 3, 4, 5

# Initial parameters (MLEs)

lambda <- 100  # -100, -10, -1, -0.1, 0.1, 1, 10, 100

init_par <- c(theta_hat[1:4], lambda)



# T_set <- seq(6, 30, 6)   # 6, 12, 18, 24, 30 weeks (example)

# Optimization with ICAOD
Kv <- 6   # maximum number of design points (points will be merged afterwards)

# c-optimal design
# c-optimal design
res_c <- locally(
  fimfunc = FIM_SNNA_2D,
  lx = Lowerx, ux = Upperx,
  inipars = init_par,
  iter = 1000, k = Kv,
  crtfunc = Copt_criterion,
  sensfunc = Copt_sensitivity,
  ICA.control = list(
    #   ncount = 50,
    checkfreq = 200,
    stop_rule = "equivalence",
    plot_cost = TRUE
  )
)

# Output the c-optimal design (merging nearby points)
print(l)
merge_2d_design_points(res = res_c, Kv = Kv, threshold = 1, digits = 4)
print(res_c$design$min_cost)



# D-optimal design
res_d <- locally(
  fimfunc = FIM_SNNA_2D,
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

print(l)
merge_2d_design_points(res = res_d, Kv = Kv, threshold = 1, digits = 4)
print(res_d$design$min_cost)

## Case 3
# Table 8 (left part, efficiencies)

l <- 5  # the number of observation times: 3, 4, 5


lambda <- 100  # -100, -10, -1, -0.1, 0.1, 1, 10, 100
init_par <- c(theta_hat[1:4], lambda)

xsn <- c(57.8693,
         70,
         6,
         3)
# supports of the optimal designs for ADDTM-SN
wsn <- c(0.6698,
         0.3302) # weights of the optimal designs for ADDTM-SN

xn <- c(54.5196, 
        70,
        6,
        3.5358) # supports of the optimal designs for ADDTM-N
wn <- c(0.7523,
        0.2477) # weights of the optimal designs for ADDTM-N
EFF <- Copt_criterion(x = xsn, w = wsn, param = init_par, fimfunc = FIM_SNNA_2D)/ 
  Copt_criterion(x = xn, w = wn, param = init_par, fimfunc = FIM_SNNA_2D)
print(EFF)


# Table 8 (right part, efficiencies)
l <- 5  # the number of observation times: 3, 4, 5


lambda <- 100  # -100, -10, -1, -0.1, 0.1, 1, 10, 100

init_par <- c(theta_hat[1:4], lambda)

xsn <- c(50,
         56.7795,
         70,
         3,
         6,
         6) # supports of the optimal designs for ADDTM-SN
wsn <- c(0.1792,
         0.3822,
         0.4386) # weights of the optimal designs for ADDTM-SN

xn <- c(50,
        56.7808,
        70,
        3,
        6,
        6) # supports of the optimal designs for ADDTM-N
wn <- c(0.1788,
        0.3822,
        0.4390) # weights of the optimal designs for ADDTM-N
EFF <- exp((Dopt_criterion(x = xsn, w = wsn, param = init_par, fimfunc = FIM_SNNA_2D) -
              Dopt_criterion(x = xn, w = wn, param = init_par, fimfunc = FIM_SNNA_2D))/5)
print(EFF)