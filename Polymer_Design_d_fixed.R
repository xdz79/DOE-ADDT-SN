# In this procedure we aim to find the optimal design for predicting the p-th quantile
# 
# Using ICA to find the optimal design for SNNAM
# We assume that v = (x, d); x is the design variable and d is fixed.

rm(list = ls())
source("merge_function.R")
library(sn)
library(ICAOD)
library(xtable)

Lowerx <- c(25); Upperx <- c(80) # the design region



# Information matrix for the quadratic linear regression model


# nominal values
#
param0 <- c(-4.965450e+00,  1.370726e+03,  1.454751e+01, -5.956617e+03,
            2.584227e-02,  1.249420e+00)
b <- sqrt(2 / pi)
xbar <- mean(1/(273.5 + c(50, 65, 80)))
# xbar <- mean(1/(273.5 + 50))
phi <- c(exp(param0[1] + param0[2] * xbar), exp(param0[3] + param0[4] * xbar), 
         param0[2], param0[4], param0[5],  param0[6])

# w <- c(1/4, 1/2, 1/4)



FIM_SNNAM <- function(x, w, param) {
  # x : the design vector
  # x1 = [x1, x2, ..., xm], x2 = [d1, d2, ..., dm]
  # w : the weight vector
  # xori <- [50, 65, 80, 36, 36, 36]
  # 
  
  gam1 <- param[1]
  gam2 <- param[2]
  gam3 <- param[3]
  gam4 <- param[4]
  sigm <- param[5]
  lamb <- param[6]
  
  if (abs(lamb) < 200) {
    f.a <- function(y, lamb, k) y^k * dsn(y, 0, 1, lamb) * (zeta(1, lamb * y))^2
    err <- .Machine$double.eps^0.5
    a0 <- integrate(f.a, -Inf, Inf, lamb = lamb, k = 0, 
                    rel.tol = err)$value
    a1 <- integrate(f.a, -Inf, Inf, lamb = lamb, k = 1, 
                    rel.tol = err)$value
    a2 <- integrate(f.a, -Inf, Inf, lamb = lamb, k = 2, 
                    rel.tol = err)$value
  } else {
    u <- 1 + 8 * (lamb/pi)^2
    b <- sqrt(2/pi)
    a0 <- 1.019149098 * b^2/sqrt(u)
    a1 <- 1.020466516 * (-lamb * b^3/sqrt(u^3 * (1 + lamb^2/u)))
    a2 <- 1.009258704 * b^2/sqrt(u)^3
  }
  #
  
  
  
  SXXT <- matrix(0, ncol = 4, nrow = 4) 
  SX <- matrix(0, ncol = 1, nrow = 4)
  for (i in 1:length(x)) {
    FXi <- matrix(0, ncol = 4, nrow = 4)
    fxi <- matrix(0, ncol = 1, nrow = 4)
    
    xsi <- 1/(273.15 + x[i])
    alphai <- exp(gam1 + gam2 * xsi)
    betai <- exp(gam3 + gam4 * xsi)
    for (j in 1:l) {
      tij <- T_set[j]
      etaij <- exp(-betai * (tij)^alphai)
      fxij <- c(-alphai * betai * etaij * tij^alphai * log(tij),
                -alphai * betai * etaij * tij^alphai * log(tij) * xsi,
                -betai * etaij * tij^alphai,
                -betai * etaij * tij^alphai * xsi)
      FXi <- FXi + outer(fxij, fxij)
      fxi <- fxi + fxij
    } # end j
    SXXT <- SXXT + w[i] * FXi # FIM for theta
    SX <- SX + w[i] * fxi  
  } # end i
  
  I_theta2 <- ((1 + lamb^2 * a0) / (l * sigm^2)) * SXXT
  I_sig_theta <- (1/(l * sigm^2) * (b * lamb^3 * (1 - a0 * (1 + lamb^2))/(1 + lamb^2)^1.5 + lamb^2 * a1)) * SX
  I_lamb_theta <- -1/(l * sigm) * (b * a0 * lamb^2/(1 + lamb^2)^1.5 + lamb * a1) * SX
  
  I_sig2 <- (2 + lamb^2 * a2 + (a0 * b^2 * lamb^4 - b^2 * lamb^2)/(1 + lamb^2) -
               2 * b^2 * lamb^4/(1 + lamb^2)^2 - 2 * a1 * b * lamb^3/(1 + lamb^2)^0.5) / (sigm^2)
  I_sig_lamb <- -(lamb * a2 - a1 * b * lamb^4/(1 + lamb^2)^1.5 - a0 * b^2 * lamb^3/(1 + lamb^2)^2
                  + b^2 * lamb * (1 + 2 * lamb^2)/(1 + lamb^2)^3) / sigm
  I_lamb2 <- a2 + 2 * a1 * b * lamb/(1 + lamb^2)^1.5 + b^2 * (a0 * lamb^2 - 1)/(1 + lamb^2)^3
  
  FIM <- rbind(
    cbind(I_theta2, I_sig_theta, I_lamb_theta),
    cbind(t(I_sig_theta), I_sig2, I_sig_lamb),
    cbind(t(I_lamb_theta), I_sig_lamb, I_lamb2)
  )
  return(FIM)
}
# 


# Information matrix after reparameterization
# Input:
#   x      : the vector of design points (temperature values)
#   w      : the corresponding weight vector (same length as x)
#   param  : the new parameter vector c(alpha0, beta0, gamma2, gamma4, sigma, lambda)
FIM_SNNAM_reparam <- function(x, w, param) {
  # param is the new parameter vector: c(alpha0, beta0, gamma2, gamma4, sigma, lambda)
  alpha0 <- param[1]
  beta0  <- param[2]
  gamma2 <- param[3]
  gamma4 <- param[4]
  sigma  <- param[5]
  lambda <- param[6]
  
  # Compute the original parameters gamma1 and gamma3
  gamma1 <- log(alpha0) - gamma2 * xbar
  gamma3 <- log(beta0)  - gamma4 * xbar
  
  # Construct the original parameter vector
  theta <- c(gamma1, gamma2, gamma3, gamma4, sigma, lambda)
  
  # Call the original information matrix function (global variables l and T_set already exist)
  FIM_old <- FIM_SNNAM(x = x, w = w, param = theta)
  
  # Construct the Jacobian matrix J (6 x 6)
  J <- matrix(0, 6, 6)
  J[1, 1] <- 1 / alpha0
  J[1, 3] <- -xbar
  J[2, 3] <- 1
  J[3, 2] <- 1 / beta0
  J[3, 4] <- -xbar
  J[4, 4] <- 1
  J[5, 5] <- 1
  J[6, 6] <- 1
  
  # Transformation
  FIM_new <- t(J) %*% FIM_old %*% J
  return(FIM_new)
}


num_deriv_quantile <- function(lambda, p, h = 1e-5) {
  (qsn(p, alpha = lambda + h, , solver = "RFB") - 
     qsn(p, alpha = lambda - h, solver = "RFB")) / (2 * h)
}

# C-optimal design --------------------------------------------------------
# 
# C-optimal design criterion under the new parameterization
Copt_reparam <- function(x, w, param, fimfunc) {
  yf <- 0.6
  p <- 0.1
  alpha0 <- param[1]
  beta0  <- param[2]
  gamma2 <- param[3]
  gamma4 <- param[4]
  sigma  <- param[5]
  lambda <- param[6]
  
  # Use the use-condition x0 and the center xbar (need to be defined globally)
  x0 <- 1/(273.15 + 25)
  # Compute alpha0_use and beta0_use under the use condition
  alpha0_use <- alpha0 * exp(gamma2 * (x0 - xbar))
  beta0_use  <- beta0  * exp(gamma4 * (x0 - xbar))
  
  # Skew-normal quantile and related quantities
  zlambp <- qsn(p, alpha = lambda, solver = "RFB")
  b <- sqrt(2/pi)
  mu_z <- b * lambda / sqrt(1 + lambda^2)
  dmuz <- b * (1 + lambda^2)^(-1.5)
  
  # The quantile t_p
  tpt <- (-log(yf + sigma * mu_z - sigma * zlambp + 1e-12) / beta0_use)^(1/alpha0_use)
  
  # Construct the gradient vector Cv (order: alpha0, beta0, gamma2, gamma4, sigma, lambda)
  Cv <- matrix(0, nrow = 6, ncol = 1)
  Cv[1] <- -tpt * log(tpt) / alpha0
  Cv[2] <- -tpt / (alpha0_use * beta0)
  Cv[3] <- tpt * log(tpt) * (xbar - x0)
  Cv[4] <- tpt / alpha0_use * (xbar - x0)
  denom <- alpha0_use * log(yf + sigma * mu_z - sigma * zlambp + 1e-12) * 
    (yf + sigma * mu_z - sigma * zlambp + 1e-12)
  Cv[5] <- tpt * (mu_z - zlambp) / denom
  dz_lamp <- num_deriv_quantile(lambda, p = p, h = 1e-5)
  Cv[6] <- sigma * tpt * (dmuz - dz_lamp) / denom
  
  Cv_norm <- sqrt(sum(Cv^2))
  Cv <- Cv / Cv_norm 
  
  # Compute the criterion value
  fim <- fimfunc(x = x, w = w, param = param)
  
  # eps <- 1e-8 * max(eigen(fim)$values)   # adaptive
  eps <- 1e-8
  # Or a fixed small value such as 1e-6 (if the scale of fim is not large)
  fim_reg <- fim + diag(eps, nrow(fim))
  fim_inv <- solve(fim_reg)   # or use MASS::ginv
  Ccri <- t(Cv) %*% fim_inv %*% Cv
  Ccri <- Ccri * Cv_norm^2
  return(Ccri)
}

# Sensitivity function under the new parameterization
Copt_sens_reparam <- function(xi_x, x, w, param, fimfunc) {
  # Compute the information matrix of the current design and its inverse
  fim <- fimfunc(x = x, w = w, param = param)
  # eps <- 1e-8 * max(eigen(fim)$values)   # adaptive
  eps <- 1e-8 
  # Or a fixed small value such as 1e-6 (if the scale of fim is not large)
  fim_reg <- fim + diag(eps, nrow(fim))
  M_inv <- solve(fim_reg)   # or use MASS::ginv
  
  # Information matrix of a single-point design
  M_x <- fimfunc(x = xi_x, w = 1, param = param)
  
  # Same gradient computation as in Copt_reparam
  yf <- 0.6
  p <- 0.1
  alpha0 <- param[1]
  beta0  <- param[2]
  gamma2 <- param[3]
  gamma4 <- param[4]
  sigma  <- param[5]
  lambda <- param[6]
  
  x0 <- 1/(273.15 + 25)
  alpha0_use <- alpha0 * exp(gamma2 * (x0 - xbar))
  beta0_use  <- beta0  * exp(gamma4 * (x0 - xbar))
  
  zlambp <- qsn(p, alpha = lambda, solver = "RFB")
  b <- sqrt(2/pi)
  mu_z <- b * lambda / sqrt(1 + lambda^2)
  dmuz <- b * (1 + lambda^2)^(-1.5)
  
  tpt <- (-log(yf + sigma * mu_z - sigma * zlambp + 1e-12) / beta0_use)^(1/alpha0_use)
  
  Cv <- matrix(0, nrow = 6, ncol = 1)
  Cv[1] <- -tpt * log(tpt) / alpha0
  Cv[2] <- -tpt / (alpha0_use * beta0)
  Cv[3] <- tpt * log(tpt) * (xbar - x0)
  Cv[4] <- tpt / alpha0_use * (xbar - x0)
  denom <- alpha0_use * log(yf + sigma * mu_z - sigma * zlambp + 1e-12) * 
    (yf + sigma * mu_z - sigma * zlambp + 1e-12)
  Cv[5] <- tpt * (mu_z - zlambp) / denom
  dz_lamp <- num_deriv_quantile(lambda, p = p, h = 1e-5)
  Cv[6] <- sigma * tpt * (dmuz - dz_lamp) / denom
  
  Cv_norm <- sqrt(sum(Cv^2))
  Cv <- Cv / Cv_norm 
  # cat("Cv norm:", Cv_norm, "\n")
  # cat("cond(fim):", kappa(fim_reg), "\n")
  # Sensitivity function value
  Ccri_sen <- t(Cv) %*% M_inv %*% M_x %*% M_inv %*% Cv - t(Cv) %*% M_inv %*% Cv
  return(Ccri_sen)
}
# 
# Table 1 (left)
# ### case 1: different l and d 
# 
l <- 5  # the number of observation times: 5, 10, 15
d <- 6  # the distance between inspections: 6, 12, 24, 36
# equal inspection interval
T_set <- seq(d, l * d, d)
# given equal inspection interval
# T_set <- c(8, 25, 75, 130, 180)

Kv <- 6


resC <- locally(fimfunc = FIM_SNNAM_reparam,
                lx = Lowerx, ux = Upperx, inipars = phi,
                iter = 1000, k = Kv,
                crtfunc = Copt_reparam,
                sensfunc = Copt_sens_reparam,
                ICA.control = list(
                  ncount = 50,        
                  checkfreq = 200,     
                  stop_rule = "equivalence",
                  # stoptol = 0.99,      # can try 0.999 or higher
                  plot_cost = TRUE      # visualize the decrease of the cost function
                ))

print(d)
print(l)
merge_design_points(res = resC, Kv = Kv, threshold = 0.1, digits = 4)
print(resC$design$min_cost)

##



# Table 3 （design part）
# ### Case 2: fixed l and d, but different lambda 
# 
l <- 5  # the number of observation times: 5, 10, 15
d <- 36  # the distance between inspections: 6, 12, 24, 36
lamb_d <- -150  #-150, -15, -1.5, -0.15, 0.15, 1.5, 15，150

param0 <- c(-4.965450e+00,  1.370726e+03,  1.454751e+01, -5.956617e+03,
            2.584227e-02,  lamb_d)
b <- sqrt(2 / pi)
xbar <- mean(1/(273.5 + c(50, 65, 80)))
# xbar <- mean(1/(273.5 + 50))
phi <- c(exp(param0[1] + param0[2] * xbar), exp(param0[3] + param0[4] * xbar), 
         param0[2], param0[4], param0[5],  param0[6])


# equal inspection interval
T_set <- seq(d, l * d, d)
# given equal inspection interval
# T_set <- c(8, 25, 75, 130, 180)

Kv <- 6
resC <- locally(fimfunc = FIM_SNNAM_reparam,
                lx = Lowerx, ux = Upperx, inipars = phi,
                iter = 1000, k = Kv,
                crtfunc = Copt_reparam,
                sensfunc = Copt_sens_reparam,
                ICA.control = list(
                  ncount = 50,        # increase the number of countries (exploration ability)
                  checkfreq = 500,     # check the stopping rule more frequently
                  stop_rule = "equivalence",
                  # stoptol = 0.99,      # can try 0.999 or higher
                  plot_cost = TRUE      # visualize the decrease of the cost function
                ))

print(lamb_d)
merge_design_points(res = resC, Kv = Kv, threshold = 0.1, digits = 4)
print(resC$design$min_cost)




# 
# # Table 3 (efficiencies)
# To compute the efficiencies of optimal designs derived under the ADDT model with normal errors

# Lowerx <- c(25); Upperx <- c(80) # the design region


l <- 5  # the number of observation times:  5, 10, 15
d <- 36  # the distance between inspections: 6, 12, 24, 36

lamb_d <- 150  ##-150, -15, -1.5, -0.15, 0.15, 1.5, 15,150
param0 <- c(-4.965450e+00,  1.370726e+03,  1.454751e+01, -5.956617e+03,
            2.584227e-02,  lamb_d)
b <- sqrt(2 / pi)
xbar <- mean(1/(273.5 + c(50, 65, 80)))
# xbar <- mean(1/(273.5 + 50))
phi <- c(exp(param0[1] + param0[2] * xbar), exp(param0[3] + param0[4] * xbar), 
         param0[2], param0[4], param0[5],  param0[6])

# equal inspection interval
T_set <- seq(d, l * d, d)

xsn <- c(25,
         57.2755,
         80
) # supports of the optimal designs for ADDTM-SN
wsn <- c(0.5504,
         0.3131,
         0.1365) # weights of the optimal designs for ADDTM-SN

xn <- c(25,
        57.1752,
        80
)      # supports of the optimal designs for ADDTM-N
wn <- c(0.5631,
        0.3052,
        0.1317) # weights of the optimal designs for ADDTM-N
EFF <- Copt_reparam(x = xsn, w = wsn, param = phi, fimfunc = FIM_SNNAM_reparam)/
  Copt_reparam(x = xn, w = wn, param = phi, fimfunc = FIM_SNNAM_reparam)
print(EFF)


# library(xtable)
# xtable(t(aggregate(w ~ x, data = DW, FUN = sum)), digits = 4)
# print(d)
# print(l)
# print(resC$design$min_cost)




## D-opt


# D-optimal design --------------------------------------------------------
# 
Dopt <- function(x, w, param, fimfunc) {
  Dcri <- -log(det(fimfunc(x, w = w, param = param)))
  return(Dcri)
}
## the sensitivity function
# xi_x is a design that puts all its mass on x in the definition of the sensitivity function
# x is a vector of design points
Dopt_sens <- function(xi_x, x, w, param, fimfunc) {
  fim <- fimfunc(x, w = w, param = param)
  M_inv <- matrixcalc::matrix.inverse(fim)
  M_x <- fimfunc(xi_x, w = 1, param = param)
  
  Dcri_sen <- matrixcalc::matrix.trace(M_inv %*% M_x) - length(param)
  return(Dcri_sen)
}

# 
# Table 1 (right)
# ### case D1: different l and d 
# 
l <- 5  # the number of observation times: 5, 10, 15
d <- 6  # the distance between inspections: 6, 12, 24, 36
# equal inspection interval
T_set <- seq(d, l * d, d)
# given equal inspection interval
# T_set <- c(8, 25, 75, 130, 180)

Kv <- 6
resD <- locally(fimfunc = FIM_SNNAM_reparam,
                lx = Lowerx, ux = Upperx, inipars = phi,
                iter = 1000, k = Kv,
                crtfunc = Dopt,
                sensfunc = Dopt_sens,
                ICA.control = list(checkfreq = 500,
                                   plot_cost = FALSE,
                                   rseed = 1, 
                                   stop_rule = "equivalence",
                                   stoptol = .999))

print(d)
print(l)
merge_design_points(res = resD, Kv = Kv, threshold = 0.1, digits = 4)
print(resD$design$min_cost)


# Table 3(right part, design)
# ### DDcase 2: fixed l and d, but different lambda 
# 
l <- 5  # the number of observation times: 1, 5, 10, 15
d <- 36  # the distance between inspections: 6, 12, 24, 36
lamb_d <- -150  ##-150, -15, -1.5, -0.15, 0.15, 1.5, 15,150

param0 <- c(-4.965450e+00,  1.370726e+03,  1.454751e+01, -5.956617e+03,
            2.584227e-02,  lamb_d)
b <- sqrt(2 / pi)
xbar <- mean(1/(273.5 + c(50, 65, 80)))
# xbar <- mean(1/(273.5 + 50))
phi <- c(exp(param0[1] + param0[2] * xbar), exp(param0[3] + param0[4] * xbar), 
         param0[2], param0[4], param0[5],  param0[6])

# equal inspection interval
T_set <- seq(d, l * d, d)
# given equal inspection interval
# T_set <- c(8, 25, 75, 130, 180)

Kv <- 6
resD <- locally(fimfunc = FIM_SNNAM_reparam,
                lx = Lowerx, ux = Upperx, inipars = phi,
                iter = 1000, k = Kv,
                crtfunc = Dopt,
                sensfunc = Dopt_sens,
                ICA.control = list(checkfreq = 500,
                                   plot_cost = FALSE,
                                   rseed = 1, 
                                   stop_rule = "equivalence",
                                   stoptol = .999))

print(d)
print(l)
merge_design_points(res = resD, Kv = Kv, threshold = 0.1, digits = 4)
print(resD$design$min_cost)
# 


# 
# # Table 3 (right part, efficiencies) -----------------------------------------------------------------------
# 
# # To compute the efficiencies of optimal designs derived under the ADDT model with normal errors
# 
# # Lowerx <- c(25); Upperx <- c(80) # the design region
# 


# 
# ### DDcase 2: fixed l and d, but different lambda 
# 
l <- 5  # the number of observation times: 1, 5, 10, 15
d <- 36  # the distance between inspections: 6, 12, 24, 36
lamb_d <- 150  ## -150,-15, -1.5, -0.15, 0.15, 1.5, 15,150

param0 <- c(-4.965450e+00,  1.370726e+03,  1.454751e+01, -5.956617e+03,
            2.584227e-02,  lamb_d)
b <- sqrt(2 / pi)
xbar <- mean(1/(273.5 + c(50, 65, 80)))
# xbar <- mean(1/(273.5 + 50))
phi <- c(exp(param0[1] + param0[2] * xbar), exp(param0[3] + param0[4] * xbar), 
         param0[2], param0[4], param0[5],  param0[6])

# equal inspection interval
T_set <- seq(d, l * d, d)


xsn <- c(25,
         57.0626,
         80
) # supports of the optimal designs for ADDTM-SN
wsn <- c(0.3333,
         0.2650,
         0.4017) # weights of the optimal designs for ADDTM-SN

xn <- c(25,
        56.9561,
        80
)      # supports of the optimal designs for ADDTM-N
wn <- c(0.2828,
        0.2853,
        0.4319) # weights of the optimal designs for ADDTM-N
EFF <- exp((Dopt(x = xsn, w = wsn, param = phi, fimfunc = FIM_SNNAM_reparam) -
              Dopt(x = xn, w = wn, param = phi, fimfunc = FIM_SNNAM_reparam))/6)
print(EFF)