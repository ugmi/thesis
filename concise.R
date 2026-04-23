library(survival)
library(icenReg)

set.seed(1984)

# Helper functions -------------------------------------------------------------

#' Generate lower and upper limits for interval-censored data
#' 
#' @param exit Survival time at exit in days
#' @param width Interval size in days for each individual
#' @return A data frame with two columns: lower and upper censoring limit
generate_intervals <- function(exit, width) {
  offset <- floor(runif(length(exit), 1, width))
  L <- pmax(0, exit - offset)
  U <- ifelse(L == 0 & width != 1, exit + width - offset, L + width)
  df <- data.frame("lower" = L, "upper" = U, "obs" = U, "mid" = U - (U - L)/2)
  return(df)
}


#' Generate the data set
#' 
#' @param n Number of observations
#' @param pars A named list of dist. parameter values
#' @param width Vector of interval widths
#' @return A list of generated data frames with interval-censored times
generate_df <- function(n, pars, width) {
  exit <- ceiling(rweibull(n, pars[["shape"]], pars[["scale"]]))
  df.intervals <- lapply(width, function(x)
    cbind(exit, generate_intervals(exit, x)))
  return(df.intervals)
}


# Transfrom parameters from survreg to pweibull parametrisation
transform_pars <- function(surv.pars) {
  return(c(1/exp(surv.pars[[2]]), exp(surv.pars[[1]])))
}


#' Estimate survival for different time points and distribution parameters
#' 
#' @param D Time points for which to estimate survival
#' @param df Data frame containing interval-censored survival times
#' @return A list of arrays of survival and parameter estimates
get_estimates <- function(D, df) {
  status <- rep(1, nrow(df))
  
  # Fit non-parametric models
  KM.mid <- survfit(Surv(mid, status) ~ 1, data = df)
  KM.obs <- survfit(Surv(obs, status) ~ 1, data = df)
  NPMLE <- ic_np(df[, c("lower", "upper")])
  scurv <- getSCurves(NPMLE)
  
  # Fit parametric models
  mid <- survreg(Surv(mid, status) ~ 1, data = df)
  ign <- survreg(Surv(obs, status) ~ 1, data = df)
  ic <- survreg(Surv(lower, upper, type = "interval2") ~ 1, data = df)
  
  # Estimates from model fit differ from base R parametrisation, 
  # thus we need to transform them.
  mid.pars <- transform_pars(mid$icoef)
  ign.pars <- transform_pars(ign$icoef)
  ic.pars <- transform_pars(ic$icoef)
  
  # Estimate and save probabilities into an array
  pD <- array(c(
    summary(KM.mid, times = D)[[6]],
    summary(KM.obs, times = D)[[6]],
    sapply(D, function(d) 
      scurv$S_curves$baseline[scurv$Tbull_ints[, 2] >= d][1]),
    1 - pweibull(D, mid.pars[[1]], mid.pars[[2]]),
    1 - pweibull(D, ign.pars[[1]], ign.pars[[2]]),
    1 - pweibull(D, ic.pars[[1]], ic.pars[[2]])
  ),
  dim = c(length(D), 6), 
  dimnames = list(D, c("KM.MID", "KM.IGN", "NPMLE", "MID", "IGN", "IC"))
  )
  
  # Save parameters into an array.
  pars <- array(c(mid.pars, ign.pars, ic.pars), 
                dim = c(2, 3),
                dimnames = list(c("shape", "scale"), c("MID", "IGN", "IC")))
  
  return(list(pD, pars))
}


#' Repeated simulation and parameter estimation
#' 
#' @param iter Number of simulated data sets
#' @param days Vector of survival times in days to estimate
#' @param df.pars Named list of parameters to generate the data sets
#' @return List of arrays containing est. survival and dist. parameter values
repeat_estimation <- function(iter, days, df.pars) {
  # Empty array to save survival probability estimates
  model.names <- c("KM.MID", "KM.IGN", "NPMLE", "MID", "IGN", "IC")
  est.prob <- array(
    NA, 
    dim = c(iter, length(days), 6, 3),
    dimnames = list(1:iter, days, model.names, 1:3)
  )
  
  # Empty array to save distribution parameter estimates
  est.pars <- array(
    NA, 
    dim = c(iter, 2, 3, 3),
    dimnames = list(1:iter, c("shape", "scale"), c("MID", "IGN", "IC"), 1:3)
  )
  
  # Repeat simulation and estimation `iter` times
  for(i in 1:iter) {
    df.lst <- do.call(generate_df, df.pars)
    for(j in 1:3) {
      est <- get_estimates(days, df.lst[[j]])
      est.prob[i,,,j] <- est[[1]]
      est.pars[i,,,j] <- est[[2]]
    }
  }
  
  return(list("prob" = est.prob, "pars" = est.pars))
}


#' Calculate average bias
#' 
#' @param thetas Array of estimates
#' @param true.vals Vector of true parameter values
#' @return Array of average bias estimates
average_bias <- function(thetas, true.vals) {
  bias <- apply(thetas, c(1, 3, 4), function(v) v - true.vals)
  avg.bias <- apply(bias, c(1, 3, 4), mean)
  return(avg.bias)
}


#' Plot the bias of survival probabilities
#' 
#' @param bias Array of bias estimates
#' @param title Character vector, title for the graph
#' @return NULL
plot_bias <- function(bias, title) {
  # Extract needed info from given objects
  labs <- dimnames(bias)
  days <- as.numeric(labs[[1]])
  ymax <- max(max(bias), max(-bias)) + 0.0005
  pchs <- c(20, 8, 21, 22, 24, 25)
  cols <- c("darkorange", "red", "brown", "forestgreen", "royalblue", "cyan3")
  
  # Make the plot
  plot(days, bias[, 1], pch = 20, col = cols[1], ylim = c(-ymax, ymax),
       xaxt = "n", xlab = "t, in days", ylab = "Bias", main = title)
  axis(1, at = days, labels = TRUE)
  for(i in 2:6) 
    points(days, bias[, i], pch = pchs[i], col = cols[i])
  lines(c(0, 2000), c(0, 0), col = "grey", lwd = 0.5)
  legend("bottom", pch = pchs, bty = "n", col = cols, legend = labs[[2]], ncol = 3)
}


# Global parameters ------------------------------------------------------------
widths <- c("1" = 1, "2" = 30, "3" = 365)  # interval widths

# Parameters for plotting
title1 <- "Bias for est. S(t), \ninterval width = 1 day"
title30 <- "Bias for est. S(t), \ninterval width = 30 days"
title365 <- "Bias for est. S(t), \ninterval width = 365 days"

# Short-term survival ----------------------------------------------------------
days1 <- c(30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 333, 365)
pars1 <- c("shape" = 0.7, "scale" = 200)

df.pars1 <- list("n" = 2000, "pars" = pars1, "width" = widths)
est1 <- repeat_estimation(1000, days1, df.pars1)

par.bias1 <- average_bias(est1[["pars"]], pars1)
prob.bias1 <- average_bias(est1[["prob"]], 1 - pweibull(days1, 0.7, 200))

plot_bias(prob.bias1[,,1], title1)
plot_bias(prob.bias1[,,2], title30)
plot_bias(prob.bias1[,,3], title365)

par.bias1[,,1]
par.bias1[,,2]
par.bias1[,,3]

# Long-term survival -----------------------------------------------------------
days5 <- c(365, 730, 1095, 1461, 1826)
pars5 <- c("shape" = 0.6, "scale" = 2000)

df.pars5 <- list("n" = 2000, "pars" = pars5, "width" = widths)
est5 <- repeat_estimation(1000, days5, df.pars5)

par.bias5 <- average_bias(est5[["pars"]], pars5)
prob.bias5 <- average_bias(est5[["prob"]], 1 - pweibull(days5, 0.6, 2000))

plot_bias(prob.bias5[,,1], title1)
plot_bias(prob.bias5[,,2], title30)
plot_bias(prob.bias5[,,3], title365)

par.bias5[,,1]
par.bias5[,,2]
par.bias5[,,3]

