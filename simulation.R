library(survival)
library(flexsurv)
library(icenReg)

# Generating survival times ----------------------------------------------------

#' Generate time-to-event data from a parametric distribution
#' 
#' @param n Number of observations
#' @param dist Name of the distribution
#' @param pars A named vector of dist. parameter values
#' @return A data frame with one column: survival time
generate_times <- function(n, dist, pars) {
  if(dist == "weibull") {
    t <- rweibull(n, pars[["shape"]], pars[["scale"]])
  } else if(dist == "gamma") {
    t <- rgamma(n, pars[["shape"]], pars[["rate"]])
  } else if(dist == "exponential") {
    t <- rexp(n, pars[["rate"]])
  }
  return(data.frame("time" = t))
}


#' Generate status, entry and exit dates given survival time in days
#' 
#' @param days Vector of survival times in days
#' @param duration Duration of enrollment in days
#' @param followup Duration of follow-up from end of enrollment in days
#' @return A data frame with two columns: status and survival days at exit
generate_right_censoring <- function(days, duration, followup) {
  # Generate entry times
  t0 <- if(duration > 0) ceiling(runif(length(days), 0, duration)) else 0
  alive <- days > t0 + followup
  
  # Create the data frame
  df <- data.frame(
    "status" = as.factor(ifelse(alive, "Alive", "Dead")),
    "exit" = ifelse(alive, t0 + followup, days)
  )
  
  return(df)
}


#' Generate lower and upper limits for interval-censored data
#' 
#' @param exit Survival time at exit in days
#' @param width Interval size in days for each individual
#' @param infinite Binary indicator if the upper limit is infinite
#' @return A data frame with two columns: lower and upper censoring limit
generate_intervals <- function(exit, infinite, width) {
  L <- ifelse(
    width == 0 | infinite, 
    exit, 
    exit - floor(runif(length(infinite), 1, pmin(exit, width)))
  )
  U <- ifelse(infinite, Inf, L + width)
  obs <- ifelse(!infinite, U, exit)
  mid <- ifelse(infinite, obs, obs - width/2)
  df <- data.frame("lower" = L, "upper" = U, "obs" = obs, "mid" = mid)
  return(df)
}


#' Generate the data set
#' 
#' @param n Number of observations
#' @param dist Name of the distribution
#' @param pars A named list of dist. parameter values
#' @param duration Duration of enrollment in days
#' @param followup Duration of follow-up from end of enrollment in days
#' @param width Vector of interval widths
#' @return A list containing the true survival time and 
#' generated data frames with interval-censored times
generate_df <- function(n, dist, pars, duration, followup, width) {
  df <- generate_times(n, dist, pars)
  df.rcensor <- generate_right_censoring(ceiling(df$time), duration, followup)
  alive <- df.rcensor$status == "Alive"
  df.intervals <- lapply(width, function(x)
    cbind(df.rcensor, generate_intervals(df.rcensor$exit, alive, x)))
  
  return(list(df, df.intervals))
}


# Estimation -------------------------------------------------------------------

#' Estimated survival for different time points and distribution parameters
#' 
#' @param D Time points for which to estimate survival
#' @param df Data frame containing interval-censored survival times
#' @param dist Parametric distribution to use for estimation
#' @return A list of arrays of survival and parameter estimates
get_estimates <- function(D, df, dist) {
  KM.mid <- survfit(Surv(mid, status == "Dead") ~ 1, data = df)
  KM.obs <- survfit(Surv(obs, status == "Dead") ~ 1, data = df)
  NPMLE <- ic_np(df[, c("lower", "upper")])
  scurv <- getSCurves(NPMLE)
  
  mid <- flexsurvreg(Surv(mid, status == "Dead") ~ 1, data = df, dist = dist)
  ign <- flexsurvreg(Surv(obs, status == "Dead") ~ 1, data = df, dist = dist)
  ic <- flexsurvreg(Surv(lower, upper, type = "interval2") ~ 1, data = df, dist = dist)
  
  pD <- array(c(
    summary(KM.mid, times = D)[[6]],
    summary(KM.obs, times = D)[[6]],
    sapply(D, function(d) scurv$S_curves$baseline[scurv$Tbull_ints[, 2] >= d][1]),
    summary(mid, type = "survival", t = D, ci = FALSE, se = FALSE)[[1]][[2]],
    summary(ign, type = "survival", t = D, ci = FALSE, se = FALSE)[[1]][[2]],
    summary(ic, type = "survival", t = D, ci = FALSE, se = FALSE)[[1]][[2]]
  ),
  dim = c(length(D), 6), 
  dimnames = list(D, c("KM.MID", "KM.IGN", "NPMLE", "MID", "IGN", "IC"))
  )
  
  pars <- array(exp(c(coef(mid), coef(ign), coef(ic))), 
                dim = c(length(mid$dlist$pars), 3),
                dimnames = list(mid$dlist$pars, c("MID", "IGN", "IC")))
  
  return(list(pD, pars))
}


repeat_estimation <- function(iter, days, df.pars, dist.pars) {
  est.prob <- array(
    NA,
    dim = c(iter, length(days), 6, 3),
    dimnames = list(
      1:iter,
      days,
      c("KM.MID", "KM.IGN", "NPMLE", "MID", "IGN", "IC"),
      1:3
    )
  )
  est.pars <- array(
    NA,
    dim = c(iter, length(dist.pars), 3, 3),
    dimnames = list(
      1:iter, 
      dist.pars, 
      c("MID", "IGN", "IC"), 
      1:3
    )
  )
  
  for(i in 1:iter) {
    df.lst <- do.call(generate_df, df.pars)
    for(j in 1:3) {
      est <- get_estimates(days, df.lst[[2]][[j]], df.pars[["dist"]])
      est.prob[i,,,j] <- est[[1]]
      est.pars[i,,,j] <- est[[2]]
    }
  }
  return(list("prob" = est.prob, "pars" = est.pars))
}


average_bias <- function(thetas, true.vals) {
  bias <- apply(thetas, c(1, 3, 4), function(v) v - true.vals)
  avg.bias <- apply(bias, c(1, 3, 4), mean)
  return(avg.bias)
}


## Plotting --------------------------------------------------------------------

plot_bias <- function(bias, cols, width) {
  n.mod <- dim(bias)[2]
  labs <- dimnames(bias)
  days <- as.numeric(labs[[1]])
  ymax <- max(max(bias), max(-bias)) + 0.0005
  pchs <- c(20, 8, 21, 22, 24, 25)
  plot(days, bias[,1], ylim = c(-ymax, ymax), pch = 20, col = cols[1],
       xaxt = "n", xlab = "Days", ylab = "Bias", main = paste0("width=", width))
  axis(1, at = days, labels = TRUE)
  for(i in 2:n.mod) points(days, bias[,i], pch = pchs[i], col = cols[i])
  lines(c(0, 2000), c(0, 0), col = "grey", lwd = 0.5)
  legend("bottom", pch = pchs[1:n.mod], bty = "n", col = cols[1:n.mod],
         legend = labs[[2]], ncol = 3)
}

### FILE END