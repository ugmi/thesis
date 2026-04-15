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
  } else if(dist == "lognormal") {
    t <- rlnorm(n, pars[["meanlog"]], pars[["sdlog"]])
  } else if(dist == "loglogistic") {
    t <- rllogis(n, pars[["shape"]], pars[["scale"]])
  } else if(dist == "gompertz") {
    t <- rgompertz(n, pars[["shape"]], pars[["rate"]])
  } else {
    stop("distribution not available")
  }
  return(t)
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
  # Create interval limits
  offset <- floor(runif(length(infinite), 1, width))
  L <- ifelse(width == 0 | infinite, exit, pmax(0, exit - offset))
  U <- ifelse(infinite,
              Inf,
              ifelse(L == 0 & width != 1, exit + width - offset, L + width))
  
  # Derive additional times
  obs <- ifelse(!infinite, U, exit)  # "observed" value
  mid <- ifelse(infinite, obs, obs - (U - L)/2)  # midpoint imputation
  
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
  times <- generate_times(n, dist, pars)
  df.rcensor <- generate_right_censoring(ceiling(times), duration, followup)
  alive <- df.rcensor$status == "Alive"
  df.intervals <- lapply(width, function(x)
    cbind(df.rcensor, generate_intervals(df.rcensor$exit, alive, x)))
  
  return(list(times, df.intervals))
}


# Estimation -------------------------------------------------------------------

#' Estimate survival for different time points and distribution parameters
#' 
#' @param D Time points for which to estimate survival
#' @param df Data frame containing interval-censored survival times
#' @param dist Parametric distribution to use for estimation
#' @return A list of arrays of survival and parameter estimates
get_estimates <- function(D, df, dist) {
  # Fit non-parametric models
  KM.mid <- survfit(Surv(mid, status == "Dead") ~ 1, data = df)
  KM.obs <- survfit(Surv(obs, status == "Dead") ~ 1, data = df)
  NPMLE <- ic_np(df[, c("lower", "upper")])
  scurv <- getSCurves(NPMLE)
  
  # Fit parametric models
  mid <- flexsurvreg(Surv(mid, status == "Dead") ~ 1, data = df, dist = dist)
  ign <- flexsurvreg(Surv(obs, status == "Dead") ~ 1, data = df, dist = dist)
  ic <- flexsurvreg(Surv(lower, upper, type = "interval2") ~ 1, 
                    data = df, dist = dist)
  
  # Estimate and save probabilities into an array
  pD <- array(c(
    summary(KM.mid, times = D)[[6]],
    summary(KM.obs, times = D)[[6]],
    sapply(D, function(d) 
      scurv$S_curves$baseline[scurv$Tbull_ints[, 2] >= d][1]),
    summary(mid, type = "survival", t = D, ci = FALSE, se = FALSE)[[1]][[2]],
    summary(ign, type = "survival", t = D, ci = FALSE, se = FALSE)[[1]][[2]],
    summary(ic, type = "survival", t = D, ci = FALSE, se = FALSE)[[1]][[2]]
  ),
  dim = c(length(D), 6), 
  dimnames = list(D, c("KM.MID", "KM.IGN", "NPMLE", "MID", "IGN", "IC"))
  )
  
  # Estimates from model fit differ from base R parametrisation, 
  # thus we need to transform them before saving into an array.
  f1 <- mid$dlist$inv.transforms[[1]]
  f2 <- mid$dlist$inv.transforms[[2]]
  pars <- array(c(f1(coef(mid)[[1]]), f2(coef(mid)[[2]]),
                  f1(coef(ign)[[1]]), f2(coef(ign)[[2]]),
                  f1(coef(ic)[[1]]), f2(coef(ic)[[2]])), 
                dim = c(2, 3),
                dimnames = list(mid$dlist$pars, c("MID", "IGN", "IC")))
  
  return(list(pD, pars))
}


#' Repeated simulation and parameter estimation
#' 
#' @param iter Number of simulated data sets
#' @param days Vector of survival times in days to estimate
#' @param df.pars Named list of parameters to generate the data sets
#' @param dist.pars Named vector of distribution parameters
#' @return List of arrays containing est. survival and dist. parameter values
repeat_estimation <- function(iter, days, df.pars, dist.pars) {
  # Empty array to save survival probability estimates
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
  
  # Empty array to save distribution parameter estimates
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
  
  # Repeat simulation and estimation `iter` times
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


## Plotting --------------------------------------------------------------------

#' Plot the bias of survival probabilities
#' 
#' @param bias Array of bias estimates
#' @param cols Vector of colors corresponding to each model to use for plotting
#' @param title Character vector, title for the graph
#' @return NULL
plot_bias <- function(bias, cols, title) {
  # Extract needed info from given objects
  n.mod <- dim(bias)[2]
  labs <- dimnames(bias)
  days <- as.numeric(labs[[1]])
  ymax <- max(max(bias), max(-bias)) + 0.0005
  pchs <- c(20, 8, 21, 22, 24, 25)
  
  # Make the plot
  plot(
    days,
    bias[, 1],
    ylim = c(-ymax, ymax),
    pch = 20,
    col = cols[1],
    xaxt = "n",
    xlab = "t, in days",
    ylab = "Bias",
    main = title
  )
  axis(1, at = days, labels = TRUE)
  for(i in 2:n.mod) 
    points(days, bias[,i], pch = pchs[i], col = cols[i])
  lines(c(0, 2000), c(0, 0), col = "grey", lwd = 0.5)
  legend(
    "bottom",
    pch = pchs[1:n.mod],
    bty = "n",
    col = cols[1:n.mod],
    legend = labs[[2]],
    ncol = 3
  )
}

### FILE END