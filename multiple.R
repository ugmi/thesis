# A script to analyse and generate multiple data sets.

library(survival)
library(icenReg)

# Helper functions -------------------------------------------------------------

#' Get the date of the first day of the given date's month
#' 
#' @param ymd Date
#' @return Date
month_start <- function(ymd) as.Date(gsub("[0-9]{2}$", "01", ymd))


#' Get the date of the last day of the given date's month
#' 
#' @param ymd Date
#' @return Date
month_end <- function(ymd) month_start(ymd + 32) - 1


#' Stratified distribution function with no explanatory variables
#' 
#' @param v Random vector
#' @param x Treatment indicator
#' @return Vector of cumulative probabilities
pdist_strat <- function(v, x, pars) {
  lst <- sapply(names(pars), function(tr) 
    pweibull(v[x == tr], shape = pars[[tr]][[1]], scale = pars[[tr]][[2]]))
  return(do.call(c, lst))
}


#' Stratified hazard function with no explanatory variables
#' 
#' @param v Random vector
#' @param x Treatment indicator
#' @return Vector of hazards
hazard_strat <- function(v, x, pars) {
  lst <- sapply(names(pars), function(tr)
    dweibull(v[x == tr], pars[[tr]][[1]], pars[[tr]][[2]]) / (1 - pdist_strat(v[x == tr], tr, pars)))
  return(do.call(c, lst))
}


#' Generate time-to-event data from the Weibull distribution
#' 
#' @param n A named list providing the registry sizes
#' @param proportions A named list of named lists of the proportions of people receiving each treatment for each registry
#' @param pars A named list of Weibull dist. parameter values for each treatment
#' @return A data frame with four columns: survival time, probability, treatment type, registry
generate_weibull <- function(n, proportions, pars) {
  # Save names for later use
  regs <- names(proportions)
  treatments <- names(pars)
  
  # Calculation of sizes
  N <- sapply(regs, function(i)
    sapply(treatments, function(tr)
      round(n[[i]] * proportions[[i]][[tr]])))
  
  # Generate time-to-event data
  t.vec <- sapply(treatments, function(i) 
    rweibull(sum(N[i,]), shape = pars[[i]][[1]], scale = pars[[i]][[2]]))
  
  # Create the combined data set
  df <- do.call(rbind, lapply(treatments, function(tr)
    data.frame(
      "time" = t.vec[[tr]],
      "survival" = pweibull(t.vec[[tr]], pars[[tr]][[1]], pars[[tr]][[2]], lower.tail = TRUE),
      "treatment" = as.factor(tr),
      "registry" = as.factor(rep(regs, N[tr, ]))
    )))
  
  return(df)
}


#' Generate status, entry and exit dates given survival time in days
#' 
#' @param days Vector of survival time in days
#' @param start Study start date
#' @param duration Duration of follow-up from the study start in days
#' @param prob Probability of being lost-to-follow-up for each individual
#' @return A data frame with columns for status, entry and exit times
generate_dates <- function(days, start, duration, prob) {
  # Generate entry times and 
  N <- length(days)
  lost <- rbinom(N, 1, prob)
  t0 <- floor(runif(N, 0, duration - 1))
  
  # Create the data frame
  df <- data.frame(
    "entry" = start + t0,
    "status" = as.factor(ifelse(
      lost, "Lost to follow up", ifelse(t0 + days > duration, "Alive", "Dead")
    )))
  df["exit"] <- start + ifelse(
    df$status == "Alive",
    duration,
    t0 + ifelse(!lost, days, floor(pmin(duration - t0, days)*runif(N)))
  )
  return(df)
}


#' Generate lower and upper limits for interval-censored data
#' 
#' @param exit Date of exit
#' @param exact Binary indicator whether the date should be left exact
#' @param lower Vector of lower limits for the lower limit
#' @param upper Vector of upper limits for the upper limit
#' @param infinite Binary indicator if the upper limit is infinite
#' @return A data frame with columns for upper and lower interval limits
generate_intervals <- function(exit, exact, lower, upper, infinite) {
  L <- as.Date(ifelse(exact, exit, pmax(month_start(exit), lower)))
  U <- as.Date(ifelse(exact, exit, ifelse(infinite, Inf, pmin(month_end(exit), upper))))
  return(data.frame("lower" = L, "upper" = U))
}


#' Generate the data set
#' 
#' @param n A named list providing the registry sizes
#' @param proportions A named list of named lists of the proportions of people receiving each treatment for each registry
#' @param pars A named list of Weibull dist. parameter values for each treatment
#' @param start Study start date
#' @param duration Duration of follow-up from the study start in days
#' @param prob Probability of being lost-to-follow-up in each registry
#' @return A generated study data set
generate_df <- function(n, proportions, pars, start, duration, prob) {
  df <- generate_weibull(n, proportions, pars)
  df_dates <- generate_dates(floor(df$time), start, duration, prob[df$registry])
  df_intervals <- generate_intervals(
    df_dates$exit,
    exact = df$registry == "1" & df_dates$status == "Dead",
    lower = df_dates$entry,
    upper = start + duration,
    infinite = df_dates$status != "Dead"
  )
  return(do.call(cbind, list(df, df_dates, df_intervals)))
}


#' Fit a distribution using the proposed approach
#' 
#' @param df A data frame with study data
#' @param dist The name of the parametric distribution to fit
#' @return A fitted cumulative distribution function of the modelling variable
proposed_approach <- function(df, dist) {
  fy <- function(l, u, b, x) Fz(u, b, x) - Fz(l, b, x)
  logL <- function(b, l, u, x) -sum(log(fy(l, u, b, x)))
  
  # Create coarsening intervals
  l <- as.numeric(ifelse(df$lower - df$entry == 0, 0, df$lower - df$entry - 1))
  u <- ifelse(df$status == "Dead", df$upper - df$entry + 0.5, Inf)
  # Covariate vector
  x <- df$treatment
  
  if (dist %in% c("exp", "exponential")) {
    Fz <- function(z, b, x) {
      (x == "A")*pexp(z, b[1]) + (x == "B")*pexp(z, b[2]) + 
        (x == "C")*pexp(z, b[3])
    }
    mle <- nlm(logL, c(1/100, 1/100, 1/100), l, u, x, hessian = TRUE)
  } else if (dist == "gamma") {
    Fz <- function(z, b, x) {
      (x == "A")*pgamma(z, b[1], b[2]) + (x == "B")*pgamma(z, b[3], b[4]) + 
        (x == "C")*pgamma(z, b[5], b[6])
    }
    mle <- nlm(logL, c(1, 1/100, 1, 1/100, 1, 1/100), l, u, x, hessian = TRUE)
  } else if (dist == "weibull") {
    Fz <- function(z, b, x) {
      (x == "A")*pweibull(z, exp(b[1]), exp(b[2])) + 
        (x == "B")*pweibull(z, exp(b[3]), exp(b[4])) + 
        (x == "C")*pweibull(z, exp(b[5]), exp(b[6]))
    }
    mle <- nlm(logL, c(-0.2, 6, -0.2, 6, -0.2, 6), l, u, x, hessian = TRUE)
  }
  stopifnot("The optimizer failed." = mle$minimum > 1)
  f <- function(z, x) Fz(z, mle$estimate, x)
  attr(f, "AIC") <- 2*mle$minimum + 2*length(mle$estimate)
  attr(f, "mle") <- mle$estimate
  return(f)
}


#' Fit a distribution using the generalised version of the proposed approach
#' 
#' @param df A data frame with study data
#' @return A fitted cumulative distribution function of the modelling variable
generalized_approach <- function(df) {
  # Define the pre-transform
  s <- function(z) qlogis(pweibull(z, 0.8, 500))
  # Create the coarsening limits
  sl <- s(as.numeric(ifelse(df$lower - df$entry == 0, 0, df$lower - df$entry - 1)))
  su <- s(ifelse(df$status == "Dead", df$upper - df$entry + 0.5, Inf))
  # Define the modelling distribution and polynomial approximation
  Fz <- function(z, b, x) plogis(g(z, b, x))
  g <- function(z, b, x) {
    (b[1] + b[2]*(x == "B") + b[3]*(x == "C")) + 
      z*(b[4] + b[5]*(x == "B") + b[6]*(x == "C"))
  }
  fy <- function(l, u, b, x) Fz(u, b, x) - Fz(l, b, x)
  logL <- function(b, l, u, x) -sum(log(fy(l, u, b, x)))
  # Estimate the parameters
  mle <- nlm(logL, rep(1, 6), u=su, l=sl, x=df$treatment)
  stopifnot("The optimizer failed." = mle$minimum > 1)
  # Save the fitted function
  f <- function(z, x) plogis(g(s(z), mle$estimate, x))
  attr(f, "AIC") <- 2*mle$minimum + 2*length(mle$estimate)
  return(f)
}


# Data generation --------------------------------------------------------------

set.seed(505)

estimate_pars <- function(df, parameters, days) {
  # Proposed approach
  pdist_weibull <- proposed_approach(df, "weibull")
  pdist_gamma <- proposed_approach(df, "gamma")
  
  # Generalized proposed approach
  pdist_ge <- generalized_approach(df)
  
  # Accelerated failure time model
  df["days.lower"] <- as.numeric(df$lower - df$entry)
  df["days.upper"] <- as.numeric(df$upper - df$entry)
  df["days.lower"] <- ifelse(df$days.lower == 0 & df$days.upper == 0, 0.5, df$days.lower)
  df["days.upper"] <- ifelse(df$days.upper == 0, 0.5, df$days.upper)
  fit.par <- ic_par(cbind(days.lower, days.upper) ~ treatment, data = df, 
                    model = "aft", dist = "weibull", weights = NULL)
  weibull_aft <- function(t, x) {
    coeffs <- exp(fit.par$coefficients)
    (x == "A")*pweibull(t, shape = coeffs[[1]], scale = coeffs[[2]]) +
      (x == "B")*pweibull(t, coeffs[[1]], coeffs[[2]]/exp(log(coeffs[[3]])/coeffs[[1]])) +
      (x == "C")*pweibull(t, coeffs[[1]], coeffs[[2]]/exp(log(coeffs[[4]])/coeffs[[1]]))
  }
  
  # Imputation
  df["mid"] <- as.Date(ifelse(
    df$upper == df$lower | df$status != "Dead",
    df$lower,
    df$lower + floor((df$upper - df$lower)/2)))
  df["days.mid"] <- as.numeric(df$mid - df$entry)
  
  # Parametric Weibull model
  fit.weibull <- survreg(Surv(days.mid, status == "Dead") ~ strata(treatment), 
                         data = df[df$days.mid != 0,], dist = "weibull")
  weibull_parametric <- function(t, x) {
    #  rweibull shape = 1/(survreg scale) = 1/exp(survreg Log(scale))
    #  rweibull scale = exp(survreg intercept)
    coeffs <- exp(fit.weibull$icoef)
    (x == "A")*pweibull(t, shape = 1/coeffs[[2]], scale = coeffs[[1]]) +
      (x == "B")*pweibull(t, 1/coeffs[[3]], coeffs[[1]]) +
      (x == "C")*pweibull(t, 1/coeffs[[4]], coeffs[[1]])
  }
  
  # Create tables
  AIC <- c("weibull" = attr(pdist_weibull, "AIC"), 
           "gamma" = attr(pdist_gamma, "AIC"), 
           "generalised" = attr(pdist_ge, "AIC"), 
           "imp.weibull" = -2*fit.weibull$loglik[1] + 2*fit.weibull$df,
           "aft.weibull" = -2*fit.par$llk + 2*length(fit.par$coefficients))
  
  get_table <- function(days, tr, parameters) {
    stopifnot(tr %in% names(parameters))
    cumulative <- pdist_strat(days, tr, parameters)
    names(cumulative) <- days
    tbl.surv <- data.frame(
      "weibull" = cumulative - pdist_weibull(days, tr),
      "gamma" = cumulative - pdist_gamma(days, tr),
      "generalised" = cumulative - pdist_ge(days, tr),
      "imp.weibull" = cumulative - weibull_parametric(days, tr),
      "aft.weibull" = cumulative - weibull_aft(days, tr),
      row.names = days
    )
    attr(tbl.surv, "true") <- 1 - cumulative
    return(tbl.surv)
  }
  
  tbl <- lapply(names(parameters), function(tr) get_table(days, tr, parameters))
  names(tbl) <- names(parameters)
  attr(tbl, "AIC") <- AIC
  return(tbl)
}

# Description of the setting
sizes <- c("1" = 1200, "2" = 900)
treat.prop <- list("1" = list("A" = 0.3, "B" = 0.2, "C" = 0.5), 
                   "2" = list("A" = 0.47, "B" = 0.36, "C" = 0.17))
parameters <- list("A" = c(0.8, 500), "B" = c(0.9, 500), "C" = c(0.7, 500))
lost.prop <- c("1" = 0.014, "2" = 0.018)
duration <- 1095
start <- as.Date("2001-01-01")

# Find average bias
days <- as.numeric(month_end(start + 1 + c(0, 1, 2, 3, 4, 5, 11, 17, 23, 35, 47, 59)*31) - start)
ndays <- length(days)
tbl.empty <- data.frame(
  "weibull" = rep(0, ndays), 
  "gamma" = rep(0, ndays), 
  "generalised" = rep(0, ndays), 
  "imp.weibull" = rep(0, ndays), 
  "aft.weibull" = rep(0, ndays)
)
tbl.surv <- list("A" = tbl.empty, "B" = tbl.empty, "C" = tbl.empty, "AIC" = rep(0, 5))
n.iter <- 100
# Add a progress bar to help estimate how long the for loop would take
pb <- txtProgressBar(min = 0, max = n.iter, style = 3, width = 50, char = "=")
for(i in 1:n.iter) {
  df.obs <- generate_df(sizes, treat.prop, parameters, start, duration, lost.prop)
  x <- estimate_pars(df.obs, parameters)
  tbl.surv$A <- tbl.surv$A + x$A/100
  tbl.surv$B <- tbl.surv$B + x$B/100
  tbl.surv$C <- tbl.surv$C + x$C/100
  tbl.surv$AIC <- tbl.surv$AIC + attr(x, "AIC")/100
  setTxtProgressBar(pb, i)
}






