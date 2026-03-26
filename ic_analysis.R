# Code for examining some of the properties of interval censoring.

library(survival)
set.seed(1984)

## DATA GENERATION FUNCTIONS ---------------------------------------------------
#' Generate time-to-event data from the Weibull distribution
#' 
#' @param n A named vector providing the register sizes
#' @param proportions A named list of named lists of the proportions of people 
#' receiving each treatment for each register
#' @param pars A named list of Weibull dist. parameter values for each treatment
#' @return A data frame with four columns: 
#' survival time, survival probability, treatment type, register
generate_weibull <- function(n, proportions, pars) {
  # Save names for later use
  regs <- names(proportions)
  treatments <- names(pars)
  
  if(length(treatments) > 1) {
    # Calculation of sizes
    k <- sapply(regs, function(i)
      sapply(treatments, function(tr) 
        round(n[[i]] * proportions[[i]][[tr]])))
    
    # Generate time-to-event data
    t.vec <- sapply(treatments, function(i) 
      rweibull(sum(k[i,]), shape = pars[[i]][[1]], scale = pars[[i]][[2]]))
    
    # Create the combined data set
    df <- do.call(rbind, lapply(treatments, function(tr)
      data.frame(
        "time" = t.vec[,tr],
        "survival" = 1 - pweibull(t.vec[,tr], pars[[tr]][[1]], pars[[tr]][[2]]),
        "treatment" = as.factor(tr),
        "register" = as.factor(rep(regs, k[tr, ]))
      )))
  } else {
    t.vec <- rweibull(sum(n), shape = pars[[1]][[1]], scale = pars[[1]][[2]])
    df <- data.frame(
      "time" = t.vec,
      "survival" = 1 - pweibull(t.vec, pars[[1]][[1]], pars[[1]][[2]]),
      "treatment" = treatments,
      "register" = rep(names(n), n)
    )
  }
  return(df)
}


#' Generate status, entry and exit dates given survival time in days
#' 
#' @param days Vector of survival times in days
#' @param duration Duration of enrollment in days
#' @param followup Duration of follow-up from end of enrollment in days
#' @param prob Probability of being lost-to-follow-up for each individual
#' @return A data frame with two columns: status and survival days until exit
generate_censoring <- function(days, duration, followup, prob) {
  # Generate entry times
  K <- length(days)
  lost <- rbinom(K, 1, prob)
  t0 <- if(duration > 0) ceiling(runif(K, 0, duration)) else 0
  
  # Create the data frame
  df <- data.frame(
    "status" = as.factor(ifelse(
      lost, "Lost", ifelse(days > t0 + followup, "Alive", "Dead")
    )))
  df["exit"] <- ifelse(
    df$status == "Alive",
    t0 + followup,
    ifelse(!lost, days, ceiling(runif(K, 0, pmin(t0 + followup, days))))
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
  return(data.frame("lower" = L, "upper" = U))
}


#' Generate the data set
#' 
#' @param n A named vector providing the register sizes
#' @param proportions A named list of named lists of the proportions of people 
#' receiving each treatment for each register
#' @param pars A named list of Weibull dist. parameter values for each treatment
#' @param duration Duration of enrollment in days
#' @param followup Duration of follow-up from end of enrollment in days
#' @param prob A named vector of probabilities of being lost-to-follow-up 
#' in each register
#' @return The generated data frame
generate_df <- function(n, proportions, pars, duration, followup, prob, cens) {
  df <- generate_weibull(n, proportions, pars)
  df.censoring <- generate_censoring(
    ceiling(df$time), duration, followup, prob[df$register])
  widths <- as.numeric(cens[df$register])
  not.dead <- df.censoring$status != "Dead"
  df.intervals <- generate_intervals(df.censoring$exit, not.dead, widths)
  df.observed <- data.frame("obs" = ifelse(
    !not.dead, df.intervals$upper, df.censoring$exit))
  df.imputed <- data.frame("mid" = ifelse(
    not.dead, df.observed$obs, df.observed$obs - widths/2))
  return(do.call(
    cbind, list(df, df.censoring, df.intervals, df.observed, df.imputed)))
}

## MODEL FITTING ---------------------------------------------------------------
#' Fit a distribution using the proposed approach
#' 
#' @param df A data frame with study data
#' @param dist The name of the parametric distribution to fit
#' @return A fitted cumulative distribution function of the modelling variable
proposed_approach <- function(df, dist) {
  logL <- function(b, l, u) -sum(log(Fz(u, b) - Fz(l, b)))
  
  l <- df$lower
  u <- df$upper
  
  if (dist %in% c("exp", "exponential")) {
    Fz <- function(z, b) pexp(z, exp(b[1]))
    mle <- nlm(logL, c(-4.5), l, u)
  } else if (dist == "gamma") {
    Fz <- function(z, b) pgamma(z, b[1], b[2])
    mle <- nlm(logL, c(1, 1/100), l, u)
  } else if (dist == "weibull") {
    Fz <- function(z, b) pweibull(z, exp(b[1]), exp(b[2]))
    mle <- nlm(logL, c(-0.2, 6), l, u)
  }
  stopifnot("The optimizer failed." = mle$minimum > 1)
  f <- function(z) Fz(z, mle$estimate)
  attr(f, "AIC") <- 2*mle$minimum + 2*length(mle$estimate)
  attr(f, "mle") <- mle$estimate
  return(f)
}


#' Fit a distribution using the generalized version of the proposed approach
#' 
#' @param df A data frame with study data
#' @return A fitted cumulative distribution function of the modelling variable
generalized_approach <- function(df, pre.transf) {
  # Define the pre-transform
  s <- function(z) qlogis(pre.transf(z))
  # Create the coarsening limits
  sl <- s(df$lower)
  su <- s(df$upper)
  # Define the modelling distribution and polynomial approximation
  Fz <- function(z, b) plogis(g(z, b))
  g <- function(z, b) b[1] + b[2]*z
  logL <- function(b, l, u) -sum(log(Fz(u, b) - Fz(l, b)))
  # Estimate the parameters
  mle <- nlm(logL, c(1, 1), u=su, l=sl)
  stopifnot("The optimizer failed." = mle$minimum > 1)
  # Save the fitted function
  f <- function(z, x) plogis(g(s(z), mle$estimate))
  attr(f, "AIC") <- 2*mle$minimum + 2*length(mle$estimate)
  return(f)
}


get_estimates <- function(df, D) {
  KM.mid <- survfit(Surv(mid, status == "Dead") ~ 1, data = df)
  KM.obs <- survfit(Surv(obs, status == "Dead") ~ 1, data = df)
  
  WEIB.mid <- survreg(Surv(mid, status == "Dead") ~ 1, data = df)
  WEIB.obs <- survreg(Surv(obs, status == "Dead") ~ 1, data = df)
  WEIB.int <- survreg(Surv(lower, upper, type = "interval2") ~ 1, data = df)
  
  PA.GAMMA <- proposed_approach(df, "gamma")
  PA.W <- proposed_approach(df, "weibull")
  PA.GEN.GAMMA <- generalized_approach(df, PA.GAMMA)
  PA.GEN.W <- generalized_approach(df, PA.W)
  
  
  pD <- array(c(
    summary(KM.mid, times = D)[[6]],
    summary(KM.obs, times = D)[[6]],
    1 - pweibull(D, 1 / WEIB.mid$scale, exp(WEIB.mid$coefficients)),
    1 - pweibull(D, 1 / WEIB.obs$scale, exp(WEIB.obs$coefficients)),
    1 - pweibull(D, 1 / WEIB.int$scale, exp(WEIB.int$coefficients))
  ),
  dim = c(length(D), 5), 
  dimnames = list(D, c("KM.MID", "KM.IGN", "W.MID", "W.IGN", "W.IC"))
  )
  
  pars <- array(c(1 / WEIB.mid$scale, exp(WEIB.mid$coefficients), 
                  1 / WEIB.obs$scale, exp(WEIB.obs$coefficients), 
                  1 / WEIB.int$scale, exp(WEIB.int$coefficients)), 
                dim = c(2, 3),
                dimnames = list(c("shape", "scale"), c("MID", "IGN", "IC")))

  PA.pD <- array(c(summary(KM.mid, times = D)[[6]], 1 - c(PA.GAMMA(D), PA.W(D), PA.GEN.GAMMA(D), PA.GEN.W(D))),
                 dim = c(length(D), 5),
                 dimnames = list(D, c("KM.MID", "PA.GAMMA", "PA.W", "GEN.GAMMA", "GEN.W")))
  
  return(list(pD, pars, PA.pD))
}

## SUMMARIES -------------------------------------------------------------------

repeat_estimation <- function(iter, days, df.pars) {
  est.prob <- array(
    NA,
    dim = c(iter, length(days), 5, 3),
    dimnames = list(
      1:iter,
      days,
      c("KM.MID", "KM.IGN", "W.MID", "W.IGN", "W.IC"),
      1:3
    )
  )
  est.pars <- array(
    NA,
    dim = c(iter, 2, 3, 3),
    dimnames = list(
      1:iter, 
      c("shape", "scale"), 
      c("W.MID", "W.IGN", "W.IC"), 
      1:3
    )
  )
  
  est.prob.PA <- array(
    NA,
    dim = c(iter, length(days), 5, 3),
    dimnames = list(
      1:iter,
      days,
      c("KM.MID", "PA.GAMMA", "PA.W", "GEN.GAMMA", "GEN.W"),
      1:3
    )
  )
  
  for(i in 1:iter) {
    df <- do.call(generate_df, df.pars)
    for(j in 1:3) {
      est <- get_estimates(df[df$register == as.character(j),], days)
      est.prob[i,,,j] <- est[[1]]
      est.pars[i,,,j] <- est[[2]]
      est.prob.PA[i,,,j] <- est[[3]]
    }
  }
  return(list("prob" = est.prob, "pars" = est.pars, "PA" = est.prob.PA))
}


average_bias <- function(thetas, true.vals) {
  bias <- apply(thetas, c(1, 3, 4), function(v) v - true.vals)
  avg.bias <- apply(bias, c(1, 3, 4), mean)
  return(avg.bias)
}

## PLOTTING --------------------------------------------------------------------

plot_bias <- function(bias, cols) {
  labs <- dimnames(bias)
  days <- as.numeric(labs[[1]])
  ymax <- max(max(bias), max(-bias)) + 0.0005
  plot(days, bias[,1], ylim = c(-ymax, ymax), pch = 20, col = cols[1],
       xaxt = "n", xlab = "Days", ylab = "Bias")
  axis(1, at = days, labels = TRUE)
  points(days, bias[,2], pch = 8, col = cols[2])
  points(days, bias[,3], pch = 21, col = cols[3])
  points(days, bias[,4], pch = 24, col = cols[4])
  points(days, bias[,5], pch = 25, col = cols[5])
  lines(c(0, 2000), c(0, 0), col = "grey", lwd = 0.5)
  legend("bottom", pch = c(20, 8, 21, 24, 25), bty = "n", col = cols,
         legend = labs[[2]], ncol = 2)
}

## PARAMETERS ------------------------------------------------------------------
# Setting for the first simulation: 5-year survival
pars1 <- list("A" = c(0.6, 2000))
followup1 <- 1826

# Setting for the second simulation: 1-year survival
pars2 <- list("A" = c(0.7, 200))
followup2 <- 365

# Common parameters for the simulations
n.iter <- 1000
days1 <- c(365, 730, 1095, 1461, 1826)
days2 <- c(30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 333, 365)
P1 <- 1 - pweibull(days1, pars1[["A"]][1], pars1[["A"]][2])
P2 <- 1 - pweibull(days2, pars2[["A"]][1], pars2[["A"]][2])
# Registers provide interval-censored data with different interval widths.
widths <- c("1" = 1, "2" = 30, "3" = 365)
lost.prop <- c("1" = 0, "2" = 0, "3" = 0)  # no loss to follow-up
treat.prop <- list("1" = list("A" = 1), "2" = list("A" = 1), "3" = list("A" = 1))
duration <- 365  # enrollment duration

m.cols <- c("darkorange", "red", "forestgreen", "royalblue", "cyan3")

## SIMULATION, N=500 -----------------------------------------------------------
N <- 2000 
sizes <- c("1" = N, "2" = N, "3" = N)

# first simulation
df1.pars <- list(sizes, treat.prop, pars1, duration, followup1, lost.prop, widths)

est1 <- repeat_estimation(n.iter, days1, df1.pars)
par1.bias <- average_bias(est1[["pars"]], pars1[["A"]])
prob1.bias <- average_bias(est1[["prob"]], P1)
PA1.bias <- average_bias(est1[["PA"]], P1)

plot_bias(prob1.bias[,,1], m.cols)
plot_bias(prob1.bias[,,2], m.cols)
plot_bias(prob1.bias[,,3], m.cols)

# second simulation
df2.pars <- list(sizes, treat.prop, pars2, duration, followup2, lost.prop, widths)

est2 <- repeat_estimation(n.iter, days2, df2.pars)
par2.bias <- average_bias(est2[["pars"]], pars2[["A"]])
prob2.bias <- average_bias(est2[["prob"]], P2)
PA2.bias <- average_bias(est2[["PA"]], P2)

plot_bias(prob2.bias[,,1], m.cols)
plot_bias(prob2.bias[,,2], m.cols)
plot_bias(prob2.bias[,,3], m.cols)
