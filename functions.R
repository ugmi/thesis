# Helper functions for the full simulation.

library(survival)
library(icenReg)
library(flexsurv)
library(msm)

## DATA GENERATION -------------------------------------------------------------
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
#' @param exact Binary indicator whether the survival time is exact
#' @param infinite Binary indicator if the upper limit is infinite
#' @return A data frame with two columns: lower and upper censoring limit
generate_intervals <- function(exit, exact, infinite) {
  L <- ifelse(
    exact | infinite, 
    exit, 
    exit - floor(runif(length(exact), 0, pmin(exit, 30)))
  )
  U <- ifelse(exact, exit, ifelse(infinite, Inf, L + 30))
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
generate_df <- function(n, proportions, pars, duration, followup, prob) {
  df <- generate_weibull(n, proportions, pars)
  df.censoring <- generate_censoring(
    ceiling(df$time), duration, followup, prob[df$register])
  df.intervals <- generate_intervals(
    df.censoring$exit,
    exact = df$register == "1" & df.censoring$status == "Dead",
    infinite = df.censoring$status != "Dead"
  )
  df.observed <- data.frame("obs" = ifelse(
    df$register == 1, 
    df.censoring$exit, 
    ifelse(df.censoring$status == "Dead", df.intervals$upper, df.censoring$exit)
    ))
  df.imputed <- data.frame("mid" = ifelse(
    df$register == "2" & df.censoring$status == "Dead", 
    df.observed$obs - 15, 
    df.observed$obs)
  )
  return(do.call(
    cbind, list(df, df.censoring, df.intervals, df.observed, df.imputed)))
}


## PROPOSED APPROACH -----------------------------------------------------------
#' Fit a distribution using the proposed approach
#' 
#' @param df A data frame with study data
#' @param dist The name of the parametric distribution to fit
#' @return A fitted cumulative distribution function of the modelling variable
proposed_approach <- function(df, dist) {
  fy <- function(l, u, b, x) Fz(u, b, x) - Fz(l, b, x)
  logL <- function(b, l, u, x) -sum(log(fy(l, u, b, x)))
  
  # Create coarsening intervals
  l <- ifelse(df$lower - df$upper == 0, df$lower - 1, df$lower)
  u <- df$upper
  # Covariate vector
  x <- df$treatment
  
  if (dist %in% c("exp", "exponential")) {
    Fz <- function(z, b, x) {
      (x == "A")*pexp(z, b[1]) + (x == "B")*pexp(z, b[2])
    }
    mle <- nlm(logL, c(1/100, 1/100), l, u, x, hessian = TRUE)
  } else if (dist == "gamma") {
    Fz <- function(z, b, x) {
      (x == "A")*pgamma(z, b[1], b[2]) + (x == "B")*pgamma(z, b[3], b[4])
    }
    mle <- nlm(logL, c(1, 1/100, 1, 1/100), l, u, x, hessian = TRUE)
  } else if (dist == "weibull") {
    Fz <- function(z, b, x) {
      (x == "A")*pweibull(z, exp(b[1]), exp(b[2])) + 
        (x == "B")*pweibull(z, exp(b[3]), exp(b[4]))
    }
    mle <- nlm(logL, c(-0.2, 6, -0.2, 6), l, u, x, hessian = TRUE)
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
  s <- function(z) qlogis(pweibull(z, 0.5, 2000))
  # Create the coarsening limits
  sl <- s(ifelse(df$lower - df$upper == 0, df$lower - 1, df$lower))
  su <- s(df$upper)
  # Define the modelling distribution and polynomial approximation
  Fz <- function(z, b, x) plogis(g(z, b, x))
  g <- function(z, b, x) {
    (b[1] + b[2]*(x == "B")) + z*(b[3] + b[4]*(x == "B"))
  }
  fy <- function(l, u, b, x) Fz(u, b, x) - Fz(l, b, x)
  logL <- function(b, l, u, x) -sum(log(fy(l, u, b, x)))
  # Estimate the parameters
  mle <- nlm(logL, rep(1, 4), u=su, l=sl, x=df$treatment)
  stopifnot("The optimizer failed." = mle$minimum > 1)
  # Save the fitted function
  f <- function(z, x) plogis(g(s(z), mle$estimate, x))
  attr(f, "AIC") <- 2*mle$minimum + 2*length(mle$estimate)
  return(f)
}


wilsonCI <- function(p.hat, n, alpha=0.05) {
  z <- qnorm(1 - alpha/2)
  lower <- (p.hat + z^2/(2*n))/(1 + z^2/n) - 
    z*sqrt(p.hat*(1 - p.hat)/n + z^2/(4*n^2))/(1 + z^2/n)
  upper <- (p.hat + z^2/(2*n))/(1 + z^2/n)+
    z*sqrt(p.hat*(1 - p.hat)/n + z^2/(4*n^2))/(1 + z^2/n)
  return(c(lower, upper))
}


## MODEL FITTING ---------------------------------------------------------------
fit_models <- function(df, ph = FALSE) {
  
  ## MODELS ASSUMING DAILY DATA
  # 1. Kaplan-Meier estimator
  fit1 <- survfit(Surv(exit, status == "Dead") ~ treatment,
                  data = df,
                  conf.type = "log-log")
  
  # 2. Cox proportional hazards model
  fit2 <- coxph(Surv(exit, status == "Dead") ~ treatment, data = df)
  
  # 3. Weibull AFT/PH model
  fit3 <- survreg(Surv(exit, status == "Dead") ~ treatment,
                  data = df,
                  dist = "weibull")
  
  # 4. Weibull AFT/PH stratified model
  fit4 <- survreg(Surv(exit, status == "Dead") ~ strata(treatment),
                  data = df,
                  dist = "weibull")
  
  
  ## MDOELS IGNORING CENSORING
  # 5. Kaplan-Meier estimator
  fit5 <- survfit(Surv(obs, status == "Dead") ~ treatment,
                  data = df,
                  conf.type = "log-log")
  
  # 6. Cox proportional hazards model
  fit6 <- coxph(Surv(obs, status == "Dead") ~ treatment, data = df)
  
  # 7. Weibull AFT/PH model
  fit7 <- survreg(Surv(obs, status == "Dead") ~ treatment,
                  data = df,
                  dist = "weibull")
  
  # 8. Weibull AFT/PH stratified model
  fit8 <- survreg(Surv(obs, status == "Dead") ~ strata(treatment),
                  data = df,
                  dist = "weibull")
  
  
  ## MODELS FOR DATA WITH MIDPOINT IMPUTATION
  # 9. Kaplan-Meier estimator
  fit9 <- survfit(Surv(mid, status == "Dead") ~ treatment,
                  data = df,
                  conf.type = "log-log")
  
  # 10. Cox proportional hazards model
  fit10 <- coxph(Surv(mid, status == "Dead") ~ treatment, data = df)
  
  # 11. Weibull AFT/PH model
  fit11 <- survreg(Surv(mid, status == "Dead") ~ treatment,
                   data = df,
                   dist = "weibull")
  
  # 12. Weibull AFT/PH stratified model
  fit12 <- survreg(Surv(mid, status == "Dead") ~ strata(treatment),
                   data = df,
                   dist = "weibull")
  
  
  ## MODELS FOR INTERVAL-CENSORED DATA
  # 13. Non-Parametric Maximum Likelihood Estimator (NPMLE) aka Turnbull’s Estimator
  fit13 <- ic_np(Surv(lower, upper, type = "interval2") ~ treatment, df)
  
  # 14. Semi-Parametric proportional hazards model
  # only fit if the hazards are proportional since otherwise won't converge
  if(ph) {
    fit14 <- ic_sp(Surv(lower, upper, type = "interval2") ~ treatment,
                   data = df,
                   model = "ph")
  } else {
    fit14 <- NA
  }
  
  # 15. Semi-Parametric proportional odds model
  fit15 <- ic_sp(Surv(lower, upper, type = "interval2") ~ treatment,
                 data = df,
                 model = "po")
  
  # 16. Weibull AFT/PH
  fit16 <- survreg(Surv(lower, upper, type = "interval2") ~ treatment,
                   data = df,
                   dist = "weibull")
  
  # 17. Weibull AFT/PH stratified
  fit17 <- survreg(
    Surv(lower, upper, type = "interval2") ~ strata(treatment),
    data = df,
    dist = "weibull"
  )
  
  
  ## FLEXIBLE PARAMETRIC MODELS FOR INTERVAL-CENSORED DATA
  # 18-20. Splines
  suppressWarnings(fit18 <- flexsurvspline(
    Surv(lower, upper, type = "interval2") ~ treatment, 
    data = df, 
    scale = "hazard", 
    k = 3, 
    anc = list(gamma1 = ~ treatment)
  ))
  suppressWarnings(fit19 <- flexsurvspline(
    Surv(lower, upper, type = "interval2") ~ treatment, 
    data = df, 
    scale = "hazard", 
    k = 3
  ))
  suppressWarnings(fit20 <- flexsurvspline(
    Surv(lower, upper, type = "interval2") ~ treatment, 
    data = df, 
    scale = "hazard", 
    k = 2, 
    anc = list(gamma1 = ~ treatment)
  ))
  
  
  # 21-22. Parametric
  fit21 <- flexsurvreg(Surv(lower, upper, type = "interval2") ~ treatment, 
                       data = df, 
                       dist = "weibull")
  fit22 <- flexsurvreg(Surv(lower, upper, type = "interval2") ~ treatment, 
                       data = df, 
                       dist = "weibull", 
                       anc = list(shape = ~ treatment))
  
  
  ## PROPOSED APPROACH
  suppressWarnings(fit23 <- proposed_approach(df, "weibull"))
  suppressWarnings(fit24 <- proposed_approach(df, "gamma"))
  suppressWarnings(fit25 <- generalized_approach(df))
  
  fit <- mget(paste0("fit", 1:25))
  names(fit) <- c("KM.DAILY", "COX.DAILY", "W.DAILY", "W.STR.DAILY",
                  "KM.IGN", "COX.IGN", "W.IGN", "W.STR.IGN",
                  "KM.IMP", "COX.IMP", "W.IMP", "W.STR.IMP",
                  "NPMLE", "SEMI.PH", "SEMI.PO", "W.IC", "W.STR.IC",
                  "SPL31", "SPL30", "SPL21", "W.FPM", "W.FPM1",
                  "PA.W", "PA.GAMMA", "PA.GEN")
  
  return(fit)
}


get_SE_weibull <- function(model, strata = FALSE) {
  if(strata) {
    SE.A <- deltamethod(~ exp(-(1826 / exp(x2))^(exp(-x1))), 
                        c(log(model$scale[[1]]), model$coefficients[[1]]), 
                        vcov(model)[1:2, 1:2])
    SE.B <- deltamethod(~ exp(-(1826 / exp(x2))^(exp(-x1))), 
                        c(log(model$scale[[2]]), model$coefficients[[1]]), 
                        vcov(model)[c(1, 3), c(1, 3)])
  } else {
    SE.A <- deltamethod(~ exp(-(1826 / exp(x1))^(exp(-x2))), 
                        c(model$coefficients[[1]], log(model$scale)), 
                        vcov(model)[c(1, 3), c(1, 3)])
    SE.B <- deltamethod(~ exp(-(1826 / exp(x1 + x2))^(exp(-x3))), 
                        c(model$coefficients, log(model$scale)), 
                        vcov(model))
  }
  return(c(SE.A, SE.B))
}


print_surv_table <- function(lst) {
  cols <- attr(lst, "estimates")
  rows <- attr(lst, "models")
  array(lst, dim = c(length(cols), length(rows)), dimnames = list(cols, rows))
}


get_surv_basic <- function(D, models) {
  KM <- summary(models[[1]], times = D)
  COX <- predict(
    models[[2]],
    data.frame(exit = D, obs = D, mid = D, status = "Dead", 
               treatment = c("A", "B")),
    type = "survival",
    se.fit = TRUE
  )
  COX.lower <- COX$fit - qnorm(0.975) * COX$se.fit
  COX.upper <- COX$fit + qnorm(0.975) * COX$se.fit
  COX.AIC <- AIC(models[[2]])
  AFT <- 1 - pweibull(D, 
                      shape = 1 / models[[3]]$scale,
                      scale = c(exp(models[[3]]$coefficients[[1]]),
                                exp(sum(models[[3]]$coefficients))))
  AFT.SE <- get_SE_weibull(models[[3]])
  AFT.CI.A <- AFT[1] + c(- qnorm(0.975), qnorm(0.975)) * AFT.SE[1]
  AFT.CI.B <- AFT[2] + c(- qnorm(0.975), qnorm(0.975)) * AFT.SE[2]
  AFT.AIC <- AIC(models[[3]])
  AFT.STRAT <- 1 - pweibull(D,
                            shape = 1 / models[[4]]$scale,
                            scale = exp(models[[4]]$coefficients))
  # <==> exp(-(D/exp(models[[4]]$coefficients))^(1/models[[4]]$scale))
  AFT.STRAT.SE <- get_SE_weibull(models[[4]], strata = TRUE)
  AFT.STRAT.CI.A <- AFT.STRAT[1] + c(- qnorm(0.975), qnorm(0.975)) * AFT.STRAT.SE[1]
  AFT.STRAT.CI.B <- AFT.STRAT[2] + c(- qnorm(0.975), qnorm(0.975)) * AFT.STRAT.SE[2]
  AFT.STRAT.AIC <- AIC(models[[4]])
  
  lst <- c(KM$surv[1], KM$lower[1], KM$upper[1], KM$std.err[1], 
           KM$surv[2], KM$lower[2], KM$upper[2], KM$std.err[2], NA,
           COX$fit[1], COX.lower[1], COX.upper[1], COX$se.fit[1], 
           COX$fit[2], COX.lower[2], COX.upper[2], COX$se.fit[2], COX.AIC,
           AFT[1], AFT.CI.A, AFT.SE[1], AFT[2], AFT.CI.B, AFT.SE[2], AFT.AIC,
           AFT.STRAT[[1]], AFT.STRAT.CI.A, AFT.STRAT.SE[1],
           AFT.STRAT[[2]], AFT.STRAT.CI.B, AFT.STRAT.SE[2], AFT.STRAT.AIC)
  
  attr(lst, "models") <- c("KM", "COX", "AFT", "AFT.STRAT")
  attr(lst, "estimates") <- c(
    "A", "A.lower", "A.upper", "A.SE", "B", "B.lower", "B.upper", "B.SE", "AIC")
  return(lst)
}


get_surv_interval <- function(D, models) {
  survA <- models[[1]]$scurves$A
  survB <- models[[1]]$scurves$B
  NPMLE <- c(survA$S_curves[[1]][survA$Tbull_ints[, 2] >= D][1],
             survB$S_curves[[1]][survB$Tbull_ints[, 2] >= D][1])
  if(class(models[[2]]) == "ic_ph") {
    SEMI.PH <- 1 - getFitEsts(models[[2]], 
                              newdata = data.frame(treatment = c("A", "B")), 
                              q = 1826)
  } else {
    SEMI.PH <- c(NA, NA)
  }
  if(class(models[[3]]) == "ic_po"){
    SEMI.PO <- 1 - getFitEsts(models[[3]], 
                              newdata = data.frame(treatment = c("A", "B")), 
                              q = 1826)
  } else {
    SEMI.PO <- c(NA, NA)
  }
  AFT <- 1 - pweibull(D, 
                      shape = 1 / models[[4]]$scale,
                      scale = c(exp(models[[4]]$coefficients[[1]]),
                                exp(sum(models[[4]]$coefficients))))
  AFT.SE <- get_SE_weibull(models[[4]])
  AFT.CI.A <- AFT[1] + c(- qnorm(0.975), qnorm(0.975)) * AFT.SE[1]
  AFT.CI.B <- AFT[2] + c(- qnorm(0.975), qnorm(0.975)) * AFT.SE[2]
  AFT.AIC <- AIC(models[[4]])
  AFT.STRAT <- 1 - pweibull(D,
                            shape = 1 / models[[5]]$scale,
                            scale = exp(models[[5]]$coefficients))
  AFT.STRAT.SE <- get_SE_weibull(models[[5]], strata = TRUE)
  AFT.STRAT.CI.A <- AFT.STRAT[1] + c(- qnorm(0.975), qnorm(0.975)) * AFT.STRAT.SE[1]
  AFT.STRAT.CI.B <- AFT.STRAT[2] + c(- qnorm(0.975), qnorm(0.975)) * AFT.STRAT.SE[2]
  AFT.STRAT.AIC <- AIC(models[[5]])
  
  lst <- c(NPMLE[1], NA, NA, NA, NPMLE[2], NA, NA, NA, NA,
           SEMI.PH[1], NA, NA, NA, SEMI.PH[2], NA, NA, NA, NA,
           SEMI.PO[1], NA, NA, NA, SEMI.PO[2], NA, NA, NA, NA,
           AFT[1], AFT.CI.A, AFT.SE[1], AFT[2], AFT.CI.B, AFT.SE[2], AFT.AIC,
           AFT.STRAT[[1]], AFT.STRAT.CI.A, AFT.STRAT.SE[1],
           AFT.STRAT[[2]], AFT.STRAT.CI.B, AFT.STRAT.SE[2], AFT.STRAT.AIC)
  
  attr(lst, "models") <- c("KM", "COX", "AFT", "AFT.STRAT")
  attr(lst, "estimates") <- c(
    "A", "A.lower", "A.upper", "A.SE", "B", "B.lower", "B.upper", "B.SE", "AIC")
  return(lst)
}


get_row_PA <- function(pdist, sizes, D) {
  p.A <- pdist(D, "A")
  p.B <- pdist(D, "B")
  CI.A <- wilsonCI(1-p.A, sizes[["A"]])
  CI.B <- wilsonCI(1-p.B, sizes[["B"]])
  AIC <- attr(pdist, "AIC")
  return(c(1-p.A, CI.A, NA, 1-p.B, CI.B, NA, AIC))
}


get_row_fpm <- function(model, D) {
  info <- summary(model, type = "survival", ci = TRUE, t = D, se = TRUE)
  return(c(info[[1]][2:5], info[[2]][2:5], model$AIC))
}


get_surv_flexible <- function(D, models) {
  lst <- unlist(sapply(models, get_row_fpm, D = D))
  attr(lst, "models") <- names(models)
  attr(lst, "estimates") <- c(
    "A", "A.lower", "A.upper", "A.SE", "B", "B.lower", "B.upper", "B.SE", "AIC")
  return(lst) 
}


get_table <- function(D, df) {
  N.tr <- table(df$treatment)
  models <- fit_models(df)
  surv.day <- get_surv_basic(D, models[1:4])
  surv.upp <- get_surv_basic(D, models[5:8])
  surv.mid <- get_surv_basic(D, models[9:12])
  surv.int <- get_surv_interval(D, models[13:17])
  surv.fpm <- get_surv_flexible(D, models[18:22])
  surv.pa <- c(sapply(models[23:25], get_row_PA, D = D, sizes = N.tr))
  lst <- c(surv.day, surv.upp, surv.mid, surv.int, surv.fpm, surv.pa)
  tbl <- array(
    lst, 
    dim = c(9, 25), 
    dimnames = list(attr(surv.day, "estimates"), names(models))
  )
  return(tbl)
}

## PLOTS -----------------------------------------------------------------------
plot_bias <- function(tbl, cols, title.append = "") {
  y.max <- round(max(max(-tbl, na.rm = TRUE), max(tbl, na.rm = TRUE)) + 0.0005, 3)
  n.mod <- dim(tbl)[2]
  n.sizes <- dim(tbl)[1]
  col.grey <- "#99999999"
  grey.x <- 0:(n.mod + 1)
  reps.y <- n.mod + 2
  par(mar = c(6, 5, 4, 6.5))
  plot(
    tbl[1,],
    ylim = c(-y.max, y.max),
    col = cols[1],
    pch = 20,
    cex = 0.5,
    xaxt = "n",
    xlab = "",
    las = 1,
    ylab = "Bias",
    main = paste0("Bias of 5-year survival prob.", title.append)
  )
  axis(
    1,
    at = 1:n.mod,
    labels = colnames(tbl),
    las = 2,
    cex.axis = 0.7
  )
  for (i in 2:n.sizes) 
    points(tbl[i, ], col = cols[i], pch = 20, cex = 0.5)
  lines(grey.x, rep(0.01, reps.y), lty = 2, col = col.grey, lwd = 0.5)
  lines(grey.x, rep(-0.01, reps.y), lty = 2, col = col.grey, lwd = 0.5)
  lines(grey.x, rep(0, reps.y), lty = 1, col = col.grey)
  legend(
    "right",
    inset = c(-0.25, 0),
    legend = paste0("N=", rownames(tbl)),
    col = cols,
    pch = rep(20, n.sizes),
    xpd = TRUE,
    bty = "n"
  )
}


plot_AIC <- function(tbl, subset, cols, ylab) {
  y.max <- round(max(tbl[subset, ]) + 50, -2)
  n.sub <- length(subset)
  n.sizes <- dim(tbl)[2]
  par(mar = c(6, 5, 4, 6.5))
  plot(
    tbl[subset, 1],
    col = cols[1],
    pch = 20,
    cex = 0.5,
    xaxt = "n",
    xlab = "",
    las = 1,
    ylab = ylab,
    main = "AIC by N",
    ylim = c(0, y.max)
  )
  axis(
    1,
    at = 1:n.sub,
    labels = rownames(tbl)[subset],
    las = 2,
    cex.axis = 0.7
  )
  for (i in 2:n.sizes)
    points(tbl[subset, i], col = cols[i], pch = 20, cex = 0.5)
  legend(
    "right",
    inset = c(-0.25, 0),
    legend = paste0("N=", colnames(tbl)),
    col = cols,
    pch = rep(20, n.sizes),
    xpd = TRUE,
    bty = "n"
  )
}


plot_se <- function(tbl, cols, title.append = "") {
  y.max <- round(
    max(max(-tbl, na.rm = TRUE), max(tbl, na.rm = TRUE)) + 0.0005, 3)
  n.mod <- dim(tbl)[2]
  n.sizes <- dim(tbl)[1]
  par(mar = c(6, 5, 4, 6.5))
  plot(
    tbl[1,],
    ylim = c(0, y.max),
    col = cols[1],
    pch = 20,
    cex = 0.5,
    xaxt = "n",
    xlab = "",
    las = 1,
    ylab = "SE",
    main = paste0("SE for 5-year survival prob.", title.append)
  )
  axis(
    1,
    at = 1:n.mod,
    labels = colnames(tbl),
    las = 2,
    cex.axis = 0.7
  )
  for (i in 2:n.sizes) 
    points(tbl[i, ], col = cols[i], pch = 20, cex = 0.5)
  lines(
    0:(n.mod + 1), 
    rep(0.005, n.mod + 2), 
    lty = 2, 
    col = "#99999999", 
    lwd = 0.5
  )
  legend(
    "right",
    inset = c(-0.25, 0),
    legend = paste0("N=", rownames(tbl)),
    col = cols,
    pch = rep(20, n.sizes),
    xpd = TRUE,
    bty = "n"
  )
}



