# This is a preliminary script to practice simulating survival data.

# We practice generating data and calculating the bias of different methods.
# We do not model the varying measurement precision in this script.

library(survival)
#library(dplyr)
set.seed(45)

# Generating the data ----------------------------------------------------------
N <- 1000

# First, we generate survival data assuming Weibull distribution.
t <- rweibull(N, shape = 0.8, scale = 500)
summary(t/365.25)

# We see that under this distribution,
# half of the people survive 0.915 years (roughly 334 days).
# Minimum survival was 0 days, maximum - 14.2 years.

o <- sort(t)
plot(o, 1 - pweibull(o, 0.8, 500), xlab = "Days", ylab = "Survival probability",
     main = "Weibull distribution with \nshape=0.8 and scale=500", pch = 20)

# We assume the data generated above is the "true" distribution.
# However, in real data, the measurements are discrete, often rounded to the day.
# We also often have censoring because people move away or the follow-up ends.

# We model this by modifying the generated data vector.

# Suppose our study ran for three years.
# Then, the survival times of people alive after day 1095 of the study are censored.
# We assume 1/14 of participants were lost to follow up and that it wasn't
# systematic, i.e. that it was independent from entry date and the survival time.
df <- data.frame("true" = round(t, 5), 
                 "days" = floor(t),
                 "lost.flw" = rbinom(N, 1, 1/14),
                 "entry" = floor(runif(N, 1, 1095)))
df["status"] <- ifelse(df$lost.flw, "Lost to follow up",
                       ifelse(df$entry + df$days > 1095, "Alive", "Dead"))
df["exit"] <- ifelse(df$status == "Alive", 1095,
                     df$entry + ifelse(!df$lost.flw, df$days, 
                                       floor(pmin(1095 - df$entry, df$days)*runif(N))))
df.obs <- df[c("status", "entry", "exit")]

# Modelling using proposed approach --------------------------------------------
# Define coarsening limits
l <- df.obs$exit - df.obs$entry
u <- ifelse(df.obs$status == "Dead", l + 1, Inf)

# Helper functions
fy <- function(l, u, b) Fz(u, b) - Fz(l, b)
logL <- function(b, l, u) -sum(log(fy(l, u, b)))

# Confidence interval formula
wilsonCI <- function(p.hat, n, alpha=0.05) {
  z <- qnorm(1 - alpha/2)
  lower <- (p.hat + z^2/(2*n))/(1 + z^2/n) - 
    z*sqrt(p.hat*(1 - p.hat)/n + z^2/(4*n^2))/(1 + z^2/n)
  upper <- (p.hat + z^2/(2*n))/(1 + z^2/n)+
    z*sqrt(p.hat*(1 - p.hat)/n + z^2/(4*n^2))/(1 + z^2/n)
  return(list("lower"=lower, "upper"=upper))
}

plotCI <- function(x, lower, upper, col, fill=NA) {
  if (!is.na(fill)) {
    polygon(c(rev(x), x), c(rev(lower), upper), col = fill, border = NA)
  }
  lines(x, lower, col = col, lty = "dashed")
  lines(x, upper, col = col, lty = "dashed")
}

# Try different distributions
Fz <- function(z, b) pexp(z, rate = b)
mle.exp <- nlm(logL, 1/100, l, u, hessian = TRUE)
CI.exp <- wilsonCI(1 - pexp(o, mle.exp$estimate), N)

Fz <- function(z, b) pweibull(z, shape = b[1], scale = b[2])
mle.weibull <- nlm(logL, c(0.5, 400), l, u, hessian = TRUE)
CI.weibull <- wilsonCI(1 - pweibull(o, mle.weibull$estimate[1], mle.weibull$estimate[2]), N)

Fz <- function(z, b) plogis(z, location = b[1], scale = b[2])
mle.logis <- nlm(logL, c(0, 400), l, u, hessian = TRUE)
CI.logis <- wilsonCI(1 - plogis(o, mle.logis$estimate[1], mle.logis$estimate[2]), N)

Fz <- function(z, b) pgamma(z, shape = b[1], rate = b[2])
mle.gamma <- nlm(logL, c(1, 1/100), l, u, hessian = TRUE)
CI.gamma <- wilsonCI(1 - pgamma(o, mle.gamma$estimate[1], mle.gamma$estimate[2]), N)

# Plot and compare the distributions
plot(o, 1 - pweibull(o, 0.8, 500), pch = 20, main = "Fitted survival curves",
     xlab = "Days", ylab = "Survival probability", cex = 0.05)
points(o, 1 - pexp(o, mle.exp$estimate), pch = 20, col = "#DF99FF66", cex = 0.1)
plotCI(o, CI.exp[[1]], CI.exp[[2]], "#DF99FF33")
points(o, 1 - pweibull(o, mle.weibull$estimate[1], mle.weibull$estimate[2]), 
       pch = 20, col = "#72BF4066", cex = 0.1)
plotCI(o, CI.weibull[[1]], CI.weibull[[2]], "#72BF4033")
points(o, 1 - plogis(o, mle.logis$estimate[1], mle.logis$estimate[2]), 
       pch = 20, col = "#3060FF66", cex = 0.1)
plotCI(o, CI.logis[[1]], CI.logis[[2]], "#3060FF33")
points(o, 1 - pgamma(o, mle.gamma$estimate[1], mle.gamma$estimate[2]),
       pch = 20, col = "#FF000066", cex = 0.1)
plotCI(o, CI.gamma[[1]], CI.gamma[[2]], "#FF000033")
legend("topright", lty = c(rep(1, 5), 2), bty = "n",
       col = c("black", "#DF99FF99", "#72BF4099", "#3060FF99", "#FF000099", "#00000033"), 
       legend = c("True distribution", "Exponential", "Weibull", "Logistic", "Gamma", "CI bounds"))

# We see in the plot that Weibull and Gamma distributions are close to the truth.

# Plot the Weibull distribution only
plot(o, 1 - pweibull(o, mle.weibull$estimate[1], mle.weibull$estimate[2]),
     pch = 20, col = "red",  cex=0.1, main = "Weibull survival curves",
     xlab = "Days", ylab = "Survival probability")
plotCI(o, CI.weibull[[1]], CI.weibull[[2]], "#FF000044", fill = "#FF000022")
points(o, 1 - pweibull(o, 0.8, 500), pch = 20, cex=0.1)
legend("topright", lty = c(1, 1, 2), bty = "n",
       col = c("black", "red", "#FF000044"),
       legend = c("True distribution", "Fitted distribution", "CI bounds"))

# Bias (model parameters)
c(0.8, 500) - mle.weibull$estimate  # -0.03495806 -39.29334380

# Bias (1-year survival probability)
(1 - pweibull(365.25, 0.8, 500)) - (1 - pweibull(365.25, mle.weibull$estimate[1], mle.weibull$estimate[2]))
(1 - pweibull(365.25, 0.8, 500)) - (1 - pgamma(365.25, mle.gamma$estimate[1], mle.weibull$estimate[2]))
(1 - pweibull(365.25, 0.8, 500)) - (1 - pexp(365.25, mle.exp$estimate))
# -0.02625908, 0.4593924, -0.03086711

# Conventional approaches ------------------------------------------------------
KM <- survfit(Surv(exit - entry, status == "Dead") ~ 1, data = df.obs, 
              conf.type = "log-log")
CI.size.KM <- KM$upper - KM$lower
v <- KM$time
plot(KM, main = "Kaplan-Meier survival curve", xlab = "Days", ylab = "Survival",
     col = "red", xlim = c(0, 5000))
points(o, 1 - pweibull(o, 0.8, 500), pch = 20, col = "black", cex = 0.1)
#lines(survfit(coxph(Surv(exit - entry, status == "Dead") ~ 1, data = df.obs)), col="blue")
# No need to add Cox model since we have no covariates

# Compare with some of the fitted models from before
plot(KM, main = "Kaplan-Meier vs proposed approach", xlab = "Days", ylab = "Survival")
points(o, 1 - pexp(o, mle.exp$estimate), pch = 20, col = "#DF99FF66", cex = 0.07)
plotCI(o, CI.exp[[1]], CI.exp[[2]], "#DF99FF33")
points(o, 1 - pweibull(o, mle.weibull$estimate[1], mle.weibull$estimate[2]), 
       pch = 20, col = "#72BF4066", cex = 0.07)
plotCI(o, CI.weibull[[1]], CI.weibull[[2]], "#72BF4033")
points(o, 1 - pgamma(o, mle.gamma$estimate[1], mle.gamma$estimate[2]),
       pch = 20, col = "#FF000066", cex = 0.07)
plotCI(o, CI.gamma[[1]], CI.gamma[[2]], "#FF000033")
legend("topright", lty = rep(1, 4), bty = "n",
       col = c("black", "#DF99FF99", "#72BF4099", "#FF000099", "#00000033"), 
       legend = c("Kaplan-Meier curve", "Exponential", "Weibull", "Gamma", "CI bounds"))

# Mean squared error
mean((1 - pweibull(KM$time, 0.8, 500) - KM$surv)^2)
mean((pweibull(KM$time, 0.8, 500) - pweibull(KM$time, mle.weibull$estimate[1], mle.weibull$estimate[2]))^2)
mean((pweibull(KM$time, 0.8, 500) - pgamma(KM$time, mle.gamma$estimate[1], mle.gamma$estimate[2]))^2)
mean((pweibull(KM$time, 0.8, 500) - pexp(KM$time, mle.exp$estimate))^2)

# Bias (1-year survival probability)
(1 - pweibull(365.25, 0.8, 500)) - KM$surv[290] # since KM$time[290] = 365
# -0.03037433

# Parametric models ------------------------------------------------------------
fit.weibull <- survreg(Surv(exit - entry, status == "Dead") ~ 1, 
                       data = df.obs[df.obs$exit - df.obs$entry != 0,], 
                       dist = "weibull")
plot(KM$time, 1 - pweibull(KM$time, shape = 1/fit.weibull$scale, scale = exp(fit.weibull$coefficients)),
     xlab = "Days", ylab = "Survival", pch = 20, cex = 0.1, col = "red")
points(v, 1 - pweibull(v, mle.weibull$estimate[1], mle.weibull$estimate[2]), 
       pch = 20, col = "#72BF4066", cex = 0.1)
points(v, 1 - pweibull(v, 0.8, 500), pch = 20, cex = 0.1)

# Bias (model parameters)
c(0.8, 500) - as.numeric(c(1/fit.weibull$scale, exp(fit.weibull$coefficients)))  #  -0.07256524 -39.81796054

# Likelihood formula



