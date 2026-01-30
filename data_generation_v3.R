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
     main = "Weibull distribution with \nshape=0.8 and scale=500", pch=20)

# We assume the data generated above is the "true" distribution.
# However, in real data, the measurements are discrete, often rounded to the day.
# We also often have censoring because people move away or the follow-up ends.

# We model this by modifying the generated data vector.

# Suppose our study ran for three years.
# Then, the survival times of people alive after day 1095 of the study are censored.
# We assume 1/14 of participants were lost to follow up.
df <- data.frame("true" = round(t, 5), 
                 "days" = floor(t),
                 "lost.flw" = rbinom(N, 1, 1/14),
                 "t.entry" = floor(runif(N, 1, 1095)))
df["status"] <- ifelse(df$lost.flw, "Lost to follow up",
                       ifelse(df$t.entry + df$days > 1095, "Alive", "Dead"))
df["t.exit"] <- ifelse(df$status == "Alive", 1095,
                       df$t.entry + ifelse(!df$lost.flw, df$days, 
                                            floor(pmin(1095 - df$t.entry, df$days)*runif(N))))
df.obs <- df[c("status", "t.entry", "t.exit")]

# Modelling using proposed approach --------------------------------------------
# Define coarsening limits
l <- df.obs$t.exit - df.obs$t.entry
u <- ifelse(df.obs$status == "Dead", l + 1, Inf)

# Helper functions
fy <- function(l, u, b) Fz(u, b) - Fz(l, b)
logL <- function(b, l, u) -sum(log(fy(l, u, b)))

# Try different distributions
Fz <- function(z, b) pexp(z, b)
mle.exp <- nlm(logL, 1/100, l, u, hessian = TRUE)

Fz <- function(z, b) pweibull(z, b[1], b[2])
mle.weibull <- nlm(logL, c(0.5, 400), l, u, hessian = TRUE)

Fz <- function(z, b) plogis(z, b[1], b[2])
mle.logis <- nlm(logL, c(0, 400), l, u, hessian = TRUE)

# Plot and compare the distributions
plot(o, 1 - pweibull(o, 0.8, 500), pch=20, main = "Fitted survival curves",
     xlab = "Days", ylab = "Survival probability")
points(o, 1 - pexp(o, mle.exp$estimate), pch = 20, col = "#FF00A033")
points(o, 1 - pweibull(o, mle.weibull$estimate[1], mle.weibull$estimate[2]), 
       pch = 20, col = "#80DF3033")
points(o, 1 - plogis(o, mle.logis$estimate[1], mle.logis$estimate[2]), 
       pch = 20, col = "#3060FF33")
legend("topright", lty = rep(1, 4), bty = "n",
       col = c("black", "#FF00A099", "#80DF3099", "#3060FF99"), 
       legend = c("True distribution", "Exponential", "Weibull", "Logistic"))

# Plot the Weibull distribution only
plot(o, 1 - pweibull(o, mle.weibull$estimate[1], mle.weibull$estimate[2]),
     pch = 20, col = "red",  main = "Weibull survival curves",
     xlab = "Days", ylab = "Survival probability")
points(o, 1 - pweibull(o, 0.8, 500), pch = 20)
legend("topright", lty = c(1, 1), bty = "n", col = c("black", "red"),
       legend = c("True distribution", "Fitted distribution"))

# Bias
c(0.8, 500) - mle.weibull$estimate

# Conventional approaches ------------------------------------------------------





