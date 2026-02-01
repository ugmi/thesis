# This is a preliminary script to practice simulating survival data.

# We practice generating data and calculating the bias of different methods.
# We model varying measurement precision in this script.

library(survival)
#library(dplyr)
set.seed(45)

# We create helper functions to help reformat the dates.
month_start <- function(ymd) as.Date(gsub("[0-9]{2}$", "01", ymd))
month_end <- function(ymd) month_start(ymd + 32) - 1

# Generating the data ----------------------------------------------------------
# See script `data_generation_v3.R` for more details.
N <- 1000

# Generate survival times
t <- rweibull(N, shape = 0.8, scale = 500)
o <- sort(t)

# Create the data set
# Assume the study started on 2001/01/01, which corresponds to 11323
df <- data.frame("true" = round(t, 5), 
                 "days" = floor(t),
                 "lost.flw" = rbinom(N, 1, 1/14),
                 "days.from.start" = floor(runif(N, 1, 1095)))
df["entry"] <- as.Date("2001-01-01") + df$days.from.start
df["status"] <- ifelse(df$lost.flw, "Lost to follow up",
                       ifelse(df$days.from.start + df$days > 1095, "Alive", "Dead"))
df["exit"] <- ifelse(df$status == "Alive", 1095,
                     df$days.from.start + ifelse(!df$lost.flw, df$days, 
                                                 floor(pmin(1095 - df$days.from.start, df$days)*runif(N))))
df["exit"] <- as.Date("2001-01-01") + df$exit
df["country"] <- factor(rbinom(N, 1, 1/3), levels = 0:1, 
                        labels = c("Sweden", "Norway"))
df["exit.lower"] <- as.Date(ifelse(df$country == "Sweden" | df$status != "Dead", 
                                   df$exit, 
                                   pmax(month_start(df$exit), df$entry)))
df["exit.upper"] <- as.Date(ifelse(df$country == "Sweden" | df$status != "Dead", 
                                   df$exit, 
                                   month_end(df$exit)))

# Select only observed variables
df.obs <- df[c("status", "country", "entry", "exit.lower", "exit.upper")]

# Modelling using proposed approach --------------------------------------------

# Helper functions
fy <- function(l, u, b) Fz(u, b) - Fz(l, b)
logL <- function(b, l, u) -sum(log(fy(l, u, b)))

# Create coarsening intervals
l <- as.numeric(df.obs$exit.lower - df.obs$entry)
u <- ifelse(df.obs$status == "Dead", df.obs$exit.upper - df.obs$entry + 1, Inf)

# Try different distributions
Fz <- function(z, b) pexp(z, rate = b)
mle.exp <- nlm(logL, 1/100, l, u, hessian = TRUE)

Fz <- function(z, b) pweibull(z, shape = b[1], scale = b[2])
mle.weibull <- nlm(logL, c(0.5, 400), l, u, hessian = TRUE)

Fz <- function(z, b) pgamma(z, shape = b[1], rate = b[2])
mle.gamma <- nlm(logL, c(1, 1/100), l, u, hessian = TRUE)

# Plot and compare the distributions
plot(o, 1 - pweibull(o, 0.8, 500), pch=20, main = "Fitted survival curves",
     xlab = "Days", ylab = "Survival probability", cex=0.1)
points(o, 1 - pexp(o, mle.exp$estimate), pch = 20, col = "#0000FF66", cex=0.1)
points(o, 1 - pweibull(o, mle.weibull$estimate[1], mle.weibull$estimate[2]), 
       pch = 20, col = "#72BF4066", cex=0.1)
points(o, 1 - pgamma(o, mle.gamma$estimate[1], mle.gamma$estimate[2]),
       pch = 20, col = "#FF000066", cex=0.1)
legend("topright", lty = rep(1, 4), bty = "n",
       col = c("black", "#0000FF99", "#72BF4099", "#FF000099"), 
       legend = c("True distribution", "Exponential", "Weibull", "Gamma"))

# Bias (model parameters)
c(0.8, 500) - mle.weibull$estimate  #  -0.05825332 -25.41698875

# Bias (1-year survival probability)
(1 - pweibull(365.25, 0.8, 500)) - 
  (1 - pweibull(365.25, mle.weibull$estimate[1], mle.weibull$estimate[2]))
(1 - pweibull(365.25, 0.8, 500)) - 
  (1 - pgamma(365.25, mle.gamma$estimate[1], mle.weibull$estimate[2]))
(1 - pweibull(365.25, 0.8, 500)) - 
  (1 - pexp(365.25, mle.exp$estimate))
# -0.02158715, 0.4593924, -0.02762035

# Midpoint imputation ----------------------------------------------------------
df.obs["exit.mid"] <- as.Date(ifelse(df.obs$exit.upper == df.obs$exit.lower,
                                     df.obs$exit.lower,
                                     df.obs$exit.lower + 
                                       floor((df.obs$exit.upper - df.obs$exit.lower)/2)))
df.obs["days.mid"] <- as.numeric(df.obs$exit.mid - df.obs$entry)

# Fit the Kaplan-Meier survival curve
KM <- survfit(Surv(days.mid, status == "Dead") ~ 1, data = df.obs, 
              conf.type = "log-log")
v <- KM$time  # unique time points for which we can estimate the probabilities

# Compare the Kaplan-Meier curve with the true distribution
plot(KM, main = "Kaplan-Meier survival curve", xlab = "Days", ylab = "Survival")
points(v, 1 - pweibull(v, 0.8, 500), pch = 20, cex = 0.1, col = "red")

# Compare with approaches above
plot(KM, main = "Kaplan-Meier vs proposed approach", xlab = "Days", ylab = "Survival")
points(v, 1 - pexp(v, mle.exp$estimate), pch = 20, col = "#0000FF66", cex = 0.1)
points(v, 1 - pweibull(v, mle.weibull$estimate[1], mle.weibull$estimate[2]), 
       pch = 20, col = "#72BF4066", cex = 0.1)
points(v, 1 - pgamma(v, mle.gamma$estimate[1], mle.gamma$estimate[2]),
       pch = 20, col = "#FF000066", cex = 0.1)
legend("topright", lty = rep(1, 4), bty = "n",
       col = c("black", "#0000FF99", "#72BF4099", "#FF000099"), 
       legend = c("Kaplan-Meier curve", "Exponential", "Weibull", "Gamma"))

# Mean squared error
mean((1 - pweibull(v, 0.8, 500) - KM$surv)^2)
mean((pweibull(v, 0.8, 500) - pweibull(v, mle.weibull$estimate[1], mle.weibull$estimate[2]))^2)
mean((pweibull(v, 0.8, 500) - pgamma(v, mle.gamma$estimate[1], mle.gamma$estimate[2]))^2)
mean((pweibull(v, 0.8, 500) - pexp(v, mle.exp$estimate))^2)
# 0.000501347, 0.0004116628, 0.0004963145, 0.001627313


