# This is a preliminary script to practice simulating survival data.

# We practice generating data and calculating the bias of different methods.
# We do not model the varying measurement precision in this script but
# we include modelling the effect of two different treatments.

library(survival)
#library(dplyr)
set.seed(45)

# Generating the data ----------------------------------------------------------
# See script `data_generation_v3.R` for more details.
Na <- 600
Nb <- 400

# Generate survival times
ta <- rweibull(Na, shape = 0.8, scale = 500)
tb <- rweibull(Nb, shape = 0.7, scale = 600)

plot(ta, 1 - pweibull(ta, 0.8, 500), col = "red", pch = 20, cex = 0.1)
points(tb, 1 - pweibull(tb, 0.7, 600), col = "blue", pch = 20, cex = 0.1)

# Create the data set
# We assume there are two possible treatments, A and B.
df <- data.frame("true" = round(c(ta, tb), 5), 
                 "days" = floor(c(ta, tb)),
                 "lost.flw" = rbinom(Na + Nb, 1, 1/14),
                 "entry" = floor(runif(Na + Nb, 1, 1095)))
df["treatment"] <- as.factor(rep(c("A", "B"), c(Na, Nb)))
df["status"] <- ifelse(df$lost.flw, "Lost to follow up",
                       ifelse(df$entry + df$days > 1095, "Alive", "Dead"))
df["exit"] <- ifelse(df$status == "Alive", 1095,
                     df$entry + ifelse(!df$lost.flw, df$days, 
                                       floor(pmin(1095 - df$entry, df$days)*runif(Na + Nb))))
df.obs <- df[c("status", "treatment", "entry", "exit")]

# Modelling using proposed approach --------------------------------------------

# Helper functions
fy <- function(l, u, b, x) Fz(u, b, x) - Fz(l, b, x)
logL <- function(b, l, u, x) -sum(log(fy(l, u, b, x)))

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

# Define coarsening limits
l <- df.obs$exit - df.obs$entry
u <- ifelse(df.obs$status == "Dead", l + 1, Inf)

x <- df.obs$treatment

Fz <- function(z, b, x) {
  (x == "A")*pweibull(z, b[1], b[2]) + (x == "B")*pweibull(z, b[3], b[4])
}
mle.weibull <- nlm(logL, c(0.5, 400, 0.5, 400), l, u, x, hessian = TRUE)
oa <- sort(ta)
ob <- sort(tb)
CI.A <- wilsonCI(1 - Fz(oa, mle.weibull$estimate, "A"), Na)
CI.B <- wilsonCI(1 - Fz(ob, mle.weibull$estimate, "B"), Nb)

# Plot the general curve
o <- unique(c(oa, ob))

# Plot the result by treatment
plot(ta, 1 - Fz(ta, mle.weibull$estimate, "A"), pch = 20, cex = 0.1, 
     col = "#FF000099", main = "Estimated survival \n(proposed approach)", 
     xlab = "Days", ylab = "Survival")
plotCI(oa, CI.A[[1]], CI.A[[2]], "#FF000044")
points(tb, 1 - Fz(tb, mle.weibull$estimate, "B"), pch = 20, cex = 0.1, 
       col = "#0000FF99")
plotCI(ob, CI.B[[1]], CI.B[[2]], "#0000FF44")
points(ta, 1 - pweibull(ta, 0.8, 500), col = "#660000FF", pch = 20, cex = 0.1)
points(tb, 1 - pweibull(tb, 0.7, 600), col = "#000066FF", pch = 20, cex = 0.1)
legend("topright", lty = c(1, 1, 1, 1, 2), bty = "n", 
       col = c("#660000FF", "#FF000099", "#000066FF", "#0000FF99", "grey"),
       legend = c("True A", "Estimated A", "True B", "Estimated B", "CI bounds"))

# Conventional approaches ------------------------------------------------------
KM <- survfit(Surv(exit - entry, status == "Dead") ~ treatment, data = df.obs, 
              conf.type = "log-log")

plot(KM, col = c("#FF000099", "#0000FF99"), main = "Kaplan-Meier survival curves",
     xlab = "Days", ylab = "Survival")
plotCI(KM$time[1:389], KM$lower[1:389], KM$upper[1:389], "#FF000044")
plotCI(KM$time[390:685], KM$lower[390:685], KM$upper[390:685], "#0000FF44")
points(ta, 1 - pweibull(ta, 0.8, 500), col = "#990000FF", pch = 20, cex = 0.1)
points(tb, 1 - pweibull(tb, 0.7, 600), col = "#000099FF", pch = 20, cex = 0.1)
legend("topright", lty = c(1, 1, 1, 1, 2), bty = "n", 
       col = c("#660000FF", "#FF000099", "#000066FF", "#0000FF99", "grey"),
       legend = c("True A", "Estimated A", "True B", "Estimated B", "CI bounds"))

# Compare
plot(KM, col = c("#FF000099", "#0000FF99"), main = "Survival curves",
     xlab = "Days", ylab = "Survival", xlim=c(0,500))
plotCI(KM$time[1:389], KM$lower[1:389], KM$upper[1:389], "#FF000044")
plotCI(KM$time[390:685], KM$lower[390:685], KM$upper[390:685], "#0000FF44")
points(ta, 1 - Fz(ta, mle.weibull$estimate, "A"), pch = 20, cex = 0.1, 
       col = "#DD007788")
plotCI(oa, CI.A[[1]], CI.A[[2]], "#FF007744")
points(tb, 1 - Fz(tb, mle.weibull$estimate, "B"), pch = 20, cex = 0.1, 
       col = "#7700DD88")
plotCI(ob, CI.B[[1]], CI.B[[2]], "#7700FF44")

# Cox proportional hazards model
fit.cox <- coxph(Surv(exit - entry, status == "Dead") ~ treatment, data = df.obs)
plot(survfit(fit.cox), xlab = "Days", ylab = "Survival")

# Parametric approach ----------------------------------------------------------
fit.weibull <- survreg(Surv(exit - entry, status == "Dead") ~ treatment, 
                       data = df.obs[df.obs$exit - df.obs$entry != 0,], 
                       dist = "weibull")

