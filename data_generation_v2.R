# This is a preliminary script to practice simulating survival data.

# We use an existing data set to examine how the actual data set could look like
# and to practice applying survival analysis methods.

library(survival)
library(biostat3)
library(dplyr)
set.seed(505) # for reproducibility

# Data -------------------------------------------------------------------------
# We adapt the data set on colon cancer to fit our problem setting.
# We assume that the event of interest is death from any cause.
df <- biostat3::colon[c("id", "sex", "status", "dx", "exit")] %>%
  mutate(
    days = as.numeric(exit - dx),
    censored = status %in% c("Alive", "Lost to follow-up"),
    country = factor(rbinom(15564, 1, 1/3), levels = 0:1, 
                     labels = c("Sweden", "Norway")),
    y.exit = biostat3::year(exit, trunc = TRUE),
    m.exit = months(exit) 
    )

# We visualize the study time by plotting the number of days survived
# against the date of diagnosis.
status_colors <- function(v) {
  case_match(v,
             "Alive" ~ "#30AE2240",
             c("Dead: cancer", "Dead: other") ~ "#00000060",
             "Lost to follow-up" ~ "#F0009F")
}
plot(df$dx, df$days, pch = 20, col = status_colors(df$status), 
     main = "Days survived by diagnosis date",
     xlab = "Diagnosis date",
     ylab = "Days")

# Kaplan-Meier curve for the whole data set.
plot(survfit(Surv(days, !censored) ~ 1, data = df, conf.type = "log-log"), 
     main = "Kaplan-Meier survival curve \n(combined data)",
     xlab = "Days",
     ylab = "Survival")

# Kaplan-Meier curves by country
plot(survfit(Surv(days, !censored) ~ 1, 
             data = df[df$country == "Norway",], 
             conf.type = "log-log"),
     main = "Kaplan-Meier survival curves by country",
     xlab = "Days",
     ylab = "Survival",
     col = "red")
lines(survfit(Surv(days, !censored) ~ 1, 
              data = df[df$country == "Sweden",], 
              conf.type = "log-log"),
      col = "blue")
legend("topright", lty = c(1, 1), col = c("red", "blue"), bty = "n", 
       legend = c("Norway", "Sweden"))

# Cox proportional hazards
fitcox <- coxph(Surv(days, !censored) ~ country, data = df)
summary(fitcox)  # coef -0.01664
plot(survfit(fitcox), xlab = "Days", ylab = "Survival probability")

# Modelling the problem --------------------------------------------------------
# Since we are assuming we have interval-censored data, we need to define the
# start- and end-points of the intervals.

# We create helper functions to help reformat the dates.
month_start <- function(ymd) as.Date(gsub("[0-9]{2}$", "01", ymd))
month_end <- function(ymd) month_start(ymd + 32) - 1

# Modify data sets so that data from Norway is only accurate to the month, 
# while data from Sweden is accurate to the day.
norway <- df %>% filter(country == "Norway") %>%
  mutate(
    exit.start = as.Date(paste0(1, m.exit, y.exit), "%d%b%Y"),
    exit.end = month_end(exit.start),
    days.min = as.numeric(exit.start - dx),
    days.max = as.numeric(exit.end - dx)
  )

sweden <- df %>% filter(country == "Sweden") %>%
  mutate(exit.start = exit, exit.end = exit, days.min = days, days.max = days)

dfmix <- rbind(sweden, norway) %>% select(-c(days, exit)) %>%
  mutate(
    exit.end = as.Date(ifelse(censored, NA, exit.end)),
    days.max = ifelse(censored, NA, days.max)
  )

# Midpoint imputation ----------------------------------------------------------
dfmix <- dfmix %>% 
  mutate(
    days.mid = ifelse(country == "Sweden", days.min, 
                      ifelse(!is.na(days.max), floor((days.max + days.min)/2), 
                             days.min)),
    exit.mid = dx + days.mid)

# Compare Kaplan-Meier curves
plot(survfit(Surv(days, !censored) ~ country, data = df, conf.type = "log-log"), 
     main = "Kaplan-Meier survival curve \n(combined data)",
     xlab = "Days",
     ylab = "Survival")
lines(survfit(Surv(days.mid, !censored) ~ country, data = dfmix, conf.type = "log-log"), 
      col = "red")

# Compare Cox proportional hazard models
fitcox.mid <- coxph(Surv(days.mid, !censored) ~ country, data = dfmix)
summary(fitcox.mid)  # coef -0.01269
plot(survfit(fitcox.mid), xlab = "Days", ylab = "Survival probability")
lines(survfit(fitcox), col="red")

# Proposed approach ------------------------------------------------------------

x <- seq(0, 8000, by=100)
hist(df$days, breaks=c(0,x+50), freq=FALSE, ylim=c(0,0.003))


# Define coarsening limits
l <- dfmix$days.min
u <- ifelse(dfmix$censored, Inf, dfmix$days.max) + 1

# Helper functions
fy <- function(l, u, b) Fz(u, b) - Fz(l, b)
logL <- function(b, l, u) -sum(log(fy(l, u, b)))

# Try different distributions
Fz <- function(z, b) pweibull(z, b[1], b[2])
-nlm(logL, c(0.63, 1000), l, u)$minimum

Fz <- function(z, b) pexp(z, b)
-nlm(logL, 1/1000, l, u)$minimum


# Fit the distribution with the highest likelihood
mle <- nlm(logL, 1/1000, l, u, hessian=TRUE)
plot(x, 1-pexp(x, mle$estimate))

