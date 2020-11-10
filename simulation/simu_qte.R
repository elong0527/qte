library(devtools)
library(dplyr)
devtools::load_all()

logit_inv <- function(x) exp(x) / (1 + exp(x))

#############################################
# Generate Data
#############################################

n <- 200                                                    # Sample size
x <- rnorm(n)                                               # Covariate
a <- rbinom(n, size = 1, logit_inv(0.5 * x))                # Treatment group generation model

# glm(a ~ x, family = "binomial")

y_0 <- rnorm(n, mean = 1)                                   # Continous outcome when a = 0
y_1 <- rnorm(n, mean = 0.6)                                 # Continous outcome when a

time_0 <- rexp(n, 1) * 5                                    # Event time when a = 0
time_1 <- rexp(n, 0.6) * 5                                  # Event time when a = 1
censor <- rexp(n, 1) * 5                                    # Censoring time
truncation <- 1                                             # Truncation time / L


db_full <- data.frame(x, a, y_0, y_1, time_0, time_1, censor)
db_obs  <- data.frame(x, a) %>% mutate(
                                  event_time = ifelse(a == 0, time_0, time_1),
                                  censor_time = censor,
                                  obs_time = pmin(event_time, censor),
                                  status = event_time < censor,
                                  y_all = ifelse(a == 0, y_0, y_1),
                                  y_obs = ifelse(obs_time >= truncation, y_all, NA)
                                  )
fit_propensity <- glm(a ~ x, family = "binomial", data = db_obs)
score <- predict(fit_propensity, type = "response")    # propensity score
db_obs$score <- score

db_obs %>% group_by(a) %>%
           summarise(n_obs = sum(! is.na(y_obs)))


#############################################
# Modeling
#############################################
library(survival)

db <- subset(db_obs, a == 0)
# Censoring K-M Estimator
fit_km_censor <- survfit(Surv(obs_time, 1 - status) ~ 1, data = db)
db$km <- stepfun(fit_km_censor$time, c(1, fit_km_censor$surv))(pmin(db$obs_time, truncation))

fit_lm <- lm(y_all ~ x, data = db)
db$y_mean <- predict(fit_lm)
db$y_sigma <- summary(fit_lm)$sigma

fit_lm_obs <- lm(y_obs ~ x, data = db)
db$y_mean_obs  <- predict(fit_lm_obs, newdata = db)
db$y_sigma_obs <- summary(fit_lm_obs)$sigma

# Formula (1) in Section 3.1
# Based on outcome with no missing value
# Not feasible with censoring
G1_fun <- function(q, xi){
  y_std <- (q -  db$y_mean)/ db$y_sigma
  f_q <- pnorm(y_std)
  term1 <- f_q - xi
  term2 <- ( (db$y_all < q) - f_q ) / db$score
  mean(term1 + term2)
}

# Formula (5) in Section 3.2
G2_fun <- function(q, xi){
  y_std <- (q - db$y_mean_obs) / db$y_sigma_obs
  f_q <- pnorm(y_std)
  term1 <- f_q - xi
  term2 <- ( (db$y_obs < q) - f_q ) / db$score
  mean( (term1 + term2) / db$km, na.rm = TRUE )
}

uniroot(G1_fun, interval = c(-5, 5), xi = 0.5)$root
uniroot(G2_fun, interval = c(-5, 5), xi = 0.5)$root
