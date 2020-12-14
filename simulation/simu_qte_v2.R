#--------------------------------------
# Simulation for quantile treatment effects
#--------------------------------------

task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))

# Load Library
library(survival)
library(dplyr)

# Inverse logit function
inv_logit <- function(x){exp(x) / (1 + exp(x))}

# Set up Simulation Environment

# task_id <- 1
set.seed(task_id)


n <- 400                                                    # Sample size per group
x <- rnorm(n, 0.5, 1)
prob <- inv_logit(-0.25 + 0.5 * x)

a <- rbinom(n, size = 1, prob = prob)
table(a)

y <- rnorm(n, 0.2 + 0.8 * a, sd = 1)
tapply(y, a, mean)

truncation <- 1

# Event time
event_time <- rexp(n, rate = exp(-2 + 0.5 * x) )
summary(event_time)

# Censoring time
censor_time <- rexp(n, 0.2)
# censor_time <- rep(max(event_time) + 1, n)

# Observed time
obs_time = pmin(event_time, censor_time)
status = event_time < censor_time

y_obs <- ifelse(obs_time < truncation, NA, y )

table("less than L" = event_time < truncation, "event" = status )

y_actual <- ifelse(event_time < truncation, min(y) - 1, y)

# Data Sets
db <- data.frame(a, x, y, y_obs, y_actual, event_time, censor_time, obs_time, status, analysis_time = truncation)

count(db, is.na(y_obs), status)

# Propensity socre
fit_propensity <- glm(a ~ x, family = "binomial", data = db)
score <- predict(fit_propensity, type = "response")    # propensity score
db$score <- score

# Determine Rho
db <- db %>% group_by(a) %>%
             mutate(
               y_star = is.na(y_obs) & status, # Subject died with missing value
               rho = sum(y_star) / n())
db <- ungroup(db)

db_summary <- db %>% mutate(group = as.character(a)) %>%
       group_by(group) %>%
       summarise(n = n(),
                 group = unique(as.character(a)),
                 rho = sum(y_star) / n(),
                 pct_missing = mean(is.na(y_obs)),
                 pct_missing_death  = mean(is.na(y_obs) & status) * 100,
                 pct_missing_censor = mean(is.na(y_obs) & (! status)) * 100 )


#############################################
# Modeling
#############################################

# based on observed data only (remove death as well)
est_fun1 <- function(q, xi, db){
  y_std <- (q - db$y_mean) / db$y_sigma
  f_q <- pnorm(y_std)
  term1 <- f_q - xi
  term2 <- ( (db$y_obs < q) - f_q ) / db$score
  mean( (term1 + term2) / db$km, na.rm = TRUE )
}

# Formula (5)
est_fun2 <- function(q, xi,  db){
  y_std <- (q -  db$y_mean)/ db$y_sigma
  f_q <- db$rho + (1 - db$rho) * pnorm(y_std)
  # f_q <- pnorm(y_std)

  r <- (db$obs_time > db$analysis_time) | (db$obs_time <= db$analysis_time & db$status)
  term1 <- f_q - xi                                   # F(q) - q
  term2 <- ifelse(is.na(db$y_obs), 1, db$y_obs <= q)  # I(Y<=q)
  term2 <- db$y_obs <= q
  term3 <- (term2 - f_q) / db$score                   # (I(Y<=q) - F(q)) / e
  term4 <- r / db$km                                  # R / pi

  res <- term4 * (term1 + term3)

  mean(res, na.rm = TRUE)
}

# Formula (6)
est_fun3 <- function(q, xi,  db){
  y_std <- (q -  db$y_mean)/ db$y_sigma
  f_q <- db$rho + (1 - db$rho) * pnorm(y_std)

  r <- db$obs_time > db$analysis_time | (db$obs_time <= db$analysis_time & db$status)
  term1 <- f_q - xi                                   # F(q) - q
  term2 <- ifelse(is.na(db$y_obs), 1, db$y_obs <= q)  # I(Y<=q)
  term2 <- db$y_obs <= q
  term3 <- (term2 - f_q) / db$score                   # (I(Y<=q) - F(q)) / e
  term4 <- r / db$km                                  # R / pi

  res <- term4 * term3 + term1

  mean(res, na.rm = TRUE)
}

qte <- function(db, xi){

  # Censoring K-M Estimator
  # fit_km_censor <- survfit(Surv(obs_time, 1 - status) ~ 1, data = db)
  # db$km <- stepfun(fit_km_censor$time, c(1, fit_km_censor$surv))(pmin(db$obs_time, db$analysis_time))

  # Cox model
  fit_km_censor <- coxph(Surv(obs_time, 1 - status) ~ x, data = db)
  newdata <- data.frame(obs_time = pmin(db$obs_time, db$analysis_time), status = db$status, x = db$x)
  db$km <- survival:::predict.coxph(fit_km_censor, newdata = newdata, type = "survival")

  db$km[is.na(db$km)] <- 1
  # Fit a parametric model
  fit_lm_obs <- lm(y_obs ~ x, data = db)
  db$y_mean  <- predict(fit_lm_obs, newdata = db)
  db$y_sigma <- summary(fit_lm_obs)$sigma

  data.frame(
    quantile_y   = quantile(db$y, probs = xi),
    quantile_all = quantile(db$y_actual, probs = xi),
    quantile_obs = quantile(db$y_obs, probs = xi, na.rm = TRUE),
    naive = uniroot(est_fun1, interval = c(-5, 5), xi = xi, db = db)$root,
    ipw   = uniroot(est_fun2, interval = c(-5, 5), xi = xi, db = db)$root,
    aipw  = uniroot(est_fun3, interval = c(-5, 5), xi = xi, db = db)$root
  )

}

q0 <- qte(db = subset(db, a == 0), xi = 0.5)
q1 <- qte(db = subset(db, a == 1), xi = 0.5)
qdiff <- q1 - q0
res <- bind_rows(
  q0, q1, qdiff
)
rownames(res) <- NULL

res$group = c("0", "1", "diff")

res <- right_join(db_summary, res)

# Save Simulation Results
filename <- paste0(task_id,".Rdata")
save(res, db, file = filename)

#----------------------
# HPC code Submission
#----------------------

# cd /SFS/scratch/zhanyilo/qte
# module add R/4.0.2
# rm *
# qsub -t 1:1000 ~/runr.sh ~/qte/simulation/simu_qte_v2.R

#----------------------
# Simulate Summary
#----------------------
# library(dplyr)
#
# path <- "/SFS/scratch/zhanyilo/qte_v2/"
#
# result <- list()
# for(i in 1:1000){
#   load(file.path(path, paste0(i, ".Rdata")))
#   try(
#     result[[i]] <- res
#   )
# }
#
# result <- bind_rows(result)
#
# t1 <- result %>% group_by(group) %>%
#                  summarise_all(mean)
# save(db, t1, file = "simulation/simu_qte_v2.Rdata")
