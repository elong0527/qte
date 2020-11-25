#--------------------------------------
# Simulation for quantile treatment effects
#--------------------------------------

task_id <- as.integer(Sys.getenv("SGE_TASK_ID"))

# Load Library
library(survival)
library(dplyr)

# Set up Simulation Environment

# task_id <- 1
set.seed(task_id)

n <- 100                                                    # Sample size per group
a <- rep(c(0,1), each = n)                                  # Treatment group (balanced design)
x <- 0.5 * a + rnorm(2 * n)                                 # Covariate for propensity score

# Qaulitiy of life score after transformation
y_0 <- rnorm(n, mean = 0.6)                                   # Continuous outcome when a = 0
y_1 <- rnorm(n, mean = 1  )                                 # Continuous outcome when a = 1
y   <- c(y_0, y_1)
cut_off <- -0.5

y_obs <- ifelse(y < cut_off, NA, y)

truncation <- 1                                              # Truncation time L
cond_time_0  <- rexp(n, 1) * 5                               # Conditional Event time after truncation when a = 0
cond_time_1  <- rexp(n, 0.6) * 5                             # Conditional Event time after truncation  when a = 1
cond_time    <- c(cond_time_0, cond_time_1)

# Event time
event_time    <- ifelse(is.na(y_obs), runif(2 * n, max = truncation), truncation + cond_time)

# Censoring time
censor_time  <- rexp(n, 1) * 5                       # Censoring time

# Observed time
obs_time = pmin(event_time, censor_time)
status = event_time < censor_time

# Data Sets
db <- data.frame(a, x, y, y_obs, event_time, censor_time, obs_time, status, analysis_time = truncation)

count(db, is.na(y_obs), status)

# Propensity socre
fit_propensity <- glm(a ~ x, family = "binomial", data = db)
score <- predict(fit_propensity, type = "response")    # propensity score
db$score <- score

# Determine Rho
db <- db %>% group_by(a) %>%
             mutate(
               y_star = is.na(y_obs) & status, # Subject died with missing value
               rho = sum(y_star) / n)

db_summary <- db %>% mutate(group = as.character(a)) %>%
       group_by(group) %>%
       summarise(group = unique(as.character(a)),
                 rho = sum(y_star) / n,
                 pct_missing_death  = mean(is.na(y_obs) & status),
                 pct_missing_censor = mean(is.na(y_obs) & (! status)) )
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
  f_q <- pmax(db$rho, pnorm(y_std))
  # f_q <- pnorm(y_std)

  r <- db$obs_time > db$analysis_time | (db$obs_time <= db$analysis_time & db$status)
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
  f_q <- pmax(db$rho, pnorm(y_std))

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
  fit_km_censor <- survfit(Surv(obs_time, 1 - status) ~ 1, data = db)
  db$km <- stepfun(fit_km_censor$time, c(1, fit_km_censor$surv))(pmin(db$obs_time, db$analysis_time))

  # Fit a parametric model
  fit_lm_obs <- lm(y_obs ~ x, data = db)
  db$y_mean  <- predict(fit_lm_obs, newdata = db)
  db$y_sigma <- summary(fit_lm_obs)$sigma

  data.frame(
    quantile_all = quantile(db$y, probs = xi),
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
# rm *
# cp ~/qte/simulation/*.R .
# module add R/4.0.2
# qsub -t 1:2 ~/runr.sh simu_qte.R

