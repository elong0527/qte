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

# # Simplest version
# est_fun <- function(q, xi, db){
#     y_std <- (q - db$y_mean) / db$y_sigma
#     f_q <- db$rho + (1 - db$rho) * pnorm(y_std)
#     term <- f_q - xi
#     mean(term, na.rm = TRUE )
# }

# Add propensity score
est_fun <- function(q, xi, db, use_propensity = TRUE, use_km = TRUE, km_rho = c("rho", "rho_km", "rho_km_simple"), type = c("ipw", "aipw")){

  type <- match.arg(type, c("ipw", "aipw"))
  km_rho <- match.arg(km_rho, c("rho", "rho_km", "rho_km_simple"))

  y_std <- (q - db$y_mean) / db$y_sigma

  rho <- db[[km_rho]]

  f_q <- rho + (1 - rho) * pnorm(y_std)

  term <- f_q - xi

  # Propensity score weight
  if(use_propensity){
    term2 <- ifelse(is.na(db$y_obs), 1, db$y_obs <= q)  # I(Y<=q)
    term_propensity <- (term2 - f_q) / db$score                   # (I(Y<=q) - F(q)) / e
  }else{
    term_propensity <- 0
  }

  # KM weight
  if(use_km){
    term_km <- db$r / db$km
  }else{
    term_km <- 1
  }

  if(type == "ipw"){
    res <- mean( (term + term_propensity) * term_km, na.rm = TRUE )
  }

  if(type == "aipw"){
    res <- mean( term + (term_propensity) * term_km, na.rm = TRUE )
  }

  res
}


qte1 <- function(db, xi, use_propensity, use_km, km_rho, type){

  est <- uniroot(est_fun, interval = c(-5, 5), xi = xi, db = db,
                 use_propensity = use_propensity,
                 use_km = use_km,
                 km_rho = km_rho,
                 type = type)$root


  est

}

simu_data <- function(n, censoring = TRUE){
  # n <- 400

  x <- rnorm(n, 0.5, 1)
  prob <- inv_logit(-0.25 + 0.5 * x)

  a <- rbinom(n, size = 1, prob = prob)
  table(a)

  y <- rnorm(n, 0.4 + 0.8 * a + 0.3 * x, sd = 1)

  # Event time
  truncation <- 1
  event_time <- rexp(n, rate = exp(-2 + 0.5 * x) )

  # Censoring time
  if(censoring){
    censor_time <- rexp(n, rate = exp(-2 + 0.4 * x))
  }else{
    censor_time <- max(event_time) + 1
  }

  # Observed time
  obs_time = pmin(event_time, censor_time)
  status = event_time < censor_time

  death     <- (obs_time <  truncation) & status
  y_missing <- obs_time < truncation

  # Observed data
  y_obs <- ifelse(y_missing, NA, y)

  # Actual data
  y_actual <- ifelse(death, min(y, na.rm = TRUE) - 1, y)

  db <- data.frame(a, y, y_obs, y_actual, death,  y_missing, x,
                   event_time, censor_time, obs_time, status, analysis_time = truncation)
}


db_prepare <- function(db){

  # Add information
  db <- db %>% group_by(a) %>%
    mutate(
      y_star = is.na(y_obs) & death, # Subject died with missing value
      rho = sum(y_star) / n())

  # Linear model
  db <- db %>% group_by(a) %>%
    do({
      fit = lm(y_obs ~ x, data = .)
      data.frame(., y_mean = predict(fit, newdata = .),
                 y_sigma = summary(fit)$sigma )

    })

  # Propensity socre
  fit_propensity <- glm(a ~ x, family = "binomial", data = db)
  score <- predict(fit_propensity, type = "response")    # propensity score
  db$score <- score

  # Cox model
  db <- db %>% group_by(a) %>%
    do({
      tmp <- .
      fit_km_censor <- coxph(Surv(obs_time, 1 - status) ~ x, data = tmp)
      newdata <- data.frame(obs_time = pmin(tmp$obs_time, tmp$analysis_time), status = tmp$status, x = tmp$x)
      km <- survival:::predict.coxph(fit_km_censor, newdata = newdata, type = "survival")

      # Surv from Cox model
      fit_km <- coxph(Surv(obs_time, status) ~ x, data = tmp)
      newdata_rho <- data.frame(obs_time = tmp$analysis_time, status = tmp$status, x = tmp$x)
      rho_km <- 1 - survival:::predict.coxph(fit_km, newdata = newdata_rho, type = "survival")

      # Surv from KM estimator
      fit_km <- survfit(Surv(obs_time, status) ~ 1, data = tmp)
      surv   <- data.frame(time = fit_km$time, surv = fit_km$surv)
      rho_km_simple <- 1 - tail(subset(surv, time < 1)$surv, 1)

      data.frame(., km, rho_km, rho_km_simple)
    })


  # R
  db$r <- (db$obs_time > db$analysis_time) | (db$obs_time <= db$analysis_time & db$status)

  db
}



summary_simu_data <- function(db, xi){
  # Summary of statistics
  db_summary <- db %>% group_by(a) %>%
    summarise(n = n(),
              pct_missing = mean(is.na(y_obs)) * 100,
              pct_missing_death  = mean(is.na(y_obs) & status) * 100,
              pct_missing_censor = mean(is.na(y_obs) & (! status)) * 100,
              xi = xi,
              q = quantile(y_actual, probs = xi)) %>%
    ungroup() %>%
    mutate(diff = diff(q))

  db_summary
}



# Simulation for censoring situation
res <- list()
truth <- list()
for(i in 1:1){
  n <- 400
  # Simulate data
  db <- simu_data(n, censoring = TRUE)
  db <- db_prepare(db)

  # Simulation Set-up
  par1 <- expand.grid(xi = 0.5,
                     use_propensity = c(TRUE, FALSE),
                     use_km = c(TRUE, FALSE),
                     km_rho = c("rho", "rho_km", "rho_km_simple"),
                     type = c("ipw"), stringsAsFactors = FALSE)

  par2 <- expand.grid(xi = 0.5,
                      use_propensity = c(TRUE),
                      use_km = c(TRUE),
                      km_rho = c("rho", "rho_km", "rho_km_simple"),
                      type = c("aipw"), stringsAsFactors = FALSE)

  par <- subset(bind_rows(par1, par2), ! (use_propensity & ! use_km))

  # Fit simulated data
  res_tmp <- purrr:::pmap(par, function(xi, use_propensity, use_km, km_rho, type){
    print(use_propensity)
    x0 <- qte1(subset(db, a == 0), xi = xi, use_propensity = use_propensity, use_km = use_km, km_rho = km_rho, type = type)
    x1 <- qte1(subset(db, a == 1), xi = xi, use_propensity = use_propensity, use_km = use_km, km_rho = km_rho, type = type)
    diff <- x1 - x0
    data.frame(xi, use_propensity, use_km, km_rho, type, x0, x1, diff, row.names = NULL)
  })

  true_tmp <- purrr:::pmap(data.frame(xi = unique(par$xi)), function(xi){
    summary_simu_data(simu_data(10000), xi = xi)
  })

  res[[i]] <- bind_rows(res_tmp)
  truth[[i]] <- bind_rows(true_tmp)
}

res_1 <- bind_rows(res)
truth <- bind_rows(truth)
res_1$censoring <- TRUE
res_1$n <- 400

# Simulation for no-censoring situation (N = 400)
res <- list()
for(i in 1:1){
  n <- 400
  # Simulate data
  db <- simu_data(n, censoring = FALSE)
  db <- db_prepare(db)

  # Simulation Set-up
  par <- expand.grid(xi = 0.5,
                      use_propensity = c(TRUE, FALSE),
                      use_km = c(FALSE),
                      km_rho = c("rho", "rho_km", "rho_km_simple"),
                      type = c("ipw"), stringsAsFactors = FALSE)

  # Fit simulated data
  res_tmp <- purrr:::pmap(par, function(xi, use_propensity, use_km, km_rho, type){
    print(use_propensity)
    x0 <- qte1(subset(db, a == 0), xi = xi, use_propensity = use_propensity, use_km = use_km, km_rho = km_rho, type = type)
    x1 <- qte1(subset(db, a == 1), xi = xi, use_propensity = use_propensity, use_km = use_km, km_rho = km_rho, type = type)
    diff <- x1 - x0
    data.frame(xi, use_propensity, use_km, km_rho, type, x0, x1, diff, row.names = NULL)
  })

  res[[i]] <- bind_rows(res_tmp)
}

res_2 <- bind_rows(res)
truth_2 <- bind_rows(truth)
res_2$censoring <- FALSE
res_2$n <- 400

# Simulation for no-censoring situation (N = 10000)
res <- list()
for(i in 1:1){
  n <- 10000
  # Simulate data
  db <- simu_data(n, censoring = FALSE)
  db <- db_prepare(db)

  # Simulation Set-up
  par <- expand.grid(xi = 0.5,
                     use_propensity = c(TRUE, FALSE),
                     use_km = c(FALSE),
                     km_rho = c("rho", "rho_km", "rho_km_simple"),
                     type = c("ipw"), stringsAsFactors = FALSE)

  # Fit simulated data
  res_tmp <- purrr:::pmap(par, function(xi, use_propensity, use_km, km_rho, type){
    print(use_propensity)
    x0 <- qte1(subset(db, a == 0), xi = xi, use_propensity = use_propensity, use_km = use_km, km_rho = km_rho, type = type)
    x1 <- qte1(subset(db, a == 1), xi = xi, use_propensity = use_propensity, use_km = use_km, km_rho = km_rho, type = type)
    diff <- x1 - x0
    data.frame(xi, use_propensity, use_km, km_rho, type, x0, x1, diff, row.names = NULL)
  })

  res[[i]] <- bind_rows(res_tmp)
}

res_3 <- bind_rows(res)
res_3$censoring <- FALSE
res_3$n <- 10000

res <- bind_rows(res_1, res_2, res_3)


# Save Simulation Results
filename <- paste0(task_id,".Rdata")
save(res, truth, file = filename)

#----------------------
# HPC code Submission
#----------------------

# cd /SFS/scratch/zhanyilo/qte3
# module add R/4.0.2
# rm *
# qsub -t 1:1000 ~/runr.sh ~/qte/simulation/simu_qte_v3.R

#----------------------
# Simulate Summary
#----------------------
# library(dplyr)
#
# path <- "/SFS/scratch/zhanyilo/qte3/"
#
# result <- list()
# result_truth <- list()
#
# for(i in 1:1000){
#   load(file.path(path, paste0(i, ".Rdata")))
#   try({
#     result[[i]] <- res
#     result_truth[[i]] <- truth
#   })
# }
#
# result <- bind_rows(result)
# result_truth <- bind_rows(result_truth)
#
# truth_summary <- bind_rows(result_truth) %>% group_by(a, xi) %>%
#   summarise_all(mean)
#
# # Method 1
# # truth_est <- data.frame(truth_summary$q[1], truth_summary$q[2], unique(truth_summary$diff))
# # names(truth_est) <- c("x0", "x1", "diff")
#
# # Method 2
# truth_est <- bind_rows(result ) %>% subset(use_propensity == FALSE & use_km == FALSE & km_rho == "rho_km" & n == 10000) %>%
#                                     group_by(xi, use_propensity, use_km, km_rho, type, censoring, n) %>%
#                                     summarise(x0 = mean(x0), x1 = mean(x1), diff = mean(diff))
#
# res_summary <- bind_rows(result ) %>% group_by(xi, use_propensity, use_km, km_rho, type, censoring, n) %>%
#                    summarise(x0 = mean(x0), x0_true = truth_est$x0, x0_rmse = sqrt(mean(x0 - truth_est$x0)^2),
#                              x1 = mean(x1), x1_true = truth_est$x1, x1_rmse = sqrt(mean(x1 - truth_est$x1)^2),
#                              diff = mean(diff), diff_true = truth_est$diff, diff_rmse = sqrt(mean(diff - truth_est$diff)^2))
#
# res_summary
#
# save(db, res_summary, truth_summary, truth_est, file = "simulation/simu_qte_v3.Rdata")

