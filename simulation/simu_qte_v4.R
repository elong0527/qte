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
est_fun <- function(q, xi, db,
                    km_rho = c("rho", "rho_km", "rho_km_simple"),
                    type = c("or", "propensity","ipw1", "aipw1", "ipw2", "aipw2")){

  type <- match.arg(type, c("or", "propensity","ipw1", "aipw1", "ipw2", "aipw2"))
  km_rho <- match.arg(km_rho, c("rho", "rho_km", "rho_km_simple"))

  y_std <- (q - db$y_mean) / db$y_sigma

  rho <- db[[km_rho]]

  f_q <- rho + (1 - rho) * pnorm(y_std)

  term <- f_q - xi

  # KM weight
  weight_km <- ifelse(is.na(db$km), db$r, db$r / db$km)
  weight_km <- weight_km / mean(weight_km[db$r])

  # Propensity score weight
  weight_propensity <- db$r / db$score
  weight_propensity <- weight_propensity / mean(weight_propensity[db$r])

  term2 <- ifelse(is.na(db$y_obs), 1, db$y_obs <= q)  # I(Y<=q)
  # term_propensity <- (term2 * weight_km - db$r * f_q) * weight_propensity
  term_propensity1 <- (term2 - f_q) * weight_km * weight_propensity
  term_propensity2 <- ((db$y_obs <= q) - pnorm(y_std)) * weight_propensity

  if(type == "or"){
    res <- mean(term, na.rm = TRUE)
  }

  if(type == "propensity"){
    # res <- mean(term2 * weight_km * weight_propensity, na.rm = TRUE) - xi
    xi_1 <- (xi - rho) / (1 - rho)
    xi_1 <- mean(xi_1[! is.na(db$y_obs)])
    res <- mean((db$y_obs <= q) * weight_km * weight_propensity, na.rm = TRUE) - mean(xi_1)
  }

  if(type == "ipw1"){
    res <- mean( term * weight_km + term_propensity1, na.rm = TRUE )
  }

  if(type == "aipw1"){
    res <- mean( term + term_propensity1, na.rm = TRUE )
  }

  if(type == "ipw2"){
    res <- mean( term * weight_km + term_propensity2, na.rm = TRUE )
  }

  if(type == "aipw2"){
    res <- mean( term + term_propensity2, na.rm = TRUE )
  }

  res
}


qte1 <- function(db, xi, km_rho, type){

  est <- uniroot(est_fun, interval = c(-5, 5), xi = xi, db = db,
                 km_rho = km_rho,
                 type = type)$root


  est

}

simu_data <- function(n,
                      censoring = TRUE,
                      misspecify_censor = FALSE,
                      misspecify_propensity = FALSE,
                      misspecify_or = FALSE){
  # n <- 400

  x <- rnorm(n, 0.5, 1)
  x2 <- rnorm(n, 0, 1)

  if(misspecify_propensity){
    prob <- inv_logit(-0.25 + 0.5 * x + x * x2)
  }else{
    prob <- inv_logit(-0.25 + 0.5 * x)
  }

  a <- rbinom(n, size = 1, prob = prob)
  table(a)

  if(misspecify_or){
    y <- rnorm(n, 0.4 + 0.8 * a + 0.3 * x + x * x2, sd = 1)
  }else{
    y <- rnorm(n, 0.4 + 0.8 * a + 0.3 * x, sd = 1)
  }


  # Event time
  truncation <- 1
  event_time <- rexp(n, rate = exp(-2 + 0.5 * x) )

  # Censoring time
  if(censoring){

    # Censoring mis-specification
    if(misspecify_censor){
      censor_time <- rexp(n, rate = exp(-2 + 0.4 * x + x*x2))
    }else{
      censor_time <- rexp(n, rate = exp(-2 + 0.4 * x))
    }

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

  db <- data.frame(a, y, y_obs, y_actual, death,  y_missing, x, x2,
                   event_time, censor_time, obs_time, status, analysis_time = truncation)
}


db_prepare <- function(db, use_x2 = FALSE){

  # Add information
  db <- db %>% group_by(a) %>%
    mutate(
      y_star = is.na(y_obs) & death, # Subject died with missing value
      rho = sum(y_star) / n())

  # Linear model
  if(use_x2){
    db <- db %>% group_by(a) %>%
      do({
        fit = lm(y_obs ~ x + I(x*x2), data = .)
        data.frame(., y_mean = predict(fit, newdata = .),
                   y_sigma = summary(fit)$sigma )

      })
  }else{
    db <- db %>% group_by(a) %>%
      do({
        fit = lm(y_obs ~ x, data = .)
        data.frame(., y_mean = predict(fit, newdata = .),
                   y_sigma = summary(fit)$sigma )

      })

  }

  # Propensity score
  fit_propensity <- glm(a ~ x, family = "binomial", data = db)
  score <- predict(fit_propensity, type = "response")    # propensity score
  score <- ifelse(db$a == 0, 1 - score, score)
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


run_simu <- function(db, par){
  n <- nrow(db)
  res_tmp <- purrr:::pmap(par, function(xi, km_rho, type){
    print(type)
    x0 <- qte1(subset(db, a == 0), xi = xi, km_rho = km_rho, type = type)
    x1 <- qte1(subset(db, a == 1), xi = xi, km_rho = km_rho, type = type)
    diff <- x1 - x0
    data.frame(n, xi, km_rho, type, x0, x1, diff, row.names = NULL)
  })
  bind_rows(res_tmp)
}

library(purrr)
xi <- seq(0.35, 0.95, 0.05)

######## Simulation for true value ####################
par <- expand.grid(xi = xi,
                   km_rho = c("rho_km"),
                   type = c("or"), stringsAsFactors = FALSE)
db <- simu_data(10000, censoring = FALSE) %>% db_prepare(use_x2 = TRUE)

meta <- map_df(xi, summary_simu_data, db = db)
res0 <- run_simu(db, par)


# db <- simu_data(10000, censoring = FALSE) %>% db_prepare(use_x2 = TRUE)
#
# meta <- map_df(xi, summary_simu_data, db = db)
# res0 <- run_simu(db, par)

######## Simulation with different Rho estimation ######
par <- expand.grid(xi = xi,
                    km_rho = c("rho", "rho_km", "rho_km_simple"),
                    type = c("or", "propensity", "ipw1", "aipw1", "ipw2", "aipw2"), stringsAsFactors = FALSE)

censoring <- TRUE
db <- simu_data(400, censoring = censoring) %>% db_prepare()
res1 <- run_simu(db, par)

db <- simu_data(400, censoring = censoring, misspecify_censor = TRUE) %>% db_prepare()
res2 <- run_simu(db, par)

db <- simu_data(400, censoring = censoring, misspecify_propensity = TRUE) %>% db_prepare()
res3 <- run_simu(db, par)

db <- simu_data(400, censoring = censoring, misspecify_or = TRUE) %>% db_prepare()
res4 <- run_simu(db, par)

# Save Simulation Results
filename <- paste0(task_id,".Rdata")
save(meta, res0, res1, res2, res3, res4, file = filename)

#----------------------
# HPC code Submission
#----------------------

# cd /SFS/scratch/zhanyilo/qte3
# module add R/4.0.2
# rm *
# qsub -t 1:1000 ~/runr.sh ~/qte/simulation/simu_qte_v4.R

#----------------------
# Simulate Summary
#----------------------
# library(dplyr)
#
# path <- "/SFS/scratch/zhanyilo/qte5"
#
# r0 <- list()
# r1 <- list()
# r2 <- list()
# r3 <- list()
# r4 <- list()
# for(i in 1:1000){
#   load(file.path(path, paste0(i, ".Rdata")))
#   try({
#     r0[[i]] <- res0
#     r1[[i]] <- res1
#     r2[[i]] <- res2
#     r3[[i]] <- res3
#     r4[[i]] <- res4
#   })
# }
#
# r0 <- bind_rows(r0)
# r1 <- bind_rows(r1)
# r2 <- bind_rows(r2)
# r3 <- bind_rows(r3)
# r4 <- bind_rows(r4)
#
# est0 <- r0 %>% group_by(xi) %>%
#   summarise(x0_true = mean(x0), x1_true = mean(x1), diff_true = mean(diff))
#
# summary_est <- function(db){
#   db %>% left_join(est0) %>%
#     group_by(n, xi, km_rho, type) %>%
#     summarise(x0 = mean(x0), x0_true = mean(x0_true), x0_bias = x0 - x0_true, x0_rmse = sqrt(mean(x0 - x0_true)^2),
#               x1 = mean(x1), x1_true = mean(x1_true), x1_bias = x1 - x1_true, x1_rmse = sqrt(mean(x1 - x1_true)^2),
#               diff = mean(diff), diff_true = mean(diff_true), diff_bias = diff - diff_true, diff_rmse = sqrt(mean(diff - diff_true)^2))
# }
#
# est1 <- r1 %>%  summary_est()
# est2 <- r2 %>%  summary_est()
# est3 <- r3 %>%  summary_est()
# est4 <- r4 %>%  summary_est()
#
# #### Visualization
# library(ggplot2)
# library(tidyr)
#
# g_summary <- function(db, scales){
#   est <- db %>% select(- ends_with("true"), - ends_with("rmse")) %>% pivot_longer(cols = c("x0", "x1", "diff"))
#   est_true <- db %>% select(- ends_with(c("x0", "x1", "diff")), - ends_with("rmse")) %>%
#     pivot_longer(cols = c("x0_true", "x1_true", "diff_true"), values_to = "truth") %>%
#     mutate(name = gsub("_true", "", name))
#   est <- left_join(est, est_true)
#   rmse <- db %>% select(- ends_with("true"), - ends_with(c("x0", "x1", "diff"))) %>% pivot_longer(cols = ends_with("rmse"))
#   bias <- db %>% select(- ends_with("true"), - ends_with(c("x0", "x1", "diff"))) %>% pivot_longer(cols = ends_with("bias"))
#
#
#   g_est <- ggplot(data = est ) +
#     geom_point(aes(x = xi, y = value, group = type, color = type)) +
#     geom_line(aes(x = xi, y = truth)) +
#     facet_grid(name ~ km_rho, scales = scales) +
#     ylab("Estimation") +
#     theme_bw()
#
#   g_rmse <- ggplot(data = rmse ) +
#     geom_point(aes(x = xi, y = value, group = type, color = type)) +
#     facet_grid(name ~ km_rho, scales = scales) +
#     ylab("RMSE") +
#     theme_bw()
#
#   g_bias <- ggplot(data = bias ) +
#     geom_point(aes(x = xi, y = value, group = type, color = type)) +
#     facet_grid(name ~ km_rho, scales = scales) +
#     ylab("BIAS") +
#     theme_bw()
#
#   list(g_est = g_est, g_rmse = g_rmse, g_bias = g_bias)
# }
#
# g1 <- g_summary(est1, scales = "free")
# g2 <- g_summary(est2, scales = "free")
# g3 <- g_summary(est3, scales = "free")
# g4 <- g_summary(est4, scales = "free")
#
#
#
# g_summary2 <- function(db){
#
#   est <- db %>% subset(km_rho == "rho_km") %>% pivot_longer(cols = x0:diff_rmse) %>%
#                 separate(name, c("var", "metric")) %>%
#                 mutate(metric = ifelse(is.na(metric), "est", metric)) %>%
#                 subset(metric != "true") %>%
#                 mutate(metric = factor(metric, levels = c("est", "bias", "rmse")))
#
#   est_true <- db %>% select(- ends_with(c("x0", "x1", "diff")), - ends_with("rmse"), - ends_with("bias")) %>%
#     pivot_longer(cols = c("x0_true", "x1_true", "diff_true"), values_to = "truth") %>%
#     mutate(name = gsub("_true", "", name)) %>%
#     rename(var = name)
#   est <- left_join(est, est_true)
#
#   g <- ggplot(data = est) +
#           geom_point(aes(x = xi, y = value, group = type, color = type)) +
#           geom_line(aes(x = xi, y = truth, group = type), data = subset(est, metric == "est")) +
#           facet_wrap(metric ~ var, scales = "free") +
#           ylab("Metric Value") +
#           theme_bw() +
#           scale_color_brewer(palette="Dark2")
#
#   g
# }
#
#
# s1 <- g_summary2(est1)
# s2 <- g_summary2(est2)
# s3 <- g_summary2(est3)
# s4 <- g_summary2(est4)
#
# db <- simu_data(400)
# save(meta, est0, est1, est2, est3, est4, db, g1, g2, g3, g4, s1, s2, s3, s4, file = "simulation/simu_qte_v4.Rdata")
