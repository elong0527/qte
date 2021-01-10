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

  est <- uniroot(est_fun, interval = c(-20, 20), xi = xi, db = db,
                 km_rho = km_rho,
                 type = type, extendInt = "yes")$root


  est

}

simu_data <- function(n,
                      censoring = TRUE,
                      misspecify_censor = FALSE,
                      misspecify_propensity = FALSE,
                      misspecify_or = FALSE){
  # n <- 400

  x <- rnorm(n, 0.5, 1)
  # x2 <- rnorm(n, 0, 1)

  if(misspecify_propensity){
    prob <- inv_logit(-0.25 + 0.5 * x + (x^2 - 0.25))
  }else{
    prob <- inv_logit(-0.25 + 0.5 * x)
  }

  a <- rbinom(n, size = 1, prob = prob)
  table(a)

  if(misspecify_or){
    y <- rnorm(n, 0.4 + 0.8 * a + 0.3 * x + (x^2 - 0.25), sd = 1)
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
      censor_time <- rexp(n, rate = exp(-2 + 0.4 * x + (x^2 - 0.25)))
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

  db <- data.frame(a, y, y_obs, y_actual, death,  y_missing, x,
                   event_time, censor_time, obs_time, status, analysis_time = truncation)
}


db_prepare <- function(db,
                       linear_formula = y_obs ~ x,
                       propensity_formula = a ~ x,
                       surv_formula = Surv(obs_time, status) ~ x,
                       censor_formula = Surv(obs_time, 1 - status) ~ x){

  # Add information
  db <- db %>% group_by(a) %>%
    mutate(
      y_star = is.na(y_obs) & death, # Subject died with missing value
      rho = sum(y_star) / n())

  # Linear model
  db <- db %>% group_by(a) %>%
      do({
        fit = lm(linear_formula, data = .)
        data.frame(., y_mean = predict(fit, newdata = .),
                   y_sigma = summary(fit)$sigma )

      })

  # Propensity score
  fit_propensity <- glm(propensity_formula, family = "binomial", data = db)
  score <- predict(fit_propensity, type = "response")    # propensity score
  score <- ifelse(db$a == 0, 1 - score, score)
  db$score <- score

  # Cox model
  db_tmp <- list()
  db_tmp[[1]] <- db %>% subset(a == 0)
  db_tmp[[2]] <- db %>% subset(a == 1)

  for(i in 1:2){
    tmp <- db_tmp[[i]]

    if(any(db$status == 0)){
      fit_km_censor <- coxph(censor_formula, data = tmp)
      newdata <- data.frame(obs_time = pmin(tmp$obs_time, tmp$analysis_time), status = tmp$status, x = tmp$x)
      km <- survival:::predict.coxph(fit_km_censor, newdata = newdata, type = "survival")
    }else{
      km <- rep(1, nrow(tmp))
    }

    # Surv from Cox model
    fit_km <- coxph(surv_formula, data = tmp)
    newdata_rho <- data.frame(obs_time = tmp$analysis_time, status = tmp$status, x = tmp$x)
    rho_km <- 1 - survival:::predict.coxph(fit_km, newdata = newdata_rho, type = "survival")

    # Surv from KM estimator
    fit_km <- survfit(Surv(obs_time, status) ~ 1, data = tmp)
    surv   <- data.frame(time = fit_km$time, surv = fit_km$surv)
    rho_km_simple <- 1 - tail(subset(surv, time < 1)$surv, 1)

    db_tmp[[i]] <- data.frame(tmp, km, rho_km, rho_km_simple)
  }
  db <- do.call(rbind, db_tmp)



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
xi <- seq(0.35, 0.85, 0.05)
# xi <- c(0.5, 0.7)
######## Simulation for true value ####################
par <- expand.grid(xi = xi,
                   km_rho = c("rho_km"),
                   type = c("or"), stringsAsFactors = FALSE)

meta_true <- list()
res_true <- list()
db <- simu_data(10000, censoring = FALSE) %>% db_prepare()
meta_true[[1]] <- map_df(xi, summary_simu_data, db = db)
res_true[[1]] <- run_simu(db, par)

db <- simu_data(10000, censoring = FALSE, misspecify_or = TRUE) %>% db_prepare(linear_formula = y_obs ~ x + x^2)
meta_true[[2]] <- map_df(xi, summary_simu_data, db = db)
res_true[[2]] <- run_simu(db, par)

db <- simu_data(10000, censoring = FALSE, misspecify_propensity = TRUE) %>% db_prepare(propensity_formula = a ~ x + x^2)
meta_true[[3]] <- map_df(xi, summary_simu_data, db = db)
res_true[[3]] <- run_simu(db, par)

db <- simu_data(10000, censoring = FALSE, misspecify_censor = TRUE) %>% db_prepare(censor_formula = Surv(obs_time, 1 - status) ~ x + x^2)
meta_true[[4]] <- map_df(xi, summary_simu_data, db = db)
res_true[[4]] <- run_simu(db, par)



######## Simulation with different Rho estimation ######
par <- expand.grid(xi = xi,
                    km_rho = c("rho", "rho_km", "rho_km_simple"),
                    type = c("or", "propensity", "ipw1", "aipw1", "ipw2", "aipw2"), stringsAsFactors = FALSE)

censoring <- TRUE
n <- 400
res <- list()
meta <- list()
db <- simu_data(n, censoring = censoring) %>% db_prepare()
meta[[1]] <- map_df(xi, summary_simu_data, db = db)
res[[1]] <- run_simu(db, par)

db <- simu_data(n, censoring = censoring, misspecify_or = TRUE) %>% db_prepare()
meta[[2]] <- map_df(xi, summary_simu_data, db = db)
res[[2]] <- run_simu(db, par)

db <- simu_data(n, censoring = censoring, misspecify_propensity = TRUE) %>% db_prepare()
meta[[3]] <- map_df(xi, summary_simu_data, db = db)
res[[3]] <- run_simu(db, par)

db <- simu_data(n, censoring = censoring, misspecify_censor = TRUE) %>% db_prepare()
meta[[4]] <- map_df(xi, summary_simu_data, db = db)
res[[4]] <- run_simu(db, par)

# Save Simulation Results
filename <- paste0(task_id,".Rdata")
save(meta, meta_true, res, res_true, file = filename)

#----------------------
# HPC code Submission
#----------------------

# cd /SFS/scratch/zhanyilo/qte3
# module add R/4.0.2
# rm *
# qsub -t 1:1000 ~/runr.sh ~/qte/simulation/simu_qte_v5.R

# ----------------------
# Simulation Summary
# ----------------------
# library(dplyr)
#
# path <- "/SFS/scratch/zhanyilo/qte5"
#
# result <- list()
# result_true <- list()
# result_meta <- list()
# result_meta_true <- list()
# for(i in 1:1000){
#   try({
#     load(file.path(path, paste0(i, ".Rdata")))
#     result[[i]] <- bind_rows(res, .id = "scenario")
#     result_true[[i]] <- bind_rows(res, .id = "scenario")
#     result_meta[[i]] <- bind_rows(meta, .id = "scenario")
#     # result_meta_true[[i]] <- bind_rows(meta_true .id = "scenario")
#   })
# }
#
# result <- bind_rows(result)
# result_true <- bind_rows(result_true)
# result_meta <- bind_rows(result_meta)
#
# # Meta information
# result_meta <- result_meta %>% select(- xi, -q, -diff) %>%
#   group_by(scenario, a) %>%
#   summarise_all(mean)
#
# # True value
# est0 <- result_true %>% group_by(scenario, xi) %>%
#   summarise(x0_true = mean(x0), x1_true = mean(x1), diff_true = mean(diff))
#
# summary_est <- function(db){
#   db %>% left_join(est0) %>%
#     group_by(scenario, n, xi, km_rho, type) %>%
#     summarise(x0 = mean(x0), x0_true = mean(x0_true), x0_bias = x0 - x0_true, x0_rmse = sqrt(mean(x0 - x0_true)^2),
#               x1 = mean(x1), x1_true = mean(x1_true), x1_bias = x1 - x1_true, x1_rmse = sqrt(mean(x1 - x1_true)^2),
#               diff = mean(diff), diff_true = mean(diff_true), diff_bias = diff - diff_true, diff_rmse = sqrt(mean(diff - diff_true)^2))
# }
#
# est1 <- result %>%  summary_est()
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
# g1 <- subset(est1, scenario == 1) %>% g_summary(scales = "free")
# g2 <- subset(est1, scenario == 2) %>% g_summary(scales = "free")
# g3 <- subset(est1, scenario == 3) %>% g_summary(scales = "free")
# g4 <- subset(est1, scenario == 4) %>% g_summary(scales = "free")
#
#
#
# g_summary2 <- function(db){
#
#   est <- db %>% subset(km_rho == "rho_km") %>% pivot_longer(cols = x0:diff_rmse) %>%
#     separate(name, c("var", "metric")) %>%
#     mutate(metric = ifelse(is.na(metric), "est", metric)) %>%
#     subset(metric != "true") %>%
#     mutate(metric = factor(metric, levels = c("est", "bias", "rmse")))
#
#   est_true <- db %>% select(- ends_with(c("x0", "x1", "diff")), - ends_with("rmse"), - ends_with("bias")) %>%
#     pivot_longer(cols = c("x0_true", "x1_true", "diff_true"), values_to = "truth") %>%
#     mutate(name = gsub("_true", "", name)) %>%
#     rename(var = name)
#   est <- left_join(est, est_true)
#
#   g <- ggplot(data = est) +
#     geom_line(aes(x = xi, y = value, group = type, color = type)) +
#     geom_line(aes(x = xi, y = truth, group = type), data = subset(est, metric == "est")) +
#     facet_wrap(metric ~ var, scales = "free") +
#     ylab("Metric Value") +
#     theme_bw() +
#     scale_color_brewer(palette="Dark2")
#
#   g
# }
#
#
# s1 <- subset(est1, scenario == 1) %>% g_summary2()
# s2 <- subset(est1, scenario == 2) %>% g_summary2()
# s3 <- subset(est1, scenario == 3) %>% g_summary2()
# s4 <- subset(est1, scenario == 4) %>% g_summary2()
#
# db <- simu_data(400)
# save(result_meta, est0, est1, db, g1, g2, g3, g4, s1, s2, s3, s4, file = "simulation/simu_qte_v5.Rdata")
