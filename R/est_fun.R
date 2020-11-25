# based on observed data only (remove death as well)
est_fun1 <- function(q, xi, db){
  y_std <- (q - db$y_mean) / db$y_sigma
  f_q <- pnorm(y_std)
  term1 <- f_q - xi
  term2 <- ( (db$y_obs < q) - f_q ) / db$score
  mean( (term1 + term2) / db$km, na.rm = TRUE )
}
