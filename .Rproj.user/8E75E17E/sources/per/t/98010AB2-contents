library(survival)
library(ggplot2)
library(magrittr)
library(timereg)

theme_set(theme(panel.background = element_rect(fill = NA, colour = "black"),
                panel.grid = element_blank()))

inverse <- function(fn, min_x, max_x){
  fn_inv <- function(y){
    uniroot((function(x){fn(x) - y}), lower=min_x, upper=max_x)[1]$root
  }
  return(Vectorize(fn_inv))
}

integrate_from_0 <- function(fn, t){
  int_fn <- function(t) integrate(fn, 0, t)
  result <- sapply(t, int_fn)
  value  <- unlist(result["value",])
  msg    <- unlist(result["message",])
  value[which(msg != "OK")] = NA
  return(value)
}

random_survival_times <- function(hazard_fn, n, max_time=100){
  cumulative_density_fn <- function(t) 1 - exp(-integrate_from_0(hazard_fn, t))
  inverse_cumulative_density_fn <- inverse(cumulative_density_fn, 0, max_time)
  return(inverse_cumulative_density_fn(runif(n)))
}

ramp <- function(x) (abs(x)+x)/2

near <- function (x, y, tol = .Machine$double.eps){
  abs(x - y) < tol
}
