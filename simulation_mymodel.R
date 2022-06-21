source("fun_prepare.R")
library(optimx)
library(dplyr)
library(tibble)

digit <- 10
`%==%` <- function(a, b) near(a, b, 10^-digit)
`%^%` <- function(a, b) exp(b*log(a))
`%/%` <- function(a, b) ifelse(is.nan(a/b), 0, a/b)

lambda <- 1e-9
shape <- 5
lambda2 <- 0.01
shape2 <- 1.3
ts_age <- function(a) lambda*shape*a^(shape-1)
# plot(ts_age, 0, 100)
# random_survival_times(ts_age, 1000, 1000) %>% hist()
ts_follow <- function(t) lambda2*shape2*t^(shape2-1)
# plot(ts_follow, 0, 50)


size <- 10000
beta1 <- 0.01
beta2 <- 0.05
age <- runif(size, 40, 80)
x <- rnorm(size, 0, 3)


time <- mapply(
  function(a, c1, c2) random_survival_times(\(t) ts_follow(t)*c1 + ts_age(t + a)*c2, 1, 1000),
  age,
  exp(beta1 * x),
  exp(beta2 * x)
)

data <- data.frame(
  age = age,
  time = time,
  x = x,
  d_age = age + time
)

ct <- runif(size, 0, 60)
data$status <- as.numeric(data$time <= ct)
data$time[data$status == 0] <- ct[data$status == 0]
data$d_age <- data$age + data$time

ll <- function(par) {
  lambda <- exp(par[1])
  p1 <- exp(par[2])
  alpha <- exp(par[3])
  p2 <- exp(par[4])
  b1 <- par[5]
  b2 <- par[6]
  # lambda <- x[1]
  # p1 <- x[2]
  # alpha <- x[3]
  # p2 <- x[4]
  -with(data, sum(log(lambda * p1 * time[status == 1] ^ (p1-1) * exp(x[status == 1] * b1) + alpha * p2 * (d_age[status == 1]) ^ (p2-1) * exp(x[status == 1] * b2))) - lambda * sum(time ^ p1 * exp(x * b1)) - alpha * sum((d_age ^ p2 - age ^ p2) * exp(x * b2)))
}

gn <- function(par) {
  lambda <- exp(par[1])
  p1 <- exp(par[2])
  alpha <- exp(par[3])
  p2 <- exp(par[4])
  b1 <- par[5]
  b2 <- par[6]
  # lambda <- x[1]
  # p1 <- x[2]
  # alpha <- x[3]
  # p2 <- x[4]
  den <- with(data, lambda * p1 * time[status == 1] ^ (p1-1) * exp(x[status == 1] * b1) + alpha * p2 * d_age[status == 1] ^ (p2-1) * exp(x[status == 1] * b2))
  return(
    with(data,
      -c(
        (sum(p1 * exp(x[status == 1] * b1) * time[status == 1] ^ (p1-1) / den) - sum(time ^ p1 * exp(x * b1))) * lambda,
        (sum(lambda * exp(x[status == 1] * b1) * (time[status == 1] ^ (p1-1) + p1 * time[status == 1] ^ (p1-1) * log(time[status == 1])) / den) - lambda * sum(time ^ p1 * log(time) * exp(x * b1))) * p1,
        (sum(p2 * exp(x[status == 1] * b2) * (d_age[status == 1]) ^ (p2-1) / den) - sum((d_age ^ p2 - age ^ p2) * exp(x * b2))) * alpha,
        (sum(alpha * exp(x[status == 1] * b2) * (d_age[status == 1] ^ (p2-1) + p2 * d_age[status == 1] ^ (p2-1) * log(d_age[status == 1])) / den) - alpha * sum((d_age ^ p2 * log(d_age) - age ^ p2 * log(age)) * exp(x * b2))) * p2,
        sum(lambda * p1 * exp(x[status == 1] * b1) * time[status == 1] ^ (p1-1) * x[status == 1] / den) - lambda * sum(exp(x * b1) * time ^ p1 * x),
        sum(alpha * p2 * exp(x[status == 1] * b2) * d_age[status == 1] ^ (p2-1) * x[status == 1] / den) - alpha * sum(exp(x * b2) * (d_age ^ p2 - age ^ p2) * x)
      )
    )
  )
}

# hess <- function(par) {
#   lambda <- exp(par[1])
#   p1 <- exp(par[2])
#   alpha <- exp(par[3])
#   p2 <- exp(par[4])
#   b1 <- par[5]
#   b2 <- par[6]
#   den <- with(data, lambda * p1 * time[status == 1] ^ (p1-1) * exp(x[status == 1] * b1) + alpha * p2 * d_age[status == 1] ^ (p2-1) * exp(x[status == 1] * b2))
#   
#   return(
#     diag(
#       with(data,
#            -c(
#              (sum(-(p1 * time[status == 1] ^ (p1-1) * exp(x[status == 1] * b1))^2 / den^2)) * lambda^2 + (sum(p1 * exp(x[status == 1] * b1) * time[status == 1] ^ (p1-1) / den) - sum(time ^ p1 * exp(x * b1))) * lambda,
#              (sum((lambda * exp(x[status == 1] * b1) * time[status == 1] ^ (p1-1) * log(time[status == 1]) * (2 + p1 * log(time[status == 1])) * den - (lambda * exp(x[status == 1] * b1) * time[status == 1] ^ (p1-1) * (1 + p1 * log(time[status == 1])))^2) / den^2) - lambda * sum(exp(x * b1) * time ^ p1 * (log(time))^2)) * p1^2 + (sum(lambda * exp(x[status == 1] * b1) * (time[status == 1] ^ (p1-1) + p1 * time[status == 1] ^ (p1-1) * log(time[status == 1])) / den) - lambda * sum(time ^ p1 * log(time) * exp(x * b1))) * p1,
#              (sum(-(p2 * d_age[status == 1] ^ (p2-1) * exp(x[status == 1] * b2))^2 / den^2)) * alpha^2 + (sum(p2 * exp(x[status == 1] * b2) * d_age[status == 1] ^ (p2-1) / den) - sum((d_age ^ p2 - age ^ p2) * exp(x * b2))) * alpha,
#              (sum((alpha * exp(x[status == 1] * b2) * d_age[status == 1] ^ (p2-1) * log(d_age[status == 1]) * (2 + p2 * log(d_age[status == 1])) * den - (alpha * exp(x[status == 1] * b2) * d_age[status == 1] ^ (p2-1) * (1 + p2 * log(d_age[status == 1])))^2) / den^2) - alpha * sum(exp(x * b2) * (d_age ^ p2 * (log(d_age))^2 - age ^ p2 * (log(age))^2))) * p2^2 + (sum(alpha * exp(x[status == 1] * b2) * (d_age[status == 1] ^ (p2-1) + p2 * d_age[status == 1] ^ (p2-1) * log(d_age[status == 1])) / den) - alpha * sum((d_age ^ p2 * log(d_age) - age ^ p2 * log(age)) * exp(x * b2))) * p2,
#              sum(lambda * p1 * time[status == 1] ^ (p1-1) * exp(x[status == 1] * b1) * x[status == 1]^2 * alpha * p2 * d_age[status == 1] ^ (p2-1) * exp(x[status == 1] * b2) / den^2) - lambda * sum(time ^ p1 * exp(x * b1) * x^2),
#              sum(alpha * p2 * d_age[status == 1] ^ (p2-1) * exp(x[status == 1] * b2) * x[status == 1]^2 * lambda * p1 * time[status == 1] ^ (p1-1) * exp(x[status == 1] * b1) / den^2) - alpha * sum((d_age ^ p2 - age ^ p2) * exp(x * b2) * x^2)
#            )
#       )
#     )
#   )
# }


opt <- optimx(c(lambda = 0, p1 = 0, alpha = 0, p2 = 0, b1 = 0, b2 = 0), ll, gn, hessian = TRUE, control = list(trace = 0, all.methods = TRUE))

summary(opt, order = "value") %>%
  rownames_to_column("algorithm") %>%
  filter(convcode != 9999) %>%
  arrange(value) %>%
  select(algorithm, lambda, p1, alpha, p2, b1, b2, value, convcode) %>%
  mutate(lambda = exp(lambda), p1 = exp(p1), alpha = exp(alpha), p2 = exp(p2), b1 = b1, b2 = b2) %>%
  head(7)



f <- function(i) {
  tryCatch({size <- 1000
  beta1 <- 0.1
  beta2 <- 0.2
  age <- runif(size, 40, 80)
  x <- rnorm(size, 0, 3)
  
  
  time <- mapply(
    function(a, c1, c2) random_survival_times(\(t) ts_follow(t)*c1 + ts_age(t + a)*c2, 1, 1000),
    age,
    exp(beta1 * x),
    exp(beta2 * x)
  )
  
  data <- data.frame(
    age = age,
    time = time,
    x = x,
    d_age = age + time
  )
  
  ct <- runif(size, 0, 60)
  data$status <- as.numeric(data$time <= ct)
  data$time[data$status == 0] <- ct[data$status == 0]
  data$d_age <- data$age + data$time
  
  ll <- function(par) {
    lambda <- exp(par[1])
    p1 <- exp(par[2])
    alpha <- exp(par[3])
    p2 <- exp(par[4])
    b1 <- par[5]
    b2 <- par[6]
    # lambda <- x[1]
    # p1 <- x[2]
    # alpha <- x[3]
    # p2 <- x[4]
    -with(data, sum(log(lambda * p1 * time[status == 1] ^ (p1-1) * exp(x[status == 1] * b1) + alpha * p2 * (d_age[status == 1]) ^ (p2-1) * exp(x[status == 1] * b2))) - lambda * sum(time ^ p1 * exp(x * b1)) - alpha * sum((d_age ^ p2 - age ^ p2) * exp(x * b2)))
  }
  
  gn <- function(par) {
    lambda <- exp(par[1])
    p1 <- exp(par[2])
    alpha <- exp(par[3])
    p2 <- exp(par[4])
    b1 <- par[5]
    b2 <- par[6]
    # lambda <- x[1]
    # p1 <- x[2]
    # alpha <- x[3]
    # p2 <- x[4]
    den <- with(data, lambda * p1 * time[status == 1] ^ (p1-1) * exp(x[status == 1] * b1) + alpha * p2 * d_age[status == 1] ^ (p2-1) * exp(x[status == 1] * b2))
    return(
      with(data,
           -c(
             (sum(p1 * exp(x[status == 1] * b1) * time[status == 1] ^ (p1-1) / den) - sum(time ^ p1 * exp(x * b1))) * lambda,
             (sum(lambda * exp(x[status == 1] * b1) * (time[status == 1] ^ (p1-1) + p1 * time[status == 1] ^ (p1-1) * log(time[status == 1])) / den) - lambda * sum(time ^ p1 * log(time) * exp(x * b1))) * p1,
             (sum(p2 * exp(x[status == 1] * b2) * (d_age[status == 1]) ^ (p2-1) / den) - sum((d_age ^ p2 - age ^ p2) * exp(x * b2))) * alpha,
             (sum(alpha * exp(x[status == 1] * b2) * (d_age[status == 1] ^ (p2-1) + p2 * d_age[status == 1] ^ (p2-1) * log(d_age[status == 1])) / den) - alpha * sum((d_age ^ p2 * log(d_age) - age ^ p2 * log(age)) * exp(x * b2))) * p2,
             sum(lambda * p1 * exp(x[status == 1] * b1) * time[status == 1] ^ (p1-1) * x[status == 1] / den) - lambda * sum(exp(x * b1) * time ^ p1 * x),
             sum(alpha * p2 * exp(x[status == 1] * b2) * d_age[status == 1] ^ (p2-1) * x[status == 1] / den) - alpha * sum(exp(x * b2) * (d_age ^ p2 - age ^ p2) * x)
           )
      )
    )
  }
  
  # library(optimx)
  opt <- optimx(c(0,0,0,0,0,0), ll, gn, hessian = TRUE, method = "nlminb", control = list(trace = 0))
  # library(dplyr)
  # library(tibble)
  # summary(opt, order = "value") %>%
  #   rownames_to_column("algorithm") %>%
  #   filter(convcode != 9999) %>%
  #   arrange(value) %>%
  #   select(algorithm, p1, p2, p3, p4, value) %>%
  #   mutate(p1 = exp(p1), p2 = exp(p2), p3 = exp(p3), p4 = exp(p4)) %>%
  #   head(7)
  if (opt$convcode != 0) return(rep(NaN, 6))
  return(exp(with(opt, c(p1,p2,p3,p4,p5,p6))))}, error = function(e) return(rep(NaN, 6)))
}
tmp <- sapply(1:100, f)
apply(tmp, 1, \(x) mean(x, na.rm = TRUE))
