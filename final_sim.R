source("fun_mymodel.R")
size <- 10000

ts_f12_lambda <- 0.005
ts_f12_p <- 1.1
ts_f12_beta3 <- 0.004
ts_follow12 <- function(t) ts_f12_lambda * ts_f12_p * t ^ ( ts_f12_p - 1 )
plot(ts_follow12, 0, 50)

ts_a12_lambda <- 2e-5
ts_a12_p <- 3
ts_a12_beta3 <- 0.003
ts_age12 <- function(t) ts_a12_lambda * ts_a12_p * t ^ ( ts_a12_p - 1 )
plot(ts_age12, 0, 100)

ts_f13_lambda <- 0.01
ts_f13_p <- 0.9
ts_f13_beta1 <- 0.002
ts_f13_beta2 <- 0.003
ts_f13_beta3 <- 0.004
ts_follow13 <- function(t) ifelse(t > .Machine$double.eps, ts_f13_lambda * ts_f13_p * t ^ ( ts_f13_p - 1 ), 0) 
plot(ts_follow13, 0, 50)

ts_a13_lambda <- 1e-9
ts_a13_p <- 5
ts_a13_beta1 <- 0.001
ts_a13_beta2 <- 0.002
ts_a13_beta3 <- 0.003
ts_age13 <- function(t) ts_a13_lambda * ts_a13_p * t ^ ( ts_a13_p - 1 )
plot(ts_age13, 0, 100)

ts_d_lambda <- 0.03
ts_d_p <- 1.3
ts_d_beta1 <- 0.003
ts_d_beta2 <- 0.002
ts_d_beta3 <- 0.001
ts_disease <- function(d) ts_d_lambda * ts_d_p * d ^ ( ts_d_p - 1 )
plot(ts_disease, 0, 100)

ts_f23_lambda <- ts_f13_lambda
ts_f23_p <- ts_f13_p
ts_f23_beta1 <- ts_f13_beta1
ts_f23_beta2 <- ts_f13_beta2
ts_f23_beta3 <- ts_f13_beta3
ts_follow23 <- function(t) ts_f23_lambda * ts_f23_p * t ^ ( ts_f23_p - 1 )

ts_a23_lambda <- ts_a13_lambda
ts_a23_p <- ts_a13_p
ts_a23_beta1 <- ts_a13_beta1
ts_a23_beta2 <- ts_a13_beta2
ts_a23_beta3 <- ts_a13_beta3
ts_age23 <- function(t) ts_a23_lambda * ts_a23_p * t ^ ( ts_a23_p - 1 )

age <- runif(size, 20, 80)
x1 <- runif(size, 10, 20)
x2 <- runif(size, 10, 20)
x3 <- runif(size, 10, 20)

T1 <- mapply(
  function(a, c1, c2) random_survival_times(\(t) ts_follow12(t) * c1 + ts_age12(t + a) * c2, 1, 1000),
  age,
  exp(ts_f12_beta3 * x3),
  exp(ts_a12_beta3 * x3)
)

T1 <- ifelse(T1 <= .Machine$double.eps, .Machine$double.eps, T1)

T2 <- mapply(
  function(a, c1, c2) random_survival_times(\(t) ts_follow13(t) * c1 + ts_age13(t + a) * c2, 1, 1000),
  age,
  exp(ts_f13_beta1 * x1 + ts_f13_beta2 * x2 + ts_f13_beta3 * x3),
  exp(ts_a13_beta1 * x1 + ts_a13_beta2 * x2 + ts_a13_beta3 * x3)
)

for(i in which(T2 > T1)) {
  T2[[i]] <- random_survival_times(
    \(d) ts_disease(d) * exp(ts_d_beta1 * x1[[i]] + ts_d_beta2 * x2[[i]] + ts_d_beta3 * x3[[i]]) + 
      ts_follow23(d + T1[[i]]) * exp(ts_f23_beta1 * x1[[i]] + ts_f23_beta2 * x2[[i]] + ts_f23_beta3 * x3[[i]]) +
      ts_age23(d + age[[i]] + T1[[i]]) * exp(ts_a23_beta1 * x1[[i]] + ts_a23_beta2 * x2[[i]] + ts_a23_beta3 * x3[[i]]), 
    1
  ) + T1[[i]]
}

T2 <- ifelse(T2 <= .Machine$double.eps, .Machine$double.eps, T2)

plot(T1, T2)

data <- data.frame(
  x1 = x1,
  x2 = x2, 
  x3 = x3,
  age = age,
  time_to_disease = T1,
  time_to_death = T2
)

data$time_to_disease[data$ time_to_death <= data$time_to_disease] <- NA
ct <- runif(size, 0, 60)
data$status <- 1
data$status[is.na(data$time_to_disease) & ct < data$time_to_death] <- 0
data$time_to_death[is.na(data$time_to_disease) & ct < data$time_to_death] <- ct[is.na(data$time_to_disease) & ct < data$time_to_death]
data$status[(!is.na(data$time_to_disease)) & ct < data$time_to_disease] <- 0
data$time_to_death[(!is.na(data$time_to_disease)) & ct < data$time_to_disease] <- ct[(!is.na(data$time_to_disease)) & ct < data$time_to_disease]
data$time_to_disease[(!is.na(data$time_to_disease)) & ct < data$time_to_disease] <- NA
data$status[(!is.na(data$time_to_disease)) & data$time_to_disease <= ct & ct < data$time_to_death] <- 0
data$time_to_death[(!is.na(data$time_to_disease)) & data$time_to_disease <= ct & ct < data$time_to_death] <- ct[(!is.na(data$time_to_disease)) & data$time_to_disease <= ct & ct < data$time_to_death]

data2 <- data
data2$time_to_death[!is.na(data2$time_to_disease)] <- data2$time_to_disease[!is.na(data2$time_to_disease)]
data2$status[!is.na(data2$time_to_disease)] <- 0
data2$time_to_disease <- NULL


###
data2 <- data
data2$status <- as.numeric(data2$time_to_death <= data2$time_to_disease)
data2$time_to_death[data2$status == 0] <- data2$time_to_disease[data2$status == 0]

##
t <- proc.time()
mymodel(cbind(time_to_disease, 1) ~ I(x3 - mean(x3)), time.scales = "age", data = data2)
proc.time() - t