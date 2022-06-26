source("fun_prepare.R")
library(optimx)
library(dplyr)
library(tibble)

# lambda <- 1e-9
# shape <- 5
# lambda2 <- 0.01
# shape2 <- 1.3
# ts_age <- function(a) lambda*shape*a^(shape-1)
# # plot(ts_age, 0, 100)
# # random_survival_times(ts_age, 1000, 1000) %>% hist()
# ts_follow <- function(t) lambda2*shape2*t^(shape2-1)
# # plot(ts_follow, 0, 50)
# 
# 
# size <- 10000
# beta1 <- 0.01
# beta2 <- 0.05
# age <- runif(size, 40, 80)
# x <- rnorm(size, 0, 3)
# x <- runif(size, 0, 2)
# 
# time <- mapply(
#   function(a, c1, c2) random_survival_times(\(t) ts_follow(t)*c1 + ts_age(t + a)*c2, 1, 1000),
#   age,
#   exp(beta1 * x),
#   exp(beta2 * x)
# )
# 
# data <- data.frame(
#   age = age,
#   time = ifelse(time <= .Machine$double.eps, .Machine$double.eps, time),
#   x = x,
#   d_age = age + time
# )
# 
# ct <- runif(size, 0, 60)
# data$status <- as.numeric(data$time <= ct)
# data$time[data$status == 0] <- ct[data$status == 0]
# data$d_age <- data$age + data$time



### customize function

mymodel <- function(formula, data, time.scales, opt.method = NULL) {
  Call <- match.call()
  indx <- match(c("formula", "data"), names(Call), 0L)
  temp <- Call[c(1L, indx)]
  temp$drop.unused.levels <- TRUE
  temp[[1L]] <- quote(stats::model.frame)
  temp$formula <- terms(formula, data = data)
  m <- eval(temp, parent.frame())
  Terms <- attr(m, "terms")
  Y <- model.extract(m, "response")
  X <- model.matrix(Terms, m) %>% data.frame() %>% select(-1)
  ts <- data[,time.scales] %>% as.data.frame()
  colnames(ts) <- sapply(seq_along(time.scales), \(i) paste0("ts.", time.scales[[i]]))
  n.ts <- ncol(ts) + 1
  n.x <- ncol(X)
  n.par <- 2 * n.ts + n.x * n.ts
  ll <- function(par, n.x. = n.x, n.ts. = n.ts, Y. = Y, X. = X, ts. = ts) {
    if(any(is.na(par))) stop("NA/NaN parameters")
    new_par <- par
    new_par[c(seq(1, by = 2+n.x., len = n.ts.), seq(2, by = 2+n.x., len = n.ts.))] <- exp(par[c(seq(1, by = 2+n.x., len = n.ts.), seq(2, by = 2+n.x., len = n.ts.))]) # weibull parameter
    tmp1 <- new_par[1] * new_par[2] * Y.[Y.[,2] == 1,1] ^ (new_par[2] - 1) * exp(t(t(X.[Y.[,2] == 1,])) %*% new_par[2+1:n.x.])
    
    tmp2 <- -new_par[1] * sum(Y.[,1] ^ new_par[2] * exp(t(t(X.)) %*% new_par[2+1:n.x.]))
    for(i in 1:(n.ts.-1)) {
      tmp1 <- tmp1 + new_par[1 + (2+n.x.)*i] * new_par[2 + (2+n.x.)*i] * (Y.[Y.[,2] == 1,1] + ts.[Y.[,2] == 1,i]) ^ (new_par[2 + (2+n.x.)*i] - 1) * exp(t(t(X.[Y.[,2] == 1,])) %*% new_par[2+1:n.x.+(2+n.x.)*i])
      tmp2 <- tmp2 - new_par[1 + (2+n.x.)*i] * sum(((Y.[,1] + ts.[,i])^new_par[2 + (2+n.x.)*i]-ts.[,i]^new_par[2 + (2+n.x.)*i]) * exp(t(t(X.)) %*% new_par[2+1:n.x.+(2+n.x.)*i]))
    }
    result <- -sum(log(tmp1))-tmp2
    # print(list(par = par, value = result))
    if(is.na(result)) stop(paste0("NA result for parameter: ", list(par)))
    return(result)
  }
  par <- rep(0, n.par)
  names(par) <- c(c("ts.prime.lambda", "ts.prime.p", sapply(1:n.x, \(j) paste0("ts.prime.", colnames(X)[j]))), sapply(2:n.ts, \(i) c(paste0(colnames(ts)[i-1], ".", "lambda"), paste0(colnames(ts)[i-1], ".", "p"), sapply(1:n.x, \(j) paste0(colnames(ts)[i-1], ".", colnames(X)[j])))))
  if(is.null(opt.method)) {
    opt <- optimx(par, ll, hessian = TRUE, control = list(trace = 0, all.methods = TRUE))
    
    summary(opt, order = "value") %>%
      rownames_to_column("algorithm") %>%
      filter(convcode != 9999) %>%
      arrange(value) %>%
      select(algorithm, names(par), value, convcode) %>%
      mutate(across(ends_with(c(".lambda", ".p")), exp)) %>%
      head(7)
  } else {
    if(!any(opt.method %in% c('Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'nlm', 'nlminb', 'spg', 'ucminf', 'newuoa', 'bobyqa', 'nmkb', 'hjkb', 'Rcgmin', 'Rvmmin'))) {
      stop("Not available method. Possible method codes at the time of writing are 'Nelder-Mead', 'BFGS', 'CG', 'L-BFGS-B', 'nlm', 'nlminb', 'spg', 'ucminf', 'newuoa', 'bobyqa', 'nmkb', 'hjkb', 'Rcgmin', or 'Rvmmin'.")
    }
    opt <- optimx(par, ll, method = opt.method, hessian = TRUE, control = list(trace = 0))
    summary(opt, order = "value") %>%
      rownames_to_column("algorithm") %>%
      arrange(value) %>%
      select(algorithm, names(par), value, convcode) %>%
      mutate(across(ends_with(c(".lambda", ".p")), exp))
  }
}

# mymodel(cbind(time, status) ~ x, time.scales = "age", data = data)

