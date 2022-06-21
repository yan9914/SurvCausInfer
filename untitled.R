source("fun_prepare.R")

lambda <- 2
shape <- 1.5
hazard_fn <- function(t) lambda*shape*t^(shape-1)
#hist(random_survival_times(hazard_fn, 200))

size <- 2000L
beta1 <- 2
T1 <- rexp(size, 2.5)
T2 <- random_survival_times(hazard_fn, size)
for(i in 1:size) {
  tmp <- T2[[i]]
  if(tmp >= T1[[i]]) {
    tmp <- random_survival_times(
      function(t) hazard_fn(t + T1[[i]])*exp(beta1*t), 
      1
    ) + T1[[i]]
    while(tmp < T1[[i]]) {
      tmp <- random_survival_times(
        function(t) hazard_fn(t + T1[[i]])*exp(beta1*t), 
        1
      ) + T1[[i]]
    }
    T2[[i]] <- tmp
  }
  #print(i)
}

df <- data.frame(
  id = 1:size,
  t1 = T1,
  t2 = T2,
  status = 1
)

fit <- coxph(Surv(t2, status) ~ tt(t1), data = df, tt = function(t1, t, ...) ifelse(t > t1, t-t1, 0))
summary(fit)
ggplot(df, aes(x = t1, y = t2)) + 
  stat_density2d(aes(color = ..level..)) + 
  stat_function(fun = function(x) x, xlim = c(0,1), linetype = "dashed", color = "red")

###########################
size <- 5000L
b0 <- 2
b1 <- 4
T1 <- rexp(size, 2.5)
T2 <- random_survival_times(hazard_fn, size)
for(i in 1:size) {
  tmp <- T2[[i]]
  if(tmp >= T1[[i]]) {
    tmp <- random_survival_times(
      function(t) hazard_fn(t + T1[[i]])*exp(b0/b1*(1-exp(-b1*t))), 
      1
    ) + T1[[i]]
    while(tmp < T1[[i]]) {
      tmp <- random_survival_times(
        function(t) hazard_fn(t + T1[[i]])*exp(b0/b1*(1-exp(-b1*t))), 
        1
      ) + T1[[i]]
    }
    T2[[i]] <- tmp
  }
  #print(i)
}

df <- data.frame(
  id = 1:size,
  t1 = T1,
  t2 = T2,
  status = 1
)

prr <- function(t1, t2, status) {
  st2 <- sort(t2)
  st1 <- t1[order(t2)]
  ss <- status[order(t2)]
  ft <- sort(t2[status==1])
  n <- length(t1)
  nf <- length(ft)
  LL <- function(b0, b1) {
    num <- numeric(n)
    den <- numeric(n)
    for(i in 1:n) {
      den <- den + c(exp(b0/b1*(1-exp(-b1*ramp(head(ft, i)-st1[i])))), rep(0, n-i))
      num[[i]] <- ifelse(
        ss[[i]] == 0,
        0,
        exp(b0/b1*(1-exp(-b1*ramp(ft[i]-st1[i]))))
      )
    }
    return(sum(log(num/den)))
  }
  output <- optim(c(0,0), fn = function(x) -LL(x[1],exp(x[2])))
  output$par[[2]] <- exp(output$par[[2]])
  return(output)
}
prr(df$t1,df$t2,df$status)

###########################
size <- 1000L
b0 <- 10
b1 <- 8
T1 <- rexp(size, 2.5)
T2 <- random_survival_times(hazard_fn, size)
for(i in 1:size) {
  tmp <- T2[[i]]
  if(tmp >= T1[[i]]) {
    tmp <- random_survival_times(
      function(t) hazard_fn(t + T1[[i]]) + exp(b0/b1*(1-exp(-b1*t))), 
      1
    ) + T1[[i]]
    while(tmp < T1[[i]]) {
      tmp <- random_survival_times(
        function(t) hazard_fn(t + T1[[i]]) + exp(b0/b1*(1-exp(-b1*t))), 
        1
      ) + T1[[i]]
    }
    T2[[i]] <- tmp
  }
  #print(i)
}

df <- data.frame(
  id = 1:size,
  t1 = T1,
  t2 = T2,
  status = 1
)

sm <- 20L
sojourn <- (df$t2 - df$t1) %>% {.[.>0]} %>%
  quantile(seq(0, 1, len = sm+2)) %>%
  {.[-c(1,length(.))]}

df2 <- data.frame(matrix(nrow = 0, ncol = 4 + sm))
for(i in 1:max(df$id)) {
  t1 <- df[df$id == i,]$t1
  t2 <- df[df$id == i,]$t2
  if(t2 <= t1) {
    df2[nrow(df2)+1,] <- c(i, 0, t2, 1, rep(0, sm))
  } else {
    df2[nrow(df2)+1,] <- c(i, 0, t1, 0, rep(0, sm))
    suppressWarnings({tmp <- max(which(sojourn <= t2 - t1))})
    if(tmp == -Inf) {
      df2[nrow(df2)+1,] <- c(i, t1, t2, 1, rep(0, sm))
    } else {
      df2[nrow(df2) + 1:tmp,] <- cbind(rep(i, tmp),
                                       t1 + c(0, head(sojourn, tmp-1)),
                                       t1 + sojourn[1:tmp],
                                       rep(0, tmp),
                                       diag(nrow = tmp, ncol = sm)
      )
      if(near(sojourn[tmp], t2 - t1)) {
        df2[nrow(df2),4] <- 1
      } else {
        df2[nrow(df2) + 1,] <- c(i, t1 + sojourn[tmp], t2, 1, c(rep(0, tmp-1), 1, rep(0, sm-tmp)))
      }
    }
  }
}
fit <- aalen(as.formula(paste("Surv(X2, X3, X4)~",paste0("const(X", 5:(4+sm), ")", collapse = " + "))),
      df2,
      id = df2$X1
)
summary(fit)
plot(fit)
ggplot() +
  geom_point(aes(sojourn, coef(fit)[,1])) +
  stat_function(fun = function(x) exp(b0/b1*(1-exp(-b1*x)))) +
  geom_line(aes(sojourn, mono_fit(sojourn, coef(fit)[,1], "increasing")$fitted), color = "red")

ggplot(as.data.frame(fit$cum)) +
  geom_step(aes(time, `(Intercept)`)) +
  stat_function(fun = function(t) integrate_from_0(hazard_fn, t))

pad <- function(t1, t2, status, sm) {
  st2 <- sort(t2)
  st1 <- t1[order(t2)]
  ss <- status[order(t2)]
  ft <- sort(unique(t2[status==1]))
  n <- length(t1)
  nf <- length(ft)
  dft <- diff(c(0, ft))
  sojourn <- (df$t2 - df$t1) %>% {.[.>0]} %>%
    quantile(seq(0, 1, len = sm+2)) %>%
    {.[-c(1,length(.))]}
  base <- numeric(nf)
  for(i in 1:nf) {
    j <- st2 >= ft[[i]]
    base[[i]] <- mean(st2[j] == ft[[i]])
  }
  res <- function(base, b0, b1) {
    tmp <- 0
    for(i in 1:nf) {
      j <- st2 >= ft[[i]]
      tmp <- tmp + sum(((st2[j] == ft[[i]]) - base[[i]] - exp(b0/b1*(1-exp(-b1*ramp(ft[[i]]-st1[j])))))^2)
    }
    return(tmp)
  }
  b0_temp <- 0
  lb1_temp <- 0
  output <- optim(c(b0_temp,lb1_temp), fn = function(x) res(base, x[1], exp(x[2])))
  b0_temp2 <- output$par[[1]]
  lb1_temp2 <- output$par[[2]]
  while(max(abs(c(b0_temp-b0_temp2, lb1_temp-lb1_temp2))) > 1e-8) {print("HI")
    b0_temp <- b0_temp2
    lb1_temp <- lb1_temp2
    for(i in 1:nf) {
      j <- st2 >= ft[[i]]
      base[[i]] <- mean((st2[j] == ft[[i]]) - dft[[i]] * exp(b0_temp/exp(lb1_temp)*(1-exp(-exp(lb1_temp)*ramp(ft[[i]]-st1[j])))))
    }
    output <- optim(c(b0_temp,lb1_temp), fn = function(x) res(base, x[1], exp(x[2])))
    b0_temp2 <- output$par[[1]]
    lb1_temp2 <- output$par[[2]]
  }
  return(list(time = ft, base = base, c(b0_temp2, exp(lb1_temp2))))
}
fit <- pad(df$t1,df$t2,df$status)
ggplot() + geom_step(aes(fit$time, cumsum(fit$base))) +
  stat_function(fun = function(t) integrate_from_0(hazard_fn, t))
ggplot() + stat_function(fun = function(t) exp(fit[[3]][1]/fit[[3]][2]*(1-exp(-fit[[3]][2]*t))))
  
