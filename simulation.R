source("fun_prepare.R")
library(Epi)
library(popEpi)
library(latex2exp)

lambda0 <- 1.3
lambda1 <- 2
lambda2 <- 1.5
shape <- 1.2
hazard_fn1 <- function(t) lambda1*shape*t^(shape-1)
hazard_fn2 <- function(t) lambda2*shape*t^(shape-1)
#hist(random_survival_times(hazard_fn, 200))

size <- 2000L
# phi <- function(d) 0.1 * d
# eta <- function(r) 3 + 2 * r
phi <- function(d) 0
eta <- function(r) 1 * r

# fra_var <- 1/2
# frailty <- rgamma(size, 1/fra_var, rate = 1/fra_var)
frailty <- rep(1, size)

X1 <- rnorm(size)
X2 <- runif(size, 0, 2)
coef12 <- c(1, 2)
coef13 <- c(2.5, 1)
coef23 <- c(2, 1.5)



X <- cbind(X1, X2)
T1 <- rexp(size, lambda0 * exp(X %*% coef12) * frailty)
T2 <- sapply(
  exp(X %*% coef13) * frailty,
  \(x) random_survival_times(\(t) hazard_fn1(t) * x, 1, 1000)
)

for(i in which(T2 >= T1)) {
  T2[[i]] <- random_survival_times(
    \(t) hazard_fn2(t + T1[[i]])*exp(c(X[i,]%*%coef23) + phi(t) + eta(T1[[i]])) * frailty[[i]], 
    1
  ) + T1[[i]]
}

df <- data.frame(
  id = 1:size,
  t1 = T1,
  t2 = T2,
  X1 = X1,
  X2 = X2,
  status = 1
)



# fit <- coxph(Surv(t2, status) ~ tt(t1) + X1, data = df, tt = function(t1, t, ...) (t > t1))
# summary(fit)
# basehaz(fit)
ggplot(df, aes(x = t1, y = t2)) +
  geom_point() +
  geom_abline(slope = 1, linetype = "dashed", color = "red")


lex <- Lexis(
  entry = list(TSE = rep(0,size)), # TSE: Time Since Entry
  exit = list(TSE = t2),
  entry.status = factor(rep(0, size), label = "State1"),
  exit.status = factor(rep(0, size), label = "State3"),
  data = df
)
lex2 <- cutLexis(
  lex, 
  cut = lex$t1,
  pre = "State1",
  new.state = "State2",
  new.scale = "TSI", # TSI: Time Since Intermediate Event
  split.states = TRUE
)
summary(lex2, t = F)
boxes(lex2, boxpos = TRUE, scale.R = 1000, hmult = 1.5)
lex3 <- splitMulti(lex2, TSE = seq(0, max(lex2$TSE), length = 500))
lex3$TSI <- ifelse(is.na(lex3$TSI), 0, lex3$TSI)
lex3$TTI <- ifelse(lex3$lex.Cst == "State2", lex3$t1, 0)
nk <- 5
( a.kn <- c(0, with( subset(lex3, lex.Xst == "State3(State2)"),
                     quantile(TSE + lex.dur, (1:nk-0.5)/nk) ) ) )
( d.kn <- c(0, with( subset(lex3, lex.Cst == "State2" & lex.Xst == "State3(State2)"),
                     quantile(TTI, (1:nk-0.5)/nk) ) ) )
( sd.kn <- c(0, with( subset(lex3, lex.Cst == "State2" & lex.Xst=="State3(State2)"),
                      quantile(TSI, (1:nk-0.5)/nk) ) ) )
pm23 <- glm( cbind(lex.Xst == "State3(State2)", lex.dur) ~ Ns(TSE, knots = a.kn)
             # + Ns(TTI, knots = d.kn)
             + Ns(TSI, knots = sd.kn)
             + X1 + X2,
             # + lex.Cst,
           # + lex.Cst * X1
           # + lex.Cst * X2,
           family = poisreg, data = lex3[lex.Cst == "State2",] )
summary(pm23)

ggplot() +
  geom_line(aes(x = seq(0, a.kn[nk+1], len = 500),
                y = predict(pm23, 
                            data.frame(
                              TSE = seq(0, a.kn[nk+1], len = 500),
                              TTI = 0,
                              TSI = 0,
                              X1 = 0,
                              X2 = 0
                            )) %>% exp)) +
  stat_function(fun = hazard_fn2, col = "red")

ggplot() +
  geom_line(aes(x = seq(0, d.kn[nk+1], len = 500),
                y = predict(pm23, 
                            data.frame(
                              TSE = 0,
                              TTI = 0,
                              TSI = seq(0, sd.kn[nk+1], len = 500),
                              X1 = 0,
                              X2 = 0
                            ))))


pm13 <- glm( cbind(lex.Xst == "State3", lex.dur) ~ Ns(TSE, knots = a.kn)
             + X1 + X2,
             family = poisreg, data = lex3[lex3$lex.Cst=="State1",] )
summary(pm13)
plot(\(x) predict(pm13,
                  data.frame(
                    TSE = x,
                    X1 = 0, X2 = 0
                  )) %>% exp, 0, a.kn[nk+1])
par(new = T)
plot(hazard_fn1, 0, a.kn[nk+1], col = "red")

# ggplot() +
#   geom_line(aes(x = seq(0, quantile(lex3$Age, 0.9), 0.01), y = exp(predict(pm13, data.frame(Age = seq(0, quantile(lex3$Age, 0.9), 0.01), lex.Cst = "State1", TTI = 0, TSI = 0, X1 = 0, X2 = 0))))) +
#   stat_function(fun = hazard_fn, color = "red")
# ggplot() +
#   geom_line(aes(x = seq(0, quantile(lex3$TTI, 0.9), 0.01), y = exp(predict(pm13, data.frame(Age = 0, TTI =  seq(0, quantile(lex3$TTI, 0.9), 0.01), lex.Cst = "State2", TSI = 0, X1 = 0, X2 = 0))) %>% log)) +
#   stat_function(fun = \(t) 0.3*t +0.5, color = "red")

h1 <- function(time, t1, x1, x2) {
  pm23 %>%
    predict(
      data.frame(
        TSE = time,
        TTI = t1,
        TSI = time - t1,
        lex.Cst = "State2",
        X1 = x1, X2 = x2
      )
    ) %>%
    exp
}
h0 <- function(time, x1, x2) {
  pm13 %>%
    predict(
      data.frame(
        TSE = time,
        X1 = x1, X2 = x2
      )
    ) %>%
    exp
} 

D <- function(time, t1) {
  c1 <- integrate(\(t) h0(t,0,0), 0, t1)$value
  c2 <- integrate(\(t) h1(t,t1,0,0), t1, time)$value
  h1(time,t1,0,0)*mean(exp(pm23$coefficients["X1"]*df$X1 + pm23$coefficients["X2"]*df$X2 - c1 * exp(pm13$coefficients["X1"]*df$X1 + pm13$coefficients["X2"]*df$X2) - c2 * exp(pm23$coefficients["X1"]*df$X1 + pm23$coefficients["X2"]*df$X2))) /
    (h0(time,0,0)*mean(exp(pm13$coefficients["X1"]*df$X1 + pm13$coefficients["X2"]*df$X2 - c1 * exp(pm13$coefficients["X1"]*df$X1 + pm13$coefficients["X2"]*df$X2) - c2 * exp(pm23$coefficients["X1"]*df$X1 + pm23$coefficients["X2"]*df$X2))))
}
D <- Vectorize(D, "time")

D2 <- function(time, t1) {
  h1(time,t1,0,0)*mean(exp(pm23$coefficients["X1"]*df$X1 + pm23$coefficients["X2"]*df$X2)) /
    (h0(time,0,0)*mean(exp(pm13$coefficients["X1"]*df$X1 + pm13$coefficients["X2"]*df$X2)))
}
D2 <- Vectorize(D2, "time")

#D(1, 0.3)
#exp(predict(pm13, data.frame(Age=0, TTI = 0.3, TSI = 0.7, lex.Cst = "State2", X1 = 0, X2 = 0))-predict(pm13, data.frame(Age=0, TTI = 0, TSI =0, lex.Cst = "State1", X1 = 0, X2 = 0)))



tmp <- seq(0.3, quantile(lex3$TSE, 0.9), length = 50)

a <- proc.time()
dt1 <- D(tmp, 0.3)
dt2 <- D2(tmp, 0.3)
proc.time()-a

ggplot() +
  geom_line(aes(tmp, dt1)) +
  geom_line(aes(tmp, dt2), color = "red")


write.csv(data.frame(tmp, dt1, dt2), file = "theta3o2.csv")

completecsv <- read.csv("complete.csv")
ggplot(data = completecsv) +
  geom_line(aes(tmp, dt1)) +
  geom_line(aes(tmp, dt2), color = "red") +
  xlab("t") +
  ylab(TeX("$\\Delta D_1(t,0.3)$")) +
  ylim(c(0.5,6))
theta3o2csv <- read.csv("theta3o2.csv")
ggplot(data = theta3o2csv) +
  geom_line(aes(tmp, dt1)) +
  geom_line(aes(tmp, dt2), color = "red") +
  xlab("t") +
  ylab(TeX("$\\Delta D_1(t,0.3)$")) +
  ylim(c(0.5,6))
theta2csv <- read.csv("theta2.csv")
ggplot(data = theta2csv) +
  geom_line(aes(tmp, dt1)) +
  geom_line(aes(tmp, dt2), color = "red") +
  xlab("t") +
  ylab(TeX("$\\Delta D_1(t,0.3)$")) +
  ylim(c(0.5,6))
theta3csv <- read.csv("theta3.csv")
ggplot(data = theta3csv) +
  geom_line(aes(tmp, dt1)) +
  geom_line(aes(tmp, dt2), color = "red") +
  xlab("t") +
  ylab(TeX("$\\Delta D_1(t,0.3)$")) +
  ylim(c(0.5,6))

x2is0d1csv <- read.csv("x2is0.1.csv")
ggplot(data = x2is0d1csv) +
  geom_line(aes(tmp, dt1)) +
  geom_line(aes(tmp, dt2), color = "red") +
  xlab("t") +
  ylab(TeX("$\\Delta D_1(t,0.3)$")) +
  ylim(c(1,3.5))
x2is1csv <- read.csv("x2is1.csv")
ggplot(data = x2is1csv) +
  geom_line(aes(tmp, dt1)) +
  geom_line(aes(tmp, dt2), color = "red") +
  xlab("t") +
  ylab(TeX("$\\Delta D_1(t,0.3)$")) +
  ylim(c(1,3.5))
x2is2csv <- read.csv("x2is2.csv")
ggplot(data = x2is2csv) +
  geom_line(aes(tmp, dt1)) +
  geom_line(aes(tmp, dt2), color = "red") +
  xlab("t") +
  ylab(TeX("$\\Delta D_1(t,0.3)$")) +
  ylim(c(1,3.5))
x2is3csv <- read.csv("x2is3.csv")
ggplot(data = x2is3csv) +
  geom_line(aes(tmp, dt1)) +
  geom_line(aes(tmp, dt2), color = "red") +
  xlab("t") +
  ylab(TeX("$\\Delta D_1(t,0.3)$")) +
  ylim(c(1,3.5))
# df2 <- df
# df2$status <- as.numeric(df2$t2 < df2$t1)
# df2$t22 <- ifelse(df2$t2 < df2$t1, df2$t2, df2$t1)
# pm23 <- glm( cbind(lex.Xst=="State3",lex.dur) ~ Ns(Age,knots=a.kn)
#            + X1 + X2,
#            family=poisreg, data = lex3[lex3$lex.Cst=="State1",] )
# fit2 <- coxph(Surv(t22, status) ~ X1 + X2, data = df2)
# summary(pm23)
# summary(fit2)
# 
# exp(-integrate(\(time) exp(predict(pm23,data.frame(Age = time, X1 = 0, X2 = 0))),0,3)$value)
# basehaz(fit2, centered = FALSE) %>% {stepfun(.$time, c(0,.$hazard))} %>% do.call(list(3)) %>% {exp(-.)}

