rm(list=ls())
library(ggplot2)



## Import the S&P 500 data
library(quantmod) 
SPX <- getSymbols("^GSPC", auto.assign = FALSE, from = "2012-01-03", to = "2023-12-30")
data <- as.data.frame(SPX)
price <- as.numeric(SPX$GSPC.Adjusted)
date <- as.Date(index(SPX))
log_returns <- rep(0, length(price) - 1)

for (i in 1:length(log_returns)) {
  log_returns[i] <- log(price[i+1] / price[i])
}
plot(date[-1], log_returns, type="l", xlab="Date", ylab="Daily Log Returns", main="Daily Log Returns for S&P 500 from 03/01/2012 - 29/12/2023" )



## Question a
n <- length(log_returns)
m <- 500
llh_f <- rep(0, n-m)           
llh_g <- rep(0, n-m)

# Log-likelihood of GARCH(1,1) model under Normal innovations:
negloglik_N <- function(params, y){
  alpha0 <- exp(params[1])
  alpha1 <- exp(params[2])
  alpha2 <- exp(params[3])
  n <- length(y)
  sigma2 <- c(var(y), rep(0, n))
  loglik <- 0
  for (i in 1:n) {
    sigma2[i+1] <- alpha0 + alpha1*y[i]^2 + alpha2*sigma2[i]
    loglik <- loglik + (-1/2*log(2*pi) - 1/2*log(sigma2[i]) - 1/2*y[i]^2/sigma2[i])
  }
  return(-loglik)
}

# Log-likelihood of GARCH(1,1) model under Laplace innovations:
negloglik_L <- function(params, y){
  alpha0 <- exp(params[1])
  alpha1 <- exp(params[2])
  alpha2 <- exp(params[3])
  n <- length(y)
  sigma2 <- c(var(y), rep(0, n))
  loglik <- 0
  for (i in 1:n) {
    sigma2[i+1] <- alpha0 + alpha1*y[i]^2 + alpha2*sigma2[i]
    loglik <- loglik + log(exp(-abs(y[i]) / sqrt(sigma2[i]))/(2 * sqrt(sigma2[i]))) 
  }
  return(-loglik)
}

# Rolling window process at a fixed window size of 500:
for (i in (m+1):n) {
  start <- i - m
  end <- i - 1
  sample <- log_returns[start:end]
  y_1 <- sample[1:(m-1)]
  y_2 <- sample[1:(m-2)]
  y_3 <- sample[1:(m-3)]
  y_4 <- sample[1:(m-4)]
  y_5 <- sample[1:(m-5)]
  y_1 <- y_1[-(1:4)]
  y_2 <- y_2[-(1:3)]
  y_3 <- y_3[-(1:2)]
  y_4 <- y_4[-1]
  y <- sample[6:m]
  X <- cbind(1, y_1, y_2, y_3, y_4, y_5)
  rho_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  mu <- X %*% rho_ols
  eps <- y - mu
  
  initial_wab <- c(0,0,0)
  mle_N <- optim(initial_wab, negloglik_N, y = eps)
  wab_mle_N <- exp(mle_N$par)
  h_N <- rep(var(eps), length(y))
  for (t in 1:length(y)) {
    h_N[t+1] <- wab_mle_N[1] + wab_mle_N[2]*eps[t]^2 + wab_mle_N[3]*h_N[t]
  } 
  
  mle_L <- optim(initial_wab, negloglik_L, y = eps)
  wab_mle_L <- exp(mle_L$par)
  h_L <- rep(var(eps), length(y))
  for (t in 1:length(y)) {
    h_L[t+1] <- wab_mle_L[1] + wab_mle_L[2]*eps[t]^2 + wab_mle_L[3]*h_L[t]
  } 
  
  mu_forecast_N <- rho_ols[1] + rho_ols[2] * y[length(y)] + rho_ols[3] * y[length(y)-1] + rho_ols[4] * y[length(y)-2] + rho_ols[5] * y[length(y)-3] + rho_ols[6] * y[length(y)-4]
  h_forecast_N <- wab_mle_N[1] + wab_mle_N[2] * eps[length(eps)]^2 + wab_mle_N[3] * h_N[length(h_N)]  
  llh_f[i - m] <- dnorm(log_returns[i], mu_forecast_N, sqrt(h_forecast_N), log = T)
  
  mu_forecast_L <- rho_ols[1] + rho_ols[2] * y[length(y)] + rho_ols[3] * y[length(y)-1] + rho_ols[4] * y[length(y)-2] + rho_ols[5] * y[length(y)-3] + rho_ols[6] * y[length(y)-4]
  h_forecast_L <- wab_mle_L[1] + wab_mle_L[2] * eps[length(eps)]^2 + wab_mle_L[3] * h_L[length(h_L)]  
  llh_g[i - m] <- log(exp(-abs(log_returns[i] - mu_forecast_L) / sqrt(h_forecast_L))/(2 * sqrt(h_forecast_L)))
}

llh_diff <- llh_f-llh_g
plot(1:length(llh_diff), llh_diff, type="l", col = 'red3', xlab = "Number of rolling the window of size 500", ylab = 'Log Score Differences', main = expression("Log Score Differences of forecast distributions " * f[t] * " and " * g[t]))



## Question b
# Diebold Mariano test for Normal and Laplace innovations:
auto_cov <- function(data, k){
  n = length(data)
  som = 0
  for (i in (k+1):n){
    som = som + (data[i] - mean(data))*(data[i-k] - mean(data))
  }
  return(som/n)
}

K <- floor(2517^0.25)
Bartlett <- rep(0,(K-1))
for (i in 1:(K-1)) {Bartlett[i] <- 1 - i/K}

gamma_hat_fg <- rep(0,(K-1))
for (i in 1:(K-1)) {
  gamma_hat_fg[i] <- auto_cov((llh_f-llh_g), i)
}

HAC_fg <- auto_cov(llh_f-llh_g,0) + 2 * sum(Bartlett * gamma_hat_fg)
t_dm_fg <- mean(llh_f - llh_g) / sqrt(HAC_fg / 2517)
print(t_dm_fg)



## Question c
qLaplace <- function(mu,sigma,p){mu- sqrt(0.5)*sigma * sign(p-0.5) * log(1-2*abs(p-0.5))}
alpha <- c(0.1, 0.05, 0.01)
sd_returns <- sd(log_returns)
mean_returns <- mean(log_returns)

VaR_N <- sd_returns*qnorm(alpha)
VaR_L <- qLaplace(mean_returns,sd_returns,alpha)

coverage <- function(VaR){mean(sort(log_returns) <= -VaR)}
coverage_N <- rep(0,3)
coverage_L <- rep(0,3)

for (i in 1:3) {
  coverage_N[i] <- coverage(VaR_N[i])
}
for (i in 1:3) {
  coverage_L[i] <- coverage(VaR_L[i])
}

coverage <- data.frame(Confidence_Level = alpha, Coverage_Rate_Normal = coverage_N, Coverage_Rate_Laplace = coverage_L)
print(coverage)



## Question d
llh_t <- rep(0, n-m) 

# Log-likelihood of GARCH(1,1) model under Student-t(v) innovations:
negloglik_t <- function(params, y){
  alpha0 <- exp(params[1])
  alpha1 <- exp(params[2])
  alpha2 <- exp(params[3])
  nu <- exp(params[4])
  n <- length(y)
  sigma2 <- c(var(y), rep(0, n))
  loglik <- 0
  for (i in 1:n) {
    sigma2[i+1] <- alpha0 + alpha1*y[i]^2 + alpha2*sigma2[i]
    loglik <- loglik + log(gamma((nu+1)/2) / gamma(nu/2)) - 0.5 * log(pi * nu * sigma2[i]) - ((nu+1)/2) * log(1 + 1/nu * (y[i]^2 / sigma2[i]) )
  }
  return(-loglik)
}

# Rolling window process at a fixed window size of 500:
for (i in (m+1):n) {
  start <- i - m
  end <- i - 1
  sample <- log_returns[start:end]
  y_1 <- sample[1:(m-1)]
  y_2 <- sample[1:(m-2)]
  y_3 <- sample[1:(m-3)]
  y_4 <- sample[1:(m-4)]
  y_5 <- sample[1:(m-5)]
  y_1 <- y_1[-(1:4)]
  y_2 <- y_2[-(1:3)]
  y_3 <- y_3[-(1:2)]
  y_4 <- y_4[-1]
  y <- sample[6:m]
  X <- cbind(1, y_1, y_2, y_3, y_4, y_5)
  rho_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  mu <- X %*% rho_ols
  eps <- y - mu
  
  initial_wab <- c(0,0,0,0)
  mle_t <- optim(initial_wab, negloglik_t, y = eps)
  wab_mle_t <- exp(mle_t$par[1:3])
  nu <- exp(mle_t$par[4])
  h_t <- rep(var(eps), length(y))
  for (j in 1:length(y)) {
    h_t[j+1] <- wab_mle_t[1] + wab_mle_t[2]*eps[j]^2 + wab_mle_t[3]*h_t[j]
  } 
  
  mu_forecast_t <- rho_ols[1] + rho_ols[2] * y[length(y)] + rho_ols[3] * y[length(y)-1] + rho_ols[4] * y[length(y)-2] + rho_ols[5] * y[length(y)-3] + rho_ols[6] * y[length(y)-4]
  h_forecast_t <- wab_mle_t[1] + wab_mle_t[2] * eps[length(eps)]^2 + wab_mle_t[3] * h_t[length(h_t)]  
  llh_t[i - m] <- log(gamma((nu+1)/2) / gamma(nu/2)) - 0.5 * log(pi * nu * h_forecast_t) - ((nu+1)/2) * log(1 + 1/nu * ((log_returns[i] - mu_forecast_t)^2 / h_forecast_t) )
}

plot(1:2517, (llh_t - llh_f), type="l", col = 'blue3', xlab = "Number of rolling the window of size 500", ylab = 'Log Score Differences', main = expression("Log Score Differences of forecast distributions " * h[t] * " and " * f[t]))
plot(1:2517, (llh_t - llh_g), type="l", col = 'green4', xlab = "Number of rolling the window of size 500", ylab = 'Log Score Differences', main = expression("Log Score Differences of forecast distributions " * h[t] * " and " * g[t]))

# Diebold Mariano tests
gamma_hat_hf <- rep(0,(K-1))
for (i in 1:(K-1)) {
  gamma_hat_hf[i] <- auto_cov((llh_t-llh_f), i)
}

HAC_hf <- auto_cov(llh_t-llh_f,0) + 2 * sum(Bartlett * gamma_hat_hf)
t_dm_hf <- mean(llh_t - llh_f) / sqrt(HAC_hf / 2517)
print(t_dm_hf)


gamma_hat_hg <- rep(0,(K-1))
for (i in 1:(K-1)) {
  gamma_hat_hg[i] <- auto_cov((llh_t-llh_g), i)
}

HAC_hg <- auto_cov(llh_t-llh_g,0) + 2 * sum(Bartlett * gamma_hat_hg)
t_dm_hg <- mean(llh_t - llh_g) / sqrt(HAC_hg / 2517)
print(t_dm_hg)


## Question e




