### This the submission for Group 6: Bowen Ma (12960780), Mianyun He (13605275), Shiyi Yang (13627295)
rm(list=ls())

## Question a
rawdata <- read.csv("assignment 2.csv")
data = rawdata
data[data == "null"] <- NA 
data <- na.omit(data)
price <- as.numeric(data$Adj.Close)
date <- as.Date(data$Date)

log_returns <- rep(0, length(price) - 1)
for (i in 1:length(log_returns)) {
  log_returns[i] <- log(price[i+1] / price[i])
}
plot(date[-1], log_returns, type="l", xlab="Date", ylab="Daily Log Returns", main="Daily Log Returns from 02/01/2002 - 29/12/2023" )



## Question b
y <- log_returns
min_aic <- Inf
best_p <- NULL
for (p in c(1, 2, 3, 4, 5)) {
  result <- arma(y, order = c(p,0))
  b_arma <- result$coef
  resid_var <- var(result$resid, na.rm=TRUE)
  n <- length(y)
  l <- -(n/2) * log(2 * pi * resid_var) - (1/(2*resid_var)) * sum(result$resid^2, na.rm=TRUE) 
  n <- length(y)
  k <- length(b_arma)
  AIC <- -2 * l + 2 * k
  BIC <- -2 * l + k * log(n)
  if (AIC < min_aic) {
    min_aic <- AIC
    best_p <- p
  }
  cat(sprintf('AR(%d) model: log-likelihood =%5.3f , AIC =%7.3f, BIC =%7.3f\n', p, l, AIC, BIC))
}
cat("Best model: AR(", best_p, ") with AIC =", min_aic, "\n")



## Question c
RMSE_rolling <- rep(0,5)
N = length(log_returns)
m = 750
in_sample <- log_returns[1:m]
out_sample <- log_returns[(m+1) : N]    

# AR(1)
y_1 <- in_sample[1:(m-1)]
y <- in_sample[2:m] 
X <- cbind(1, y_1)
phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
y_new <- c(in_sample, rep(0,(N-m)))
for (i in 1:(N-m)) {
  y_new[m+i] <- phi_ols[1] * 1 + y_new[m + i - 1] * phi_ols[2] 
}
RMSE_rolling[1] <- sqrt((1/(N-m))*sum((out_sample - y_new[(m+1):N])^2))


# AR(2)
y_1 <- in_sample[1:(m-1)]
y_2 <- in_sample[1:(m-2)]
y_1 <- y_1[-1]
y <- in_sample[3:m]
X <- cbind(1, y_1, y_2)
phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
y_new <- c(in_sample, rep(0,(N-m)))
for (i in 1:(N-m)) {
  y_new[m+i] <- phi_ols[1] * 1 + y_new[m + i - 1] * phi_ols[2] + y_new[m + i - 2] * phi_ols[3]
}
RMSE_rolling[2] <- sqrt((1/(N-m))*sum((out_sample - y_new[(m+1):N])^2))


# AR(3)
y_1 <- in_sample[1:(m-1)]
y_2 <- in_sample[1:(m-2)]
y_3 <- in_sample[1:(m-3)]
y_1 <- y_1[-(1:2)]
y_2 <- y_2[-1]
y <- in_sample[4:m]
X <- cbind(1, y_1, y_2, y_3)
phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
y_new <- c(in_sample, rep(0,(N-m)))
for (i in 1:(N-m)) {
  y_new[m+i] <- phi_ols[1] * 1 + y_new[m + i - 1] * phi_ols[2] + y_new[m + i - 2] * phi_ols[3] + y_new[m + i - 3] * phi_ols[4]
}
RMSE_rolling[3] <- sqrt((1/(N-m))*sum((out_sample - y_new[(m+1):N])^2))


# AR(4)
y_1 <- in_sample[1:(m-1)]
y_2 <- in_sample[1:(m-2)]
y_3 <- in_sample[1:(m-3)]
y_4 <- in_sample[1:(m-4)]
y_1 <- y_1[-(1:3)]
y_2 <- y_2[-(1:2)]
y_3 <- y_3[-1]
y <- in_sample[5:m]
X <- cbind(1, y_1, y_2, y_3, y_4)
phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
y_new <- c(in_sample, rep(0,(N-m)))
for (i in 1:(N-m)) {
  y_new[m+i] <- phi_ols[1] * 1 + y_new[m + i - 1] * phi_ols[2] + y_new[m + i - 2] * phi_ols[3] + y_new[m + i - 3] * phi_ols[4] + y_new[m + i - 4] * phi_ols[5]
}
RMSE_rolling[4] <- sqrt((1/(N-m))*sum((out_sample - y_new[(m+1):N])^2))


# AR(5)
y_1 <- in_sample[1:(m-1)]
y_2 <- in_sample[1:(m-2)]
y_3 <- in_sample[1:(m-3)]
y_4 <- in_sample[1:(m-4)]
y_5 <- in_sample[1:(m-5)]
y_1 <- y_1[-(1:4)]
y_2 <- y_2[-(1:3)]
y_3 <- y_3[-(1:2)]
y_4 <- y_4[-1]
y <- in_sample[6:m]
X <- cbind(1, y_1, y_2, y_3, y_4, y_5)
phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
y_new <- c(in_sample, rep(0,(N-m)))
for (i in 1:(N-m)) {
  y_new[m+i] <- phi_ols[1] * 1 + y_new[m + i - 1] * phi_ols[2] + y_new[m + i - 2] * phi_ols[3] + y_new[m + i - 3] * phi_ols[4] + y_new[m + i - 4] * phi_ols[5] + y_new[m + i - 5] * phi_ols[6]
}
RMSE_rolling[5] <- sqrt((1/(N-m))*sum((out_sample - y_new[(m+1):N])^2))

print(RMSE_rolling)



## Question d
RMSE_rolling_est <- rep(0, 5)
N <- length(log_returns)
m <- 750
out_sample <- log_returns[(m + 1):N]    


# AR(1)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[j:(m + j - 1)]
  y <- in_sample[2:m]
  y_1 <- in_sample[1:(m - 1)]
  X <- cbind(1, y_1)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  if (j == 1) {
    y_forecast[j] <- phi_ols[1] * 1 + in_sample[m] * phi_ols[2]
  } else {
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2]
  }
}
RMSE_rolling_est[1] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))


# AR(2)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[j:(m + j - 1)]
  y_1 <- in_sample[1:(m - 1)]
  y_2 <- in_sample[1:(m - 2)]
  y_1 <- y_1[-1]
  y <- in_sample[3:m]
  X <- cbind(1, y_1, y_2)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  if (j == 1) {
    y_forecast[j] <- phi_ols[1] * 1 + in_sample[m] * phi_ols[2] + in_sample[m - 1] * phi_ols[3]
  } else if (j == 2) {
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + in_sample[m] * phi_ols[3]
  } else {
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3]
  }
}
RMSE_rolling_est[2] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))


# AR(3)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[j:(m+j-1)]
  y_1 <- in_sample[1:(m-1)]
  y_2 <- in_sample[1:(m-2)]
  y_3 <- in_sample[1:(m-3)]
  y_1 <- y_1[-(1:2)]
  y_2 <- y_2[-1]
  y <- in_sample[4:m]
  X <- cbind(1, y_1, y_2, y_3)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  if (j == 1) {
    y_forecast[j] <- phi_ols[1] * 1 + in_sample[m] * phi_ols[2] + in_sample[m - 1] * phi_ols[3] + in_sample[m-2] * phi_ols[4]
  } else if (j == 2) {
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + in_sample[m] * phi_ols[3] + in_sample[m-1] * phi_ols[4]
  } else if (j == 3){
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + in_sample[m] * phi_ols[4]
  } else{
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + y_forecast[j - 3] * phi_ols[4]
  }
}
RMSE_rolling_est[3] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))


# AR(4)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[j:(m+j-1)]
  y_1 <- in_sample[1:(m-1)]
  y_2 <- in_sample[1:(m-2)]
  y_3 <- in_sample[1:(m-3)]
  y_4 <- in_sample[1:(m-4)]
  y_1 <- y_1[-(1:3)]
  y_2 <- y_2[-(1:2)]
  y_3 <- y_3[-1]
  y <- in_sample[5:m]
  X <- cbind(1, y_1, y_2, y_3, y_4)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  if (j == 1) {
    y_forecast[j] <- phi_ols[1] * 1 + in_sample[m] * phi_ols[2] + in_sample[m - 1] * phi_ols[3] + in_sample[m-2] * phi_ols[4] + in_sample[m-3] * phi_ols[5]
  } else if (j == 2) {
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + in_sample[m] * phi_ols[3] + in_sample[m-1] * phi_ols[4] + in_sample[m-2] * phi_ols[5]
  } else if (j == 3){
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + in_sample[m] * phi_ols[4] + in_sample[m-1] * phi_ols[5]
  } else if (j == 4){
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + y_forecast[j - 3] * phi_ols[4] + in_sample[m] * phi_ols[5]
  } else{
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + y_forecast[j - 3] * phi_ols[4] + y_forecast[j - 4] * phi_ols[5]
  }
}
RMSE_rolling_est[4] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))


# AR(5)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[j:(m+j-1)]
  y_1 <- in_sample[1:(m-1)]
  y_2 <- in_sample[1:(m-2)]
  y_3 <- in_sample[1:(m-3)]
  y_4 <- in_sample[1:(m-4)]
  y_5 <- in_sample[1:(m-5)]
  y_1 <- y_1[-(1:4)]
  y_2 <- y_2[-(1:3)]
  y_3 <- y_3[-(1:2)]
  y_4 <- y_4[-1]
  y <- in_sample[6:m]
  X <- cbind(1, y_1, y_2, y_3, y_4, y_5)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  if (j == 1) {
    y_forecast[j] <- phi_ols[1] * 1 + in_sample[m] * phi_ols[2] + in_sample[m - 1] * phi_ols[3] + in_sample[m-2] * phi_ols[4] + in_sample[m-3] * phi_ols[5] + in_sample[m-4] * phi_ols[6]
  } else if (j == 2) {
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + in_sample[m] * phi_ols[3] + in_sample[m-1] * phi_ols[4] + in_sample[m-2] * phi_ols[5] + in_sample[m-3] * phi_ols[6]
  } else if (j == 3){
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + in_sample[m] * phi_ols[4] + in_sample[m-1] * phi_ols[5] + in_sample[m-2] * phi_ols[6]
  } else if (j == 4){
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + y_forecast[j - 3] * phi_ols[4] + in_sample[m] * phi_ols[5] + in_sample[m-1] * phi_ols[6]
  } else if (j == 5){
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + y_forecast[j - 3] * phi_ols[4] + y_forecast[j - 4] * phi_ols[5] + in_sample[m] * phi_ols[6]
  } else{
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + y_forecast[j - 3] * phi_ols[4] + y_forecast[j - 4] * phi_ols[5] + y_forecast[j - 5] * phi_ols[6]
  }
}
RMSE_rolling_est[5] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))
print(RMSE_rolling_est)



## Question e
RMSE_rolling_exp <- rep(0, 5)
N <- length(log_returns)
m <- 750
out_sample <- log_returns[(m + 1):N]    


# AR(1)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[1:(m + j - 1)]
  y <- in_sample[2:(m + j -1)]
  y_1 <- in_sample[1:(m + j - 2)]
  X <- cbind(1, y_1)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  if (j == 1) {
    y_forecast[j] <- phi_ols[1] * 1 + in_sample[m + j - 1] * phi_ols[2]
  } else {
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2]
  }
}
RMSE_rolling_exp[1] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))


# AR(2)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[1:(m + j - 1)]
  y_1 <- in_sample[1:(m + j - 2)]
  y_2 <- in_sample[1:(m + j - 3)]
  y_1 <- y_1[-1]
  y <- in_sample[3:(m + j - 1)]
  X <- cbind(1, y_1, y_2)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  if (j == 1) {
    y_forecast[j] <- phi_ols[1] * 1 + in_sample[(m+j-1)] * phi_ols[2] + in_sample[(m+j-2)] * phi_ols[3]
  } else if (j == 2) {
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + in_sample[(m+j-2)] * phi_ols[3]
  } else {
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3]
  }
}
RMSE_rolling_exp[2] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))


# AR(3)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[1:(m + j - 1)]
  y_1 <- in_sample[1:(m + j - 2)]
  y_2 <- in_sample[1:(m + j - 3)]
  y_3 <- in_sample[1:(m + j - 4)]
  y_1 <- y_1[-(1:2)]
  y_2 <- y_2[-1]
  y <- in_sample[4:(m + j - 1)]
  X <- cbind(1, y_1, y_2, y_3)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  if (j == 1) {
    y_forecast[j] <- phi_ols[1] * 1 + in_sample[(m + j - 1)] * phi_ols[2] + in_sample[(m + j - 2)] * phi_ols[3] + in_sample[(m + j - 3)] * phi_ols[4]
  } else if (j == 2) {
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + in_sample[(m + j - 2)] * phi_ols[3] + in_sample[(m + j - 3)] * phi_ols[4]
  } else if (j == 3){
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + in_sample[(m + j - 3)] * phi_ols[4]
  } else{
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + y_forecast[j - 3] * phi_ols[4]
  }
}
RMSE_rolling_exp[3] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))


# AR(4)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[1:(m+j-1)]
  y_1 <- in_sample[1:(m + j - 2)]
  y_2 <- in_sample[1:(m + j - 3)]
  y_3 <- in_sample[1:(m + j - 4)]
  y_4 <- in_sample[1:(m + j - 5)]
  y_1 <- y_1[-(1:3)]
  y_2 <- y_2[-(1:2)]
  y_3 <- y_3[-1]
  y <- in_sample[5:(m+j-1)]
  X <- cbind(1, y_1, y_2, y_3, y_4)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  if (j == 1) {
    y_forecast[j] <- phi_ols[1] * 1 + in_sample[(m+j-1)] * phi_ols[2] + in_sample[(m+j-2)] * phi_ols[3] + in_sample[(m+j-3)] * phi_ols[4] + in_sample[(m+j-4)] * phi_ols[5]
  } else if (j == 2) {
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + in_sample[(m+j-2)] * phi_ols[3] + in_sample[(m+j-3)] * phi_ols[4] + in_sample[(m+j-4)] * phi_ols[5]
  } else if (j == 3){
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + in_sample[(m+j-3)] * phi_ols[4] + in_sample[(m+j-4)] * phi_ols[5]
  } else if (j == 4){
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + y_forecast[j - 3] * phi_ols[4] + in_sample[(m+j-4)] * phi_ols[5]
  } else{
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + y_forecast[j - 3] * phi_ols[4] + y_forecast[j - 4] * phi_ols[5]
  }
}
RMSE_rolling_exp[4] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))


# AR(5)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[1:(m+j-1)]
  y_1 <- in_sample[1:(m+j-2)]
  y_2 <- in_sample[1:(m+j-3)]
  y_3 <- in_sample[1:(m+j-4)]
  y_4 <- in_sample[1:(m+j-5)]
  y_5 <- in_sample[1:(m+j-6)]
  y_1 <- y_1[-(1:4)]
  y_2 <- y_2[-(1:3)]
  y_3 <- y_3[-(1:2)]
  y_4 <- y_4[-1]
  y <- in_sample[6:(m+j-1)]
  X <- cbind(1, y_1, y_2, y_3, y_4, y_5)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  if (j == 1) {
    y_forecast[j] <- phi_ols[1] * 1 + in_sample[(m+j-1)] * phi_ols[2] + in_sample[(m+j-2)] * phi_ols[3] + in_sample[(m+j-3)] * phi_ols[4] + in_sample[(m+j-4)] * phi_ols[5] + in_sample[(m+j-5)] * phi_ols[6]
  } else if (j == 2) {
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + in_sample[(m+j-2)] * phi_ols[3] + in_sample[(m+j-3)] * phi_ols[4] + in_sample[(m+j-4)] * phi_ols[5] + in_sample[(m+j-5)] * phi_ols[6]
  } else if (j == 3){
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + in_sample[(m+j-3)] * phi_ols[4] + in_sample[(m+j-4)] * phi_ols[5] + in_sample[(m+j-5)] * phi_ols[6]
  } else if (j == 4){
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + y_forecast[j - 3] * phi_ols[4] + in_sample[(m+j-4)] * phi_ols[5] + in_sample[(m+j-5)] * phi_ols[6]
  } else if (j == 5){
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + y_forecast[j - 3] * phi_ols[4] + y_forecast[j - 4] * phi_ols[5] + in_sample[(m+j-5)] * phi_ols[6]
  } else{
    y_forecast[j] <- phi_ols[1] * 1 + y_forecast[j - 1] * phi_ols[2] + y_forecast[j - 2] * phi_ols[3] + y_forecast[j - 3] * phi_ols[4] + y_forecast[j - 4] * phi_ols[5] + y_forecast[j - 5] * phi_ols[6]
  }
}
RMSE_rolling_exp[5] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))
print(RMSE_rolling_exp)


## Combine the results from c-e in one table:
RMSE_matrix <- matrix(NA, nrow = 5, ncol = 3)

RMSE_matrix[, 1] <- c(RMSE_rolling[1], RMSE_rolling[2], RMSE_rolling[3], RMSE_rolling[4], RMSE_rolling[5])
RMSE_matrix[, 2] <- c(RMSE_rolling_est[1], RMSE_rolling_est[2], RMSE_rolling_est[3], RMSE_rolling_est[4], RMSE_rolling_est[5])
RMSE_matrix[, 3] <- c(RMSE_rolling_exp[1], RMSE_rolling_exp[2], RMSE_rolling_exp[3], RMSE_rolling_exp[4], RMSE_rolling_exp[5])

colnames(RMSE_matrix) <- c("RMSE_rolling_fix", "RMSE_rolling_re-estimate", "RMSE_rolling_expand")
rownames(RMSE_matrix) <- c("1 lag", "2 lags", "3 lags", "4 lags", "5 lags")

print(RMSE_matrix)



#============================================Reflection=========================================================

## Question d
RMSE_rolling_est <- rep(0, 5)
N <- length(log_returns)
m <- 750
out_sample <- log_returns[(m + 1):N]    


# AR(1)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[j:(m + j - 1)]
  y <- in_sample[2:m]
  y_1 <- in_sample[1:(m - 1)]
  X <- cbind(1, y_1)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  y_forecast[j] <- phi_ols[1] * 1 + in_sample[m] * phi_ols[2]
  
}
RMSE_rolling_est[1] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))
print(RMSE_rolling_est)


# AR(2)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[j:(m + j - 1)]
  y_1 <- in_sample[1:(m - 1)]
  y_2 <- in_sample[1:(m - 2)]
  y_1 <- y_1[-1]
  y <- in_sample[3:m]
  X <- cbind(1, y_1, y_2)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  y_forecast[j] <- phi_ols[1] * 1 + in_sample[m] * phi_ols[2] + in_sample[m - 1] * phi_ols[3]
}
RMSE_rolling_est[2] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))
print(RMSE_rolling_est)


# AR(3)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[j:(m+j-1)]
  y_1 <- in_sample[1:(m-1)]
  y_2 <- in_sample[1:(m-2)]
  y_3 <- in_sample[1:(m-3)]
  y_1 <- y_1[-(1:2)]
  y_2 <- y_2[-1]
  y <- in_sample[4:m]
  X <- cbind(1, y_1, y_2, y_3)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  y_forecast[j] <- phi_ols[1] * 1 + in_sample[m] * phi_ols[2] + in_sample[m - 1] * phi_ols[3] + in_sample[m-2] * phi_ols[4]
}
RMSE_rolling_est[3] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))
print(RMSE_rolling_est)


# AR(4)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[j:(m+j-1)]
  y_1 <- in_sample[1:(m-1)]
  y_2 <- in_sample[1:(m-2)]
  y_3 <- in_sample[1:(m-3)]
  y_4 <- in_sample[1:(m-4)]
  y_1 <- y_1[-(1:3)]
  y_2 <- y_2[-(1:2)]
  y_3 <- y_3[-1]
  y <- in_sample[5:m]
  X <- cbind(1, y_1, y_2, y_3, y_4)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  y_forecast[j] <- phi_ols[1] * 1 + in_sample[m] * phi_ols[2] + in_sample[m - 1] * phi_ols[3] + in_sample[m-2] * phi_ols[4] + in_sample[m-3] * phi_ols[5]
}
RMSE_rolling_est[4] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))
print(RMSE_rolling_est)


# AR(5)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[j:(m+j-1)]
  y_1 <- in_sample[1:(m-1)]
  y_2 <- in_sample[1:(m-2)]
  y_3 <- in_sample[1:(m-3)]
  y_4 <- in_sample[1:(m-4)]
  y_5 <- in_sample[1:(m-5)]
  y_1 <- y_1[-(1:4)]
  y_2 <- y_2[-(1:3)]
  y_3 <- y_3[-(1:2)]
  y_4 <- y_4[-1]
  y <- in_sample[6:m]
  X <- cbind(1, y_1, y_2, y_3, y_4, y_5)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  y_forecast[j] <- phi_ols[1] * 1 + in_sample[m] * phi_ols[2] + in_sample[m - 1] * phi_ols[3] + in_sample[m-2] * phi_ols[4] + in_sample[m-3] * phi_ols[5] + in_sample[m-4] * phi_ols[6]
}
RMSE_rolling_est[5] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))
print(RMSE_rolling_est)



## Question e
RMSE_rolling_exp <- rep(0, 5)
N <- length(log_returns)
m <- 750
out_sample <- log_returns[(m + 1):N]    


# AR(1)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[1:(m + j - 1)]
  y <- in_sample[2:(m + j -1)]
  y_1 <- in_sample[1:(m + j - 2)]
  X <- cbind(1, y_1)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  y_forecast[j] <- phi_ols[1] * 1 + in_sample[m + j - 1] * phi_ols[2]
  
}
RMSE_rolling_exp[1] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))
print(RMSE_rolling_exp)


# AR(2)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[1:(m + j - 1)]
  y_1 <- in_sample[1:(m + j - 2)]
  y_2 <- in_sample[1:(m + j - 3)]
  y_1 <- y_1[-1]
  y <- in_sample[3:(m + j - 1)]
  X <- cbind(1, y_1, y_2)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  y_forecast[j] <- phi_ols[1] * 1 + in_sample[(m+j-1)] * phi_ols[2] + in_sample[(m+j-2)] * phi_ols[3]
  
}
RMSE_rolling_exp[2] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))
print(RMSE_rolling_exp)


# AR(3)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[1:(m + j - 1)]
  y_1 <- in_sample[1:(m + j - 2)]
  y_2 <- in_sample[1:(m + j - 3)]
  y_3 <- in_sample[1:(m + j - 4)]
  y_1 <- y_1[-(1:2)]
  y_2 <- y_2[-1]
  y <- in_sample[4:(m + j - 1)]
  X <- cbind(1, y_1, y_2, y_3)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  y_forecast[j] <- phi_ols[1] * 1 + in_sample[(m + j - 1)] * phi_ols[2] + in_sample[(m + j - 2)] * phi_ols[3] + in_sample[(m + j - 3)] * phi_ols[4]
}
RMSE_rolling_exp[3] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))
print(RMSE_rolling_exp)


# AR(4)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[1:(m+j-1)]
  y_1 <- in_sample[1:(m + j - 2)]
  y_2 <- in_sample[1:(m + j - 3)]
  y_3 <- in_sample[1:(m + j - 4)]
  y_4 <- in_sample[1:(m + j - 5)]
  y_1 <- y_1[-(1:3)]
  y_2 <- y_2[-(1:2)]
  y_3 <- y_3[-1]
  y <- in_sample[5:(m+j-1)]
  X <- cbind(1, y_1, y_2, y_3, y_4)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  y_forecast[j] <- phi_ols[1] * 1 + in_sample[(m+j-1)] * phi_ols[2] + in_sample[(m+j-2)] * phi_ols[3] + in_sample[(m+j-3)] * phi_ols[4] + in_sample[(m+j-4)] * phi_ols[5]
  
}
RMSE_rolling_exp[4] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))
print(RMSE_rolling_exp)


# AR(5)
y_forecast <- rep(0, (N - m))
for (j in 1:(N - m)) {
  in_sample <- log_returns[1:(m+j-1)]
  y_1 <- in_sample[1:(m+j-2)]
  y_2 <- in_sample[1:(m+j-3)]
  y_3 <- in_sample[1:(m+j-4)]
  y_4 <- in_sample[1:(m+j-5)]
  y_5 <- in_sample[1:(m+j-6)]
  y_1 <- y_1[-(1:4)]
  y_2 <- y_2[-(1:3)]
  y_3 <- y_3[-(1:2)]
  y_4 <- y_4[-1]
  y <- in_sample[6:(m+j-1)]
  X <- cbind(1, y_1, y_2, y_3, y_4, y_5)
  phi_ols <- solve(t(X) %*% X) %*% t(X) %*% y
  y_forecast[j] <- phi_ols[1] * 1 + in_sample[(m+j-1)] * phi_ols[2] + in_sample[(m+j-2)] * phi_ols[3] + in_sample[(m+j-3)] * phi_ols[4] + in_sample[(m+j-4)] * phi_ols[5] + in_sample[(m+j-5)] * phi_ols[6]
}
RMSE_rolling_exp[5] <- sqrt((1 / (N - m)) * sum((out_sample - y_forecast)^2))
print(RMSE_rolling_exp)


## Combine the results from c-e in one table:
RMSE_matrix <- matrix(NA, nrow = 5, ncol = 3)

RMSE_matrix[, 1] <- c(RMSE_rolling[1], RMSE_rolling[2], RMSE_rolling[3], RMSE_rolling[4], RMSE_rolling[5])
RMSE_matrix[, 2] <- c(RMSE_rolling_est[1], RMSE_rolling_est[2], RMSE_rolling_est[3], RMSE_rolling_est[4], RMSE_rolling_est[5])
RMSE_matrix[, 3] <- c(RMSE_rolling_exp[1], RMSE_rolling_exp[2], RMSE_rolling_exp[3], RMSE_rolling_exp[4], RMSE_rolling_exp[5])

colnames(RMSE_matrix) <- c("RMSE_rolling_fix", "RMSE_rolling_re-estimate", "RMSE_rolling_expand")
rownames(RMSE_matrix) <- c("1 lag", "2 lags", "3 lags", "4 lags", "5 lags")

print(RMSE_matrix)

