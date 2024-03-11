## Question a
rawdata <- read.csv("assignment 1.csv")
data = rawdata
data[data == "null"] <- NA 
data <- na.omit(data)
price <- as.numeric(data$Adj.Close)
date <- as.Date(data$Date)



## Question b
#Adj.Close
plot(date, price, type = 'l', xlab='Date', ylab = 'Adjusted Closing Price', main = 'Adjusted Closing Price from 01/02/2000 - 31/01/2024')


#Log-returns
log_returns <- rep(0, length(price) - 1)
for (i in 1:length(log_returns)) {
  log_returns[i] <- log(price[i+1] / price[i])
}
plot(date[-1], log_returns, type="l", xlab="Date", ylab="Daily Log Returns", main="Daily Log Returns from 01/02/2000 - 31/01/2024" )

#Rates 
rates <- rep(0, length(price)-1)
for (i in 1:length(rates)) {
  rates[i] <- (price[i+1] / price[i]) - 1
}
plot(date[-1], rates, type="l", xlab="Date", ylab="Daily Rate of Returns", main="Daily Rate of Returns from 01/02/2000 - 31/01/2024" )

#ACF
num_lag <- 6000 #try 100 and 6000
acf_values <- numeric(1 + num_lag)
mean_log_returns <- mean(log_returns)

for (k in 0:num_lag) {
  n <- length(log_returns)
  cov_value <- sum((log_returns[(k + 1):n] - mean_log_returns) * (log_returns[1:(n - k)] - mean_log_returns)) / (n - 1)  
  var_value <- sum((log_returns - mean_log_returns)^2) / (n - 1)
  acf_values[k+1] <- cov_value / var_value
}
lag_values <- 0:num_lag
plot(lag_values, acf_values, type = 'h', col = 'blue3', xlab = 'Lag', ylab = 'ACF', main = 'Autocorrelation Function (ACF) of Daily Log Returns', ylim = c(-0.1, 1))
abline(h = 0, col = 'red3', lty = 2)
points(lag_values, acf_values, pch = 19, col = 'blue3', cex=0.5)



## Question c
#OLS
return_est <- log_returns[1:(length(log_returns)-20)]
y_1 <- return_est[1:(length(return_est)-1)]
y_2 <- return_est[1:(length(return_est)-2)]
y_1 <- y_1[-1]
y = return_est[3:length(return_est)]
X = cbind(1, y_1, y_2)
phi_ols = solve(t(X) %*% X) %*% t(X) %*% y
print(phi_ols)

#MLE
log_likelihood <- function(params , X, y, sigma) { 
  phi <- params [1]
  phi_1 <- params [2]
  phi_2 <- params [3] 
  n <- length(y)
  loglik <- -(n/2)*log(2*pi) - (n/2)*log(sigma^2) - sum((y - X%*%c(phi, phi_1, phi_2))^2)*(1/(2*sigma^2)) 
  return(-loglik)
}
initial_params <- y_2[1:3] #first three values of log-returns
mle_result <- optim(initial_params , log_likelihood , X = X, y = y, sigma=var(y), method = "BFGS")
phi_mle <- mle_result$par [1:3] 
print(phi_mle)



## Question d
library(tseries)
#ar_model <- ar(return_est, order.max = 2)
#phi_ar <- c(ar_model$x.mean, ar_model$ar)
#print(phi_ar)

arima_model <- arima(return_est, order = c(2,0,0))
phi_arima <- arima_model$coef
print(phi_arima)

#fitted AR (unnecessary to present this model, we can skip it and only focus on ARIMA in d-f)
#y_ar <- X %*% phi_ar
#plot(date[1:length(y)], y, type="l", xlab="Date", ylab="Daily Log Returns", main="Daily Log Returns from 01/02/2000 - 31/01/2024" )
#lines(date[1:length(y)], y_ar, col='yellow')

#fitted ARIMA
y_arima <- X %*% c(phi_arima[3], phi_arima[1], phi_arima[2])
plot(date[1:length(y)], y, type="l", xlab="Date", ylab="Daily Log Returns", main="Fitted Returns (ARIMA) and Actual Returns" )
lines(date[1:length(y)], y_arima, col='deeppink3')
legend("bottomleft", 
       legend=c("Actual Returns", "Fitted Returns (ARIMA)"), 
       col=c("black", "deeppink3"), 
       lty=1, 
       cex=0.8)

#fitted OLS
y_ols <- X %*% phi_ols
plot(date[1:length(y)], y, type="l", xlab="Date", ylab="Daily Log Returns", main="Fitted Returns (OLS) and Actual Returns")
lines(date[1:length(y)], y_ols, col="deepskyblue3")
legend("bottomleft", 
       legend=c("Actual Returns", "Fitted Returns (OLS)"), 
       col=c("black", "deepskyblue3"), 
       lty=1, 
       cex=0.8)

#fitted MLE
y_mle <- X %*% phi_mle
plot(date[1:length(y)], y, type="l", xlab="Date", ylab="Daily Log Returns", main="Fitted Returns (MLE) and Actual Returns")
lines(date[1:length(y)], y_mle, col="green3")
legend("bottomleft", 
       legend=c("Actual Returns", "Fitted Returns (MLE)"), 
       col=c("black", "green3"), 
       lty=1, 
       cex=0.8)



## Question e (present the result of JB-Normality test is enough, Q-Q plot unnecessary)
#ARIMA
e_arima = arima_model$residuals
qqnorm(e_arima, main = "Normal Q-Q plot (ARIMA)")
qqline(e_arima)

jb_arima <- jarque.bera.test(e_arima)
jb_arima_pval <- jb_arima$p.value
print(jb_arima_pval)

#OLS
e_ols = y - y_ols
qqnorm(e_ols, main = "Normal Q-Q plot (OLS)")
qqline(e_ols)

jb_ols<- jarque.bera.test(e_ols)
jb_ols_pval <- jb_ols$p.value
print(jb_ols_pval)

#MLE
e_mle = y - y_mle
qqnorm(e_mle, main = "Normal Q-Q plot (MLE)")
qqline(e_mle)

jb_mle <- jarque.bera.test(e_mle)
jb_mle_pval <- jb_mle$p.value
print(jb_mle_pval)
   


## Question f
#ARIMA for the forecast
y_new <- c(return_est, rep(0,20))
for (i in 1:20) {
  y_new[length(return_est)+i] <- phi_arima[3] + y_new[length(return_est) - 1 + i] * phi_arima[1] + y_new[length(return_est) - 2 + i] * phi_arima[2]
}
print(y_new[(length(return_est)+1):(length(return_est)+20)])

plot(tail(date, 20), tail(log_returns,20), type = "l", main = "Forecasted Returns (ARIMA) and Actual Log Returns over Jan 2024", xlab="Date", ylab="Daily Log Returns")
points(tail(date, 20), tail(log_returns, 20), pch = 16, col = "red3")
lines(tail(date, 20), tail(y_new, 20))
points(tail(date, 20), tail(y_new, 20), pch=16, col='blue3')
legend("bottomleft", 
       legend=c("Actual Returns", "Forecasted Returns (ARIMA)"), 
       col=c("red3", "blue3"), 
       lty=1, 
       cex=0.8)

