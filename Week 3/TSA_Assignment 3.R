### This the submission for Group 6: Bowen Ma (12960780), Mianyun He (13605275), Shiyi Yang (13627295)
rm(list=ls())
library(ggplot2)
## Question a
rawdata <- read.csv("assignment 3.csv")
data = rawdata
data[data == "null"] <- NA 
data <- na.omit(data)
price <- as.numeric(data$Adj.Close)
date <- as.Date(data$Date, format = "%Y-%m-%d")



## Question b
log_returns <- rep(0, length(price) - 1)
for (i in 1:length(log_returns)) {
  log_returns[i] <- log(price[i+1] / price[i])
}
ggplot(data.frame(Date = date[-1], Log_Return = log_returns), aes(x = Date, y = Log_Return)) +
  geom_line() +scale_x_date(
    breaks = seq(as.Date("2000-01-01"), as.Date("2020-01-01"), by = "5 years"),
    date_labels = "%Y"
  ) +
  labs(title = "Daily Log Returns from 01/02/2000 - 31/01/2024",
       x = "Date",
       y = "Daily Log Returns")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                      panel.background = element_blank(), axis.line = element_line(colour = "black"))


log_returns_sqr <- log_returns^2
ggplot(data.frame(Date=date[-1],S_log_return=log_returns_sqr),aes(x=Date,y=S_log_return))+geom_line()+
  scale_x_date(
    breaks = seq(as.Date("2000-01-01"), as.Date("2020-01-01"), by = "5 years"),
    date_labels = "%Y"
  ) +
  labs(title = "Daily Squared Log Returns from 01/02/2000 - 31/01/2024",
       x = "Date",
       y = "Daily Squared Log Returns")+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                              panel.background = element_blank(), axis.line = element_line(colour = "black"))



## Question c
return_est <- log_returns[1:(length(log_returns)-20)]
y <- c(0,return_est)
n <- length(y)

negloglik <- function(params, y){
  alpha0 <- exp(params[1])
  alpha1 <- exp(params[2])
  alpha2 <- exp(params[3])
  sigma2 <- c(var(y), rep(0, n))
  loglik <- 0
  for (i in 1:n) {
    sigma2[i+1] <- alpha0 + alpha1*y[i]^2 + alpha2*sigma2[i]
    loglik <- loglik + (-1/2*log(2*pi) - 1/2*log(sigma2[i]) - 1/2*y[i]^2/sigma2[i])
  }
  return(-loglik)
}

initial_params <- c(0,0,0)
mle_c <- optim(initial_params, negloglik, y = y)
mle_c_results <- exp(mle_c$par)
mle_c_results
mle_c_llh <- -negloglik(mle_c$par, y)      # check the value of the log-likelihood in c)
mle_c_llh




## Question d
# qqnorm(return_est, main = "Normal Q-Q plot (Daily Log Returns)")
# qqline(return_est)

skewness<-function(y){
  n<-length(y)
  result<-(sum((y-mean(y))^3/n))/(sum((y-mean(y))^2/n)^(3/2))
  print(result)
}

kurtosis<-function(y){
  n<-length(y)
  result<-(sum((y-mean(y))^4/n))/(sum((y-mean(y))^2/n)^(2))
  print(result)
}

jb_stat <- (length(return_est)/6) * (skewness(return_est)^2 + (1/4) * (kurtosis(return_est)-3)^2)
jb_stat

mu <- mean(return_est)
sigma <- sd(return_est)
x <- seq(-0.2, 0.2, length=100)
y <- dnorm(x, mean=mu, sd=sigma)
hist(return_est, main="Histogram with Density Plot", xlab="Value", freq=FALSE, xlim=c(-0.2, 0.2),ylim=c(0,50))
lines(density(return_est), col="chocolate", lwd=2)
lines(x, y, col="deepskyblue", lwd=2)
legend("topleft", 
       legend=c("Data probability density function", "Normal distribution"), 
       col=c("chocolate", "deepskyblue3"), 
       lty=1,lwd=1.5, 
       cex=0.6)

## Question e
library("rugarch")

garch_e <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                      mean.model = list(armaOrder = c(0, 0), include.mean = F), 
                      distribution.model = "norm")

#garch_e <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
#mean.model = list(armaOrder = c(0, 0), include.mean = F), 
#distribution.model = "std")

garch_e_result <- ugarchfit(garch_e, data = return_est)
garch_e_result
garch_e_coef <- garch_e_result@fit$coef
garch_e_coef
garch_e_llh <- garch_e_result@fit$LLH      # check the value of the log-likelihood in e)
garch_e_llh                           

res_sqr <- residuals(garch_e_result)^2
convar_fit <- sigma(garch_e_result)^2
ggplot(data.frame(date = date[1:length(res_sqr)], res_sqr = res_sqr, convar_fit = convar_fit), aes(x = date)) +
  geom_line(aes(y = res_sqr), color = "black", linetype = "solid") +
  geom_line(aes(y = convar_fit), color = "red", linetype = "solid") +
  labs(x = "Date", y = "Squared Residuals", title = "Fitted Conditional Variances and The Squared Residuals") +
  theme_minimal() +scale_x_date(
    breaks = seq(as.Date("2000-01-01"), as.Date("2020-01-01"), by = "5 years"),
    date_labels = "%Y"
  )

ggplot(data.frame(date = date[1:length(res_sqr)], res_sqr = res_sqr, convar_fit = convar_fit), aes(x = date)) +
  geom_line(aes(y = res_sqr, color = "Squared Residuals"), linetype = "solid") +
  geom_line(aes(y = convar_fit, color = "Fitted Conditional Variances"), linetype = "solid") +
  labs(x = "Date", y = "Squared Residuals", title = "Fitted Conditional Variances and The Squared Residuals") +
  theme_minimal() + 
  scale_x_date(
    breaks = seq(as.Date("2000-01-01"), as.Date("2020-01-01"), by = "5 years"),
    date_labels = "%Y"
  ) +
  scale_color_manual(values = c("Squared Residuals" = "black", "Fitted Conditional Variances" = "red")) +
  theme(legend.position = "top-left")


## Question f
y <- c(0,return_est)
sigma2 <- c(var(y), rep(0, n))
for (i in 1:n) {
  sigma2[i+1] <- mle_c_results[1] + mle_c_results[2]*y[i]^2 + mle_c_results[3]*sigma2[i]
}

sigma2_forecast <- c(sigma2[length(sigma2)], rep(0,20))
for (i in 1:20) {
  sigma2_forecast[i+1] <- mle_c_results[1] + (mle_c_results[2] + mle_c_results[3]) * sigma2_forecast[i]
}

plot(tail(date, 20), tail(log_returns_sqr, 20), type = "l", 
     xlab = "Date", ylab = "Squared Returns", 
     main = "Forecasted Conditional Variances and The Squared Returns")
lines(tail(date,20), tail(sigma2_forecast,20), col="deeppink3")
legend("topleft", 
       legend=c("The Actual Squared Returns", "Forecasted Conditional Variances"), 
       col=c("black", "deeppink3"), 
       lty=1, 
       cex=0.8)
  
## Question g
rm(list=ls())

# 1)
rawdata <- read.csv("assignment 3_BTC.csv")
data = rawdata
data[data == "null"] <- NA 
data <- na.omit(data)
price <- as.numeric(data$Adj.Close)
date <- as.Date(data$Date, format = "%Y-%m-%d")


# 2)
log_returns <- rep(0, length(price) - 1)
for (i in 1:length(log_returns)) {
  log_returns[i] <- log(price[i+1] / price[i])
}
ggplot(data.frame(Date = date[-1], Log_Return = log_returns), aes(x = Date, y = Log_Return)) +
  geom_line() +
  scale_x_date(
    breaks = seq(as.Date("2014-10-01"), as.Date("2022-01-01"), by = "2 years"),
    date_labels = "%Y"
  ) +
  labs(title = "Daily Log Returns from 17/09/2014 - 20/01/2023",
       x = "Date",
       y = "Daily Log Returns") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

log_returns_sqr <- log_returns^2
ggplot(data.frame(Date = date[-1], Log_Return = log_returns_sqr), aes(x = Date, y = Log_Return)) +
  geom_line() +
  scale_x_date(
    breaks = seq(as.Date("2014-10-01"), as.Date("2023-01-01"), by = "1 years"),
    date_labels = "%Y"
  ) +
  labs(title = "Daily Squared Log Returns from 17/09/2014 - 20/01/2023",
       x = "Date",
       y = "Daily Log Returns") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# 3)
return_est <- log_returns[1:(length(log_returns)-20)]
y <- c(0,return_est)
n <- length(y)

negloglik <- function(params, y){
  alpha0 <- exp(params[1])
  alpha1 <- exp(params[2])
  alpha2 <- exp(params[3])
  sigma2 <- c(var(y), rep(0, n))
  loglik <- 0
  for (i in 1:n) {
    sigma2[i+1] <- alpha0 + alpha1*y[i]^2 + alpha2*sigma2[i]
    loglik <- loglik + (-1/2*log(2*pi) - 1/2*log(sigma2[i]) - 1/2*y[i]^2/sigma2[i])
  }
  return(-loglik)
}

initial_params <- c(0,0,0)
mle_c <- optim(initial_params, negloglik, y = y)
mle_c_results <- exp(mle_c$par)
mle_c_results
mle_c_llh <- -negloglik(mle_c$par, y)      # check the value of the log-likelihood in c)
mle_c_llh



# 4)
# qqnorm(return_est, main = "Normal Q-Q plot (Daily Log Returns)")
# qqline(return_est)

skewness<-function(y){
  n<-length(y)
  result<-(sum((y-mean(y))^3/n))/(sum((y-mean(y))^2/n)^(3/2))
  print(result)
}

kurtosis<-function(y){
  n<-length(y)
  result<-(sum((y-mean(y))^4/n))/(sum((y-mean(y))^2/n)^(2))
  print(result)
}

jb_stat <- (length(return_est)/6) * (skewness(return_est)^2 + (1/4) * (kurtosis(return_est)-3)^2)
jb_stat


# 5)
library("rugarch")

garch_e <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1,1)), 
                      mean.model = list(armaOrder = c(0, 0), include.mean = F), 
                      distribution.model = "norm")

garch_e_result <- ugarchfit(garch_e, data = return_est)
garch_e_result
garch_e_coef <- garch_e_result@fit$coef
garch_e_coef
garch_e_llh <- garch_e_result@fit$LLH      # check the value of the log-likelihood in e)
garch_e_llh   

res_sqr <- residuals(garch_e_result)^2
convar_fit <- sigma(garch_e_result)^2
plot(date[1:length(res_sqr)], res_sqr, type="l", xlab="Date", ylab = "Squared Residuals", main="Fitted Conditional Variances and The Squared Residuals" )
lines(date[1:length(res_sqr)], convar_fit, col="red")
legend("topleft", 
       legend=c("Squared Residuals", "Fitted Conditional Variances"), 
       col=c("black", "red"), 
       lty=1, 
       cex=0.8)


# 6)
sigma2 <- c(var(y), rep(0, n))
for (i in 1:n) {
  sigma2[i+1] <- mle_c_results[1] + mle_c_results[2]*y[i]^2 + mle_c_results[3]*sigma2[i]
}

sigma2_forecast <- c(sigma2[length(sigma2)], rep(0,20))
for (i in 1:20) {
  sigma2_forecast[i+1] <- mle_c_results[1] + (mle_c_results[2] + mle_c_results[3]) * sigma2_forecast[i]
}

plot(tail(date, 20), tail(log_returns_sqr, 20), lwd=1.5,type = "l", 
     xlab = "Date", ylab = "Squared Returns", 
     main = "Forecasted Conditional Variances and The Squared Returns")
lines(tail(date,20), tail(sigma2_forecast,20), lwd=1.5,col="red")
legend("topleft", 
       legend=c("The Actual Squared Returns", "Forecasted Conditional Variances"), 
       col=c("black", "red"), 
       lty=1, 
       cex=0.8)




