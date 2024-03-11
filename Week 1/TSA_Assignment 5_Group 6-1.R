rm(list=ls())
library(ggplot2)

## Question a
IRX<-read.csv("^IRX.csv") 
FVX<-read.csv("^FVX.csv")
IRX[IRX=="null"]<-NA 
IRX<-na.omit(IRX)
FVX[FVX=="null"]<-NA 
FVX<-na.omit(FVX)
x_t <- as.numeric(IRX[,6])
y_t <- as.numeric(FVX[,6])

# Applying the formula in d) from Part 1:
Y <- rbind(t(x_t[2:length(x_t)]), t(y_t[2:length(y_t)]))
Z <- rbind(1, t(x_t[1:(length(x_t)-1)]), t(y_t[1:(length(x_t)-1)]))
B_hat <- Y %*% t(Z) %*% solve(Z %*% t(Z))
B_hat

# Applying OLS per equation, which yields the same results:
# Y_t_1<-cbind(1, x_t[1:(length(x_t)-1)], y_t[1:(length(x_t)-1)])
# Phi_r1 <- solve(t(Y_t_1) %*% Y_t_1) %*% t(Y_t_1) %*% x_t[2:length(x_t)]
# Phi_r2 <- solve(t(Y_t_1) %*% Y_t_1) %*% t(Y_t_1) %*% y_t[2:length(x_t)]
# Phi <- rbind(t(Phi_r1), t(Phi_r2))
# Phi



## Question b
# It can be observed that phi_12 and phi_21 are very close to zero.
# This implies that it cannot be assured that yt Granger causes xt or the other way around. 
# It will have to tested with the appropriate F-test.

Phi <- B_hat[,-1]
Phi_11 <- rep(0,10000)
Phi_12 <- rep(0,10000)
Phi_21 <- rep(0,10000)
Phi_22 <- rep(0,10000)

for (i in 1:10000) {
  if (i == 1) {
    Phi_11[i] <- Phi[1,1]
    Phi_12[i] <- Phi[1,2]
    Phi_21[i] <- Phi[2,1]
    Phi_22[i] <- Phi[2,2]
    Phi_old <- Phi
  } else{
      Phi_new <- Phi_old %*% Phi
      Phi_11[i] <- Phi_new[1,1]
      Phi_12[i] <- Phi_new[1,2]
      Phi_21[i] <- Phi_new[2,1]
      Phi_22[i] <- Phi_new[2,2]
      Phi_old <- Phi_new
    }
}

par(mfrow = c(2, 2))
plot(1:10000, Phi_11, type = "l", col = "red3", xlab = "Time Periods j", ylab = "IRF", main = expression("Impulse Response Function of " * x[t] * " on " * epsilon[1][t-j]))
plot(1:10000, Phi_12, type = "l", col = "red3", xlab = "Time Periods j", ylab = "IRF", main = expression("Impulse Response Function of " * x[t] * " on " * epsilon[2][t-j]))
plot(1:10000, Phi_21, type = "l", col = "red3", xlab = "Time Periods j", ylab = "IRF", main = expression("Impulse Response Function of " * y[t] * " on " * epsilon[1][t-j]))
plot(1:10000, Phi_22, type = "l", col = "red3", xlab = "Time Periods j", ylab = "IRF", main = expression("Impulse Response Function of " * y[t] * " on " * epsilon[2][t-j]))



## Question c
# Apply the formula in e) from Part 1
omega_hat <- (1 / (ncol(Y) - 3))  * Y %*% (diag(ncol(Y)) - t(Z) %*% solve(Z %*% t(Z)) %*% Z) %*% t(Y)
omega_hat

# Apply the formula in f) from Part 1

loglik <- function(Y, B, Z, cov){
  n <- ncol(Y)
  llh <- 0
  for (t in 1:n) {
    llh <- llh - log(det(2*pi*cov)) / 2 - t(Y[,t] - B %*% Z[,t]) %*% solve(cov) %*% (Y[,t] - B %*% Z[,t]) / 2
  }
  return(llh)
}

llh_var1 <- loglik(Y = Y, B = B_hat, Z = Z, cov = omega_hat)
print(llh_var1)



## Question d
AIC <- function(cov, p, m, n) {return (log(det(cov)) + 2 * p * m^2 / n) }
AIC_var <- rep(0,3)
llh_var <- rep(0,3)

# also include the parameters in the covariance matrix.
nr_param <- c(9, 13, 17)         

# VAR(1)
llh_var[1] <- llh_var1
AIC_var[1] <- AIC(cov = omega_hat, p = 1, m = 2, n = ncol(Y))

# VAR(2)
Y2 <- rbind(t(x_t[3:length(x_t)]), t(y_t[3:length(y_t)]))
Z2<- rbind(1, t(x_t[2:(length(x_t)-1)]), t(y_t[2:(length(x_t)-1)]), t(x_t[1:(length(x_t)-2)]), t(y_t[1:(length(x_t)-2)]))
B2_hat <- Y2 %*% t(Z2) %*% solve(Z2 %*% t(Z2))
B2_hat

omega2_hat <- (1 / (ncol(Y2) - 5))  * Y2 %*% (diag(ncol(Y2)) - t(Z2) %*% solve(Z2 %*% t(Z2)) %*% Z2) %*% t(Y2)
llh_var2 <- loglik(Y = Y2, B = B2_hat, Z = Z2, cov = omega2_hat)
llh_var[2] <- llh_var2
AIC_var[2] <- AIC(cov = omega2_hat, p = 2, m = 2, n = ncol(Y2))

# VAR(3)
Y3 <- rbind(t(x_t[4:length(x_t)]), t(y_t[4:length(y_t)]))
Z3<- rbind(1, t(x_t[3:(length(x_t)-1)]), t(y_t[3:(length(x_t)-1)]), t(x_t[2:(length(x_t)-2)]), t(y_t[2:(length(x_t)-2)]),
           t(x_t[1:(length(x_t)-3)]), t(y_t[1:(length(x_t)-3)]))
B3_hat <- Y3 %*% t(Z3) %*% solve(Z3 %*% t(Z3))
B3_hat

omega3_hat <- (1 / (ncol(Y3) - 7))  * Y3 %*% (diag(ncol(Y3)) - t(Z3) %*% solve(Z3 %*% t(Z3)) %*% Z3) %*% t(Y3)
llh_var3 <- loglik(Y = Y3, B = B3_hat, Z = Z3, cov = omega3_hat)
llh_var[3] <- llh_var3
AIC_var[3] <- AIC(cov = omega3_hat, p = 3, m = 2, n = ncol(Y3))

# results
cbind(AIC_var, llh_var, nr_param)



## Question e
# We prefer VAR(3) because it has the lowest AIC.



## Question f
# Johansen’s trace test for Cointegration on VAR(3)
LR <- rep(0,2)
Y <- cbind(x_t, y_t)
delta_Y <- diff(Y)
Y_t_1 <- Y[1:nrow(delta_Y),]
cancor_result <- cancor(delta_Y, Y_t_1)
lambda <- (cancor_result$cor)^2
LR[1] <- -(ncol(Y3)-3) * (log(1-lambda[1]) + log(1-lambda[2]))
LR[2] <- -(ncol(Y3)-3) * (log(1-lambda[2]))
print(LR)

# check the calculation of the Johansen-Statistics using package:
library(urca)
jotest <- ca.jo(data.frame(x_t, y_t), type = 'trace', K = 3, ecdet = 'none')
summary(jotest)





