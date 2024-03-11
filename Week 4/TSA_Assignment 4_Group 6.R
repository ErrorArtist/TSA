rm(list=ls())
library(ggplot2)

## Question a
rawdata <- read.csv("FEDFUNDS.csv")
data = rawdata
rate <- as.numeric(data$FEDFUNDS)
date <- as.Date(data$DATE, format = "%Y-%m-%d")
AIC <- rep(0,5)

# for each AR(p) model we use OLS and calculate the log-likelihood;
# to avoid lengthy codes, we define a function to perform OLS and calculate the AIC.

ols_llh_AIC <- function(X, y){
  n <- length(y)
  k <- ncol(X) + 1
  params <- solve(t(X) %*% X) %*% t(X) %*% y
  resid <- (y-X%*%params)
  l <- -(n/2) * log(2*pi*var(resid)) - (1/(2*var(resid))) * sum(resid^2)
  AIC <- -2 * l + 2 * k
  return(AIC)
}

# AR(1)
y_1 <- rate[1:(length(rate)-1)]
X <- cbind(1,y_1)
y <- rate[2:length(rate)]
AIC[1] <- ols_llh_AIC(X,y)

# AR(2)
y_1 <- rate[1:(length(rate)-1)]
y_2 <- rate[1:(length(rate)-2)]
y_1 <- y_1[-1]
X <- cbind(1,y_1, y_2)
y <- rate[3:length(rate)]
AIC[2] <- ols_llh_AIC(X,y)

# AR(3)
y_1 <- rate[1:(length(rate)-1)]
y_2 <- rate[1:(length(rate)-2)]
y_3 <- rate[1:(length(rate)-3)]
y_1 <- y_1[-(1:2)]
y_2 <- y_2[-1]
X <- cbind(1, y_1, y_2, y_3)
y <- rate[4:length(rate)]
AIC[3] <- ols_llh_AIC(X,y)

# AR(4)
y_1 <- rate[1:(length(rate)-1)]
y_2 <- rate[1:(length(rate)-2)]
y_3 <- rate[1:(length(rate)-3)]
y_4 <- rate[1:(length(rate)-4)]
y_1 <- y_1[-(1:3)]
y_2 <- y_2[-(1:2)]
y_3 <- y_3[-1]
X <- cbind(1, y_1, y_2, y_3, y_4)
y <- rate[5:length(rate)]
AIC[4] <- ols_llh_AIC(X,y)

# AR(5)
y_1 <- rate[1:(length(rate)-1)]
y_2 <- rate[1:(length(rate)-2)]
y_3 <- rate[1:(length(rate)-3)]
y_4 <- rate[1:(length(rate)-4)]
y_5 <- rate[1:(length(rate)-5)]
y_1 <- y_1[-(1:4)]
y_2 <- y_2[-(1:3)]
y_3 <- y_3[-(1:2)]
y_4 <- y_4[-1]
X <- cbind(1, y_1, y_2, y_3, y_4, y_5)
y <- rate[6:length(rate)]
AIC[5] <- ols_llh_AIC(X,y)
print(AIC)



## Question b
plot<-ggplot(data.frame(Date=date,Rate=rate),aes(x=Date,y=Rate))+geom_line()+geom_vline(xintercept = as.numeric(as.Date("1970-01-01")), linetype = "dashed", color = "red") +
  geom_vline(xintercept = as.numeric(as.Date("1977-01-01")), linetype = "dashed", color = "red") +
  geom_vline(xintercept = as.numeric(as.Date("1981-01-01")), linetype = "dashed", color = "red") +
  geom_vline(xintercept = as.numeric(as.Date("1990-01-01")), linetype = "dashed", color = "red") + labs(title = "Monthly FFER from 01/07/1954 - 31/01/2024 with highlited dates")+
  geom_vline(xintercept = as.numeric(as.Date("2000-01-01")), linetype = "dashed", color = "red") 
plot<-plot+scale_x_date(breaks = as.Date(c("1970-01-01","1977-01-01","1981-01-01","1990-01-01","2000-01-01")),date_labels = "%Y")
plot+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
           panel.background = element_blank(), axis.line = element_line(colour = "black"))



## Question c --- Calculate the AIC manually

# Define a function TARnegllh which returns the negative Log-Likelihood of the TAR model.
TARnegllh <- function(X, y, tau, params, nr_lag){
  n <- length(y)
  tau <- tau - nr_lag    # we lose "nr_lag" observations in the AR model, so the index of threshold should shift forward accordingly.
  P <- matrix(params, nrow = nrow(X), ncol = length(params), byrow = T)
  for (i in 1:length(tau)){
    P[1:(tau[i] - 1), ((i*(nr_lag + 1))+1):((nr_lag + 1)*(i + 1))] <- 0
  }
  X <- matrix(rep(X, (length(tau)+1)), nrow = nrow(X))
  resid <- (y-rowSums(X*P))
  llh <- -(n/2) * log(2*pi*var(resid)) - (1/(2*var(resid))) * sum(resid^2)
  return(-llh)
}

# i)
tau = which(date ==  "1981-01-01")
nr_lag <- 1
params <- rep(0, (1+nr_lag)*(length(tau)+1))
y_1 <- rate[1:(length(rate)-1)]
X <- cbind(1,y_1)
y <- rate[2:length(rate)]
llh_tar1 <- optim(params, TARnegllh, method ="BFGS", y=y, X=X, tau=tau, nr_lag=nr_lag)$value
AIC_tar1 <- -2 * (-llh_tar1) + 2 * (length(params)+1)

nr_lag <- 2
params = rep(0, (1+nr_lag)*(length(tau)+1))
y_1 <- rate[1:(length(rate)-1)]
y_2 <- rate[1:(length(rate)-2)]
y_1 <- y_1[-1]
X <- cbind(1,y_1, y_2)
y <- rate[3:length(rate)]
llh_tar2 <- optim(params, TARnegllh, method ="BFGS", y=y, X=X, tau=tau, nr_lag=nr_lag)$value
AIC_tar2 <- -2 * (-llh_tar2) + 2 * (length(params)+1)

# ii)
tau = c(which(date == "1970-01-01"), which(date == "1990-01-01"))
nr_lag <- 1
params <- rep(0, (1+nr_lag)*(length(tau)+1))
y_1 <- rate[1:(length(rate)-1)]
X <- cbind(1,y_1)
y <- rate[2:length(rate)]
llh_tar3 <- optim(params, TARnegllh, method ="BFGS", y=y, X=X, tau=tau, nr_lag=nr_lag)$value
AIC_tar3 <- -2 * (-llh_tar3) + 2 * (length(params)+1)

nr_lag <- 2
params = rep(0, (1+nr_lag)*(length(tau)+1))
y_1 <- rate[1:(length(rate)-1)]
y_2 <- rate[1:(length(rate)-2)]
y_1 <- y_1[-1]
X <- cbind(1,y_1, y_2)
y <- rate[3:length(rate)]
llh_tar4 <- optim(params, TARnegllh, method ="BFGS", y=y, X=X, tau=tau, nr_lag=nr_lag)$value
AIC_tar4 <- -2 * (-llh_tar4) + 2 * (length(params)+1)

# iii)
tau = c(which(date == "1970-01-01"), which(date == "1977-01-01"),
        which(date == "1981-01-01"), which(date == "1990-01-01"),
        which(date == "2000-01-01"))
nr_lag <- 1
params <- rep(0, (1+nr_lag)*(length(tau)+1))
y_1 <- rate[1:(length(rate)-1)]
X <- cbind(1,y_1)
y <- rate[2:length(rate)]
llh_tar5 <- optim(params, TARnegllh, method ="BFGS", y=y, X=X, tau=tau, nr_lag=nr_lag)$value
AIC_tar5 <- -2 * (-llh_tar5) + 2 * (length(params)+1)

nr_lag <- 2
params = rep(0, (1+nr_lag)*(length(tau)+1))
y_1 <- rate[1:(length(rate)-1)]
y_2 <- rate[1:(length(rate)-2)]
y_1 <- y_1[-1]
X <- cbind(1,y_1, y_2)
y <- rate[3:length(rate)]
llh_tar6 <- optim(params, TARnegllh, method ="BFGS", y=y, X=X, tau=tau, nr_lag=nr_lag)$value
AIC_tar6 <- -2 * (-llh_tar6) + 2 * (length(params)+1)

# results
AIC_TAR <- matrix(c(AIC_tar1, AIC_tar2, AIC_tar3, AIC_tar4, AIC_tar5, AIC_tar6), nrow = 3, ncol = 2, byrow = T)
rownames(AIC_TAR) <- c("(i)", "(ii)", "(iii)")
colnames(AIC_TAR) <- c("p=1", "p=2")
AIC_TAR



## Question c --- Check the results using the built-in lm function
# i)
n <- nrow(data)
data$y_t_1 <- c(NA, data$FEDFUNDS[-n])
data$y_t_2 <- c(rep(NA,2), data$FEDFUNDS[-c((n-1),n)])

TAR1 <- lm(FEDFUNDS ~ y_t_1 + y_t_1 * (DATE > "1981-01-01"), data)
TAR2 <- lm(FEDFUNDS ~ y_t_1 + y_t_2 + y_t_1 * (DATE > "1981-01-01") + y_t_2 * (DATE > "1981-01-01"), data)
c(AIC(TAR1), AIC(TAR2))

# ii)
TAR3 <- lm(FEDFUNDS ~ y_t_1 + y_t_1 * (DATE > "1970-01-01") 
           + y_t_1 * (DATE > "1990-01-01"), data)
TAR4 <- lm(FEDFUNDS ~ y_t_1 + y_t_2 + y_t_1 * (DATE > "1970-01-01") 
           + y_t_1 * (DATE > "1990-01-01") + y_t_2 * (DATE > "1970-01-01" ) + y_t_2 * (DATE > "1990-01-01"), data)
c(AIC(TAR3), AIC(TAR4))

# iii)
TAR5 <- lm(FEDFUNDS ~ y_t_1 + y_t_1 * (DATE > "1970-01-01") + y_t_1 * (DATE > "1977-01-01") 
           + y_t_1 * (DATE > "1981-01-01") + y_t_1 * (DATE > "1990-01-01") + y_t_1 * (DATE > "2000-01-01"), data)

TAR6 <- lm(FEDFUNDS ~ y_t_1 + y_t_2 + y_t_1 * (DATE > "1970-01-01")
           + y_t_1 * (DATE > "1977-01-01") + y_t_1 * (DATE > "1981-01-01" ) + y_t_1 * (DATE > "1990-01-01") + y_t_1 * (DATE > "2000-01-01") 
           + y_t_2 * (DATE > "1970-01-01") + y_t_2 * (DATE > "1977-01-01") + y_t_2 * (DATE > "1981-01-01") + y_t_2 * (DATE > "1990-01-01") + y_t_2 * (DATE > "2000-01-01") , data)

c(AIC(TAR5), AIC(TAR6))

# results
tar_i <- c(AIC(TAR1), AIC(TAR2)) 
tar_ii <- c(AIC(TAR3), AIC(TAR4))
tar_iii <- c(AIC(TAR5), AIC(TAR6))
aic_tar <- matrix(cbind(tar_i , tar_ii , tar_iii), ncol = 3, nrow = 2, byrow = FALSE)
colnames(aic_tar) <- c("(i)", "(ii)", "(ii)") 
rownames(aic_tar) <- c("p=1", "p=2")
aic_tar



## Question d
library(lubridate) 
gamma <- c(0.5, 1, 2, 5, 10)
data$diff <-interval("1981-01-01", data$DATE) %/% months(1)

aic_star1 <- rep(NA, length(gamma))
aic_star2 <- rep(NA, length(gamma)) 

j <- 1
for(i in gamma){
  data$sf <- 1/(1 + exp(-i*(data$diff)))
  STAR1 = lm(FEDFUNDS ~ y_t_1 + sf*(1 + y_t_1), data)
  STAR2 = lm(FEDFUNDS ~ y_t_1 + y_t_2 + sf*(1 + y_t_1 + y_t_2), data)
  aic_star1[j] <- AIC(STAR1) 
  aic_star2[j] <- AIC(STAR2)
  j <- j + 1 
}
aic_star <- matrix(cbind(gamma , aic_star1 , aic_star2), ncol = 3, nrow = 5, byrow = F)
colnames(aic_star) <- c("Gamma", "AIC Star(1)", "AIC Star(2)")
aic_star

