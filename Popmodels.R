AllMeans <- read.csv("AllMeans.csv")
WholeAreaMeans <- read.csv("WholeAreaMeans.csv")


lambda1 <- AllMeans$MDperKMsqFall_mean_tplus1/AllMeans$MDperKMsqFall_mean 
lambda_gmean <- prod(lambda1, na.rm = TRUE)^(1/length(lambda1))
N_0 <- WholeAreaMeans$MDperKMsqFall_mean[1]
tsteps <- length(WholeAreaMeans$MDperKMsqFall_mean)



###----------- Discrete Exponential Growth with geometric mean lambda


exp_discrete <- function(N0, lambda, steps){
  output <- numeric(steps)
  output[1] <- N0
  for(t in 1:(steps-1)){
    output[t+1] <- output[t]*lambda
  }
  return(output)
}


maf <- function(par, data){
  model <- exp_discrete(par[1],par[2], par[3])
  error <- mean(abs(model-data))
  return (data.frame(model, error))
}

mss <- function(par, data){
  model <- exp_discrete(par[1],par[2], par[3])
  error <- sum((abs(model-data))^2)
  return (data.frame(model, error))
}

result <- maf(par=c(N_0, lambda_gmean,tsteps),data = WholeAreaMeans$MDperKMsqFall_mean) #MAF:0,419
result <- mss(par=c(N_0, lambda_gmean,tsteps),data = WholeAreaMeans$MDperKMsqFall_mean) #MSS:16,59
plot(WholeAreaMeans$MDperKMsqFall_mean, cex=0.5)
lines(result$model, col="red") #does not fit the data well



###-------------Discrete logistic model

log_d <- function(par, suppar){
  output <- numeric(suppar[2])
  output[1] <- suppar[1]
  for (t in 1:(length(output)-1)){
    output[t+1] <- output[t]+output[t]*par[1]*(1-output[t]/par[2])
  }
  return(output)
}

maf <- function(par, data, suppar){
  model <- log_d(par, suppar)
  error <- mean(abs(model-data))
  return (error)
}

mss <- function(par, data, suppar){
  model <- log_d(par, suppar)
  error <- sum((abs(model-data))^2)
  return (error)
}

opt_maf <- optim(par= c(rd=0.3, K=2), fn=maf, suppar=c(N0=N_0, steps=tsteps), data=WholeAreaMeans$MDperKMsqFall_mean)
result_maf <- log_d(c(opt_maf$par[1], opt_maf$par[2]),c(N0=N_0,  steps=tsteps)) #MAF: 0.3056459, rd:0.05669221, K:2.12359422
opt_mss <- optim(par= c(rd=0.3, K=2), fn=mss, suppar=c(N0=N_0, steps=tsteps), data=WholeAreaMeans$MDperKMsqFall_mean)
result_mss <- log_d(c(opt_mss$par[1], opt_mss$par[2]),c(N0=N_0,  steps=tsteps)) #MSS:8.654571, rd:0.03246672, K:3.33934080


plot(WholeAreaMeans$MDperKMsqFall_mean, cex=0.5, main="Discrete logistic model")
lines(result_maf, col="red") 
lines(result_mss, col="blue")
legend("topleft", legend=c("MAF", "MSS"), col=c("red", "blue"), lty=1)


### ----------------- Discrete logistic Theta-Model


log_td <- function(par, suppar){
  output <- numeric(suppar[2])
  output[1] <- suppar[1]
  for (t in 1:(length(output)-1)){
    output[t+1] <- output[t]+output[t]*par[1]*(1-((output[t]/par[2])^par[3]))
  }
  return(output)
}

maf <- function(par, data, suppar){
  model <- log_td(par, suppar)
  error <- mean(abs(model-data))
  return (error)
}

mss <- function(par, data, suppar){
  model <- log_td(par, suppar)
  error <- sum((abs(model-data))^2)
  return (error)
}


opt_maf <- optim(par= c(rd=0.3, K=2, theta = 1), fn=maf, suppar=c(N0=N_0, steps=tsteps), data=WholeAreaMeans$MDperKMsqFall_mean)
result_maf <- log_td(c(opt_maf$par[1], opt_maf$par[2], opt_maf$par[3]), c(N0=N_0,  steps=tsteps)) #MAF: 0.3045028, rd:0.1180750 K:2.1629352 theta 0.4040302
opt_mss <- optim(par= c(rd=0.3, K=2, theta = 1), fn=mss, suppar=c(N0=N_0, steps=tsteps), data=WholeAreaMeans$MDperKMsqFall_mean)
result_mss <- log_td(c(opt_mss$par[1], opt_mss$par[2], opt_mss$par[3]),c(N0=N_0,  steps=tsteps)) #MSS:8.650606, rd:0.56623063 K:3.59976240 theta:0.03617367


plot(WholeAreaMeans$MDperKMsqFall_mean~WholeAreaMeans$year, cex=0.5, main="Theta-model")
lines(result_maf~WholeAreaMeans$year, col="red") 
lines(result_mss~WholeAreaMeans$year, col="blue")
legend("topleft", legend=c("MAF", "MSS"), col=c("red", "blue"), lty=1)
