AllMeans <- read.csv("AllMeans.csv")
WholeAreaMeans <- read.csv("WholeAreaMeans.csv")



#----------- Discrete Exponential Growth with geometric mean lambda

lambda1 <- AllMeans$MDperKMsqFall_mean_tplus1/AllMeans$MDperKMsqFall_mean 
lambda_gmean <- prod(lambda1, na.rm = TRUE)^(1/length(lambda1))
N_0 <- WholeAreaMeans$MDperKMsqFall_mean[1]
tsteps <- length(WholeAreaMeans$MDperKMsqFall_mean)

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



#-------------Discrete logistic model

log_discrete <- function(rd, K, N0, steps){
  output <- numeric(steps)
  output[1] <- N0
  for (t in 1:(length(output)-1)){
    output[t+1] <- output[t]+output[t]*rd*(1-output[t]/K)
  }
  return(output)
}

mss <- function(par, data, parfix){
  model <- log_discrete(par[1], par[2], parfix[1], parfix[2])
  error <- sum((abs(model-data))^2)
  return (error)
}

opt <- optim(par= c(rd=0.3, K=2), fn=mss, parfix=c(N0=N_0,  steps=tsteps), data=WholeAreaMeans$MDperKMsqFall_mean)
result <- log_discrete(opt$par[1], opt$par[2],N0=N_0,  steps=tsteps) #MSS:8.654571, rd:0.03246672, K:3.33934080

plot(WholeAreaMeans$MDperKMsqFall_mean, cex=0.5)
lines(result, col="red") 
