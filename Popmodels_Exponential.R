

###----------- Discrete Exponential Growth with geometric mean lambda

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
  error <- sum((abs(model-data))^2, na.rm=TRUE)/length(data)
  return (data.frame(model, error))
}

result <- maf(par=c(N_0, lambda_gmean,tsteps),data = WholeAreaMeans$MDperKMsqFall_mean) #MAF:0,419
result <- mss(par=c(N_0, lambda_gmean,tsteps),data = WholeAreaMeans$MDperKMsqFall_mean) #MSS:16,59
plot(WholeAreaMeans$MDperKMsqFall_mean, cex=0.5)
lines(result$model, col="red") #does not fit the data well