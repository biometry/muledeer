AllMeans <- read.csv("AllMeans.csv")
WholeAreaMeans <- read.csv("WholeAreaMeans.csv")




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
  error <- sum((abs(model-data))^2)
  return (data.frame(model, error))
}

result <- maf(par=c(N_0, lambda_gmean,tsteps),data = WholeAreaMeans$MDperKMsqFall_mean) #MAF:0,419
result <- mss(par=c(N_0, lambda_gmean,tsteps),data = WholeAreaMeans$MDperKMsqFall_mean) #MSS:16,59
plot(WholeAreaMeans$MDperKMsqFall_mean, cex=0.5)
lines(result$model, col="red") #does not fit the data well



###-------------Discrete logistic model

N_0 <- WholeAreaMeans$MDperKMsqFall_mean[1]
tsteps <- length(WholeAreaMeans$MDperKMsqFall_mean)

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
rd_maf <- result_maf[-1]/result_maf[-51]
opt_mss <- optim(par= c(rd=0.3, K=2), fn=mss, suppar=c(N0=N_0, steps=tsteps), data=WholeAreaMeans$MDperKMsqFall_mean)
result_mss <- log_d(c(opt_mss$par[1], opt_mss$par[2]),c(N0=N_0,  steps=tsteps)) #MSS:8.654571, rd:0.03246672, K:3.33934080
rd_mss <- result_mss[-1]/result_mss[-51]

plot(rd_maf~result_maf[-1], type="l")#correct:linear relationship of rd and N
cor(rd_maf,result_maf[-1])#-0.9999763
plot(rd_mss~result_mss[-1], type="l")#correct:linear
cor(rd_mss,result_mss[-1]) #-0.9999949

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
result_maf <- log_td(par=c(opt_maf$par[1], opt_maf$par[2], opt_maf$par[3]), c(N0=N_0,  steps=tsteps)) #MAF: 0.3045028, rd:0.1180750 K:2.1629352 theta 0.4040302
rd_maf <- result_maf[-1]/result_maf[-51]
opt_mss <- optim(par= c(rd=0.3, K=2, theta = 1), fn=mss, suppar=c(N0=N_0, steps=tsteps), data=WholeAreaMeans$MDperKMsqFall_mean)
result_mss <- log_td(c(opt_mss$par[1], opt_mss$par[2], opt_mss$par[3]),c(N0=N_0,  steps=tsteps)) #MSS:8.650606, rd:0.56623063 K:3.59976240 theta:0.03617367
rd_mss <- result_mss[-1]/result_mss[-51]


plot(WholeAreaMeans$MDperKMsqFall_mean~WholeAreaMeans$year, cex=0.5, main="Theta-model")
lines(result_maf~WholeAreaMeans$year, col="red") 
lines(result_mss~WholeAreaMeans$year, col="blue")
legend("topleft", legend=c("MAF", "MSS"), col=c("red", "blue"), lty=1)

plot(rd_maf~result_maf[-1], type="l")
cor(rd_maf,result_maf[-1])#-0.998955 (correct:less than model without theta)
plot(rd_mss~result_mss[-1], type="l")#correct: non-linear but cubix
cor(rd_mss,result_mss)#-0.99534


## Discrete logistic theta-model for each macrounit = on AllMeans
par(mfrow=c(2,2),oma=c(2,0,2,0))
macrounits <- levels(AllMeans$macrounit)
result <- data.frame(numeric(51),numeric(51),numeric(51),numeric(51))
rd <- data.frame(numeric(51),numeric(51),numeric(51),numeric(51))
parinfo <- data.frame(numeric(3), numeric(3), numeric(3), numeric(3))
for (i in 1:length(macrounits)){
  N0 <- AllMeans$MDperKMsqFall_mean[cond[1]]
  tsteps <- length(AllMeans$MDperKMsqFall_mean[cond])
  cond = which(AllMeans$macrounit==macrounits[i])  
  opt <- optim(par= c(rd=0.5, K=4, theta = 0.1), fn=mss, suppar=c(N0, tsteps), data=AllMeans$MDperKMsqFall_mean[cond])
  result[,i] <- log_td(c(opt$par[1], opt$par[2], opt$par[3]),c(N0=N_0,  steps=tsteps)) 
  rd[,i] <- c((result[-1,i]/result[-51,i]),NA)

  plot(AllMeans$MDperKMsqFall_mean[cond]~AllMeans$year[cond], cex=0.5, main=macrounits[i], xlab="Year",ylab="Population Density")
  lines(result[,i]~AllMeans$year[cond], col="red")
#   plot(rd[,i]~result[,i], type="l")
   parinfo[,i] <- c(opt$par)
}
title("Theta-Model - Population Densities", outer=TRUE)       
parinfo  


