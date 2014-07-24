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
  error <- sum((abs(model-data))^2, na.rm=TRUE)/length(data)
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
  error <- sum((abs(model-data))^2, na.rm=TRUE)/length(data)
  return (error)
}

opt_maf <- optim(par= c(rd=0.3, K=2), fn=maf, suppar=c(N0=N_0, steps=tsteps), data=WholeAreaMeans$MDperKMsqFall_mean)
result_maf <- log_d(c(opt_maf$par[1], opt_maf$par[2]),c(N0=N_0,  steps=tsteps)) #MAF: 0.3056459, rd:0.05669221, K:2.12359422
rd_maf <- c(((result_maf[-1]-result_maf[-tsteps])/result_maf[-tsteps]),NA)
opt_mss <- optim(par= c(rd=0.3, K=2), fn=mss, suppar=c(N0=N_0, steps=tsteps), data=WholeAreaMeans$MDperKMsqFall_mean)
result_mss <- log_d(c(opt_mss$par[1], opt_mss$par[2]),c(N0=N_0,  steps=tsteps)) #MSS:8.654571, rd:0.03246672, K:3.33934080
rd_mss <- c(((result_mss[-1]-result_mss[-tsteps])/result_mss[-tsteps]),NA)



plot(rd_maf~result_maf, type="l", col="red")#correct:linear relationship of rd and N
cor(rd_maf,result_maf, use="pairwise.complete.obs")#--1
lines(rd_mss~result_mss, type="l", col="blue")#correct:linear
cor(rd_mss,result_mss, use="pairwise.complete.obs") #-1

plot(WholeAreaMeans$MDperKMsqFall_mean, cex=0.5, main="Discrete logistic model")
lines(result_maf, col="red") 
lines(result_mss, col="blue")
legend("topleft", legend=c("MAF", "MSS"), col=c("red", "blue"), lty=1)



### Discrete logistic model for each macrounit = on AllMeans
Popdata <- data.frame(MDperKMsqFall_mean=AllMeans$MDperKMsqFall_mean, year=AllMeans$year, macrounit=AllMeans$macrounit, RepRateFall_mean =AllMeans$RepRateFall_mean)#remove NAs
Popdata <- Popdata[which(complete.cases(Popdata)==TRUE),]

par(mfrow=c(2,4),oma=c(2,0,2,0))
macrounits <- levels(Popdata$macrounit)
result_out <- list()
rd_out <- list()
parinfo <- data.frame("0-1" = numeric(3), "0-2" = numeric(3), "0-3" = numeric(3), "0-4" = numeric(3), row.names=c("rd", "K", "MSS"))


for (i in 1:length(macrounits)){
  cond = which(Popdata$macrounit==macrounits[i])  
  
  N_0 <- Popdata$MDperKMsqFall_mean[cond[1]]
  tsteps <- length(Popdata$MDperKMsqFall_mean[cond])
  opt <- optim(par= c(rd=0.3, K=2), fn=maf, method="L-BFGS-B", lower=c(0 , 0), suppar=c(N_0, tsteps), data=Popdata$MDperKMsqFall_mean[cond])
  result <- log_d(c(opt$par[1], opt$par[2]),c(N0=N_0,  steps=tsteps))
  names(result) <- Popdata$year[cond]
  result_out <- append(result_out, list(c(result))) 
  names(result) <- Popdata$year[cond]
  rd <- c(((result[-1]-result[-tsteps])/result[-tsteps]),NA)
  rd_out <- append(rd_out, list(c(rd)))
  parinfo[,i] <- c(opt$par, opt$value)
  
  
  plot(Popdata$MDperKMsqFall_mean[cond]~Popdata$year[cond], cex=0.5, main=macrounits[i], xlab="Year",ylab="Population Density")
  lines(result~Popdata$year[cond], col="red")
  plot(rd~result, type="l", lty=2, main=macrounits[i], xlab="Predicted Population Density",ylab="Predicted Reproduction rate")
  remove(result)
  remove(rd)
}
title("Basic Discrete Logistic Model - Population Densities and Reproduction Rates (MAF)", outer=TRUE)       
parinfo  




### ----------------- Discrete logistic Theta-Model

N_0 <- WholeAreaMeans$MDperKMsqFall_mean[1]
tsteps <- length(WholeAreaMeans$MDperKMsqFall_mean)

log_td <- function(par, suppar){
  output <- numeric(suppar[2])#steps
  output[1] <- suppar[1]#n0
  for (t in 1:(length(output)-1)){
    output[t+1] <- output[t]+(output[t]*par[1]*(1-((output[t]/par[2])^par[3])))
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
  error <- sum((abs(model-data))^2, na.rm=TRUE)/length(data)
  return (error)
}


opt_maf <- optim(par= c(rd=0.3, K=2, theta = 1), fn=maf, suppar=c(N0=N_0, steps=tsteps), data=WholeAreaMeans$MDperKMsqFall_mean)
result_maf <- log_td(par=c(opt_maf$par[1], opt_maf$par[2], opt_maf$par[3]), c(N0=N_0,  steps=tsteps)) #MAF: 0.3045028, rd:0.1180750 K:2.1629352 theta 0.4040302
rd_maf <- c(((result_maf[-1]-result_maf[-tsteps])/result_maf[-tsteps]),NA)
opt_mss <- optim(par= c(rd=0.3, K=2, theta = 1), fn=mss, suppar=c(N0=N_0, steps=tsteps), data=WholeAreaMeans$MDperKMsqFall_mean)
result_mss <- log_td(c(opt_mss$par[1], opt_mss$par[2], opt_mss$par[3]),c(N0=N_0,  steps=tsteps)) #MSS:8.650606, rd:0.56623063 K:3.59976240 theta:0.03617367
rd_mss <- c(((result_mss[-1]-result_mss[-tsteps])/result_mss[-tsteps]),NA)


plot(WholeAreaMeans$MDperKMsqFall_mean~WholeAreaMeans$year, cex=0.5, main="Theta-model")
lines(result_maf~WholeAreaMeans$year, col="red") 
lines(result_mss~WholeAreaMeans$year, col="blue")
legend("topleft", legend=c("MAF", "MSS"), col=c("red", "blue"), lty=1)

plot(rd_maf~result_maf, type="l", col="red")#concave shape,correct
cor(rd_maf,result_maf, use="pairwise.complete.obs")#-0.999
lines(rd_mss~result_mss, type="l", col="blue")#more concave
cor(rd_mss,result_mss, use="pairwise.complete.obs") #-0.995





###----------------------Discrete logistic theta-model for each macrounit = on AllMeans

Popdata <- data.frame(MDperKMsqFall_mean=AllMeans$MDperKMsqFall_mean, year=AllMeans$year, macrounit=AllMeans$macrounit, RepRateFall_mean =AllMeans$RepRateFall_mean)#remove NAs
Popdata <- Popdata[which(complete.cases(Popdata)==TRUE),]

par(mfrow=c(2,4),oma=c(2,0,2,0))
macrounits <- levels(Popdata$macrounit)
result_out <- list()
rd_out <- list()
parinfo <- data.frame("0-1" = numeric(4), "0-2" = numeric(4), "0-3" = numeric(4), "0-4" = numeric(4), row.names=c("rd", "K", "theta", "MSS"))


for (i in 1:length(macrounits)){
  cond = which(Popdata$macrounit==macrounits[i])  
  
  N_0 <- Popdata$MDperKMsqFall_mean[cond[1]]
  tsteps <- length(Popdata$MDperKMsqFall_mean[cond])
  opt <- optim(par= c(rd=0.3, K=2, theta = 1), fn=maf, method="L-BFGS-B", lower=c(0, 0 , 0), suppar=c(N_0, tsteps), data=Popdata$MDperKMsqFall_mean[cond])
  result <- log_td(c(opt$par[1], opt$par[2], opt$par[3]),c(N0=N_0,  steps=tsteps))
  names(result) <- Popdata$year[cond]
  result_out <- append(result_out, list(c(result))) 
  names(result) <- Popdata$year[cond]
  rd <- c(((result[-1]-result[-tsteps])/result[-tsteps]),NA)
  rd_out <- append(rd_out, list(c(rd)))
  parinfo[,i] <- c(opt$par, opt$value)
  
  
  plot(Popdata$MDperKMsqFall_mean[cond]~Popdata$year[cond], cex=0.5, main=macrounits[i], xlab="Year",ylab="Population Density")
  lines(result~Popdata$year[cond], col="red")
  plot(rd~result, type="l", lty=2, main=macrounits[i], xlab="Predicted Population Density",ylab="Predicted Reproduction rate", col="red")
  remove(result)
  remove(rd)
}
title("Discrete Logistic Theta-Model - Population Densities and Reproduction Rates (MAF)", outer=TRUE)       
parinfo  


# visualization of density dependence observed vs. modelled
par(mfrow=c(2,2))
for (i in 1:length(macrounits)){
  cond = which(Popdata$macrounit==macrounits[i])  
  plot(Popdata$RepRateFall_mean[cond]~Popdata$MDperKMsqFall_mean[cond], main=macrounits[i], type="p",cex=0.5,xlab="Observed Population Density",ylab="Observed and Predicted Reproduction rate")
  curve ((parinfo[1,i]*(1-(x/parinfo[2,i])^parinfo[3,i])), from=min(Popdata$MDperKMsqFall_mean[cond]), to=max(Popdata$MDperKMsqFall_mean[cond]), n=1000, col="red", add=TRUE)
  #lines(rd_out[[i]]~result_out[[i]], type="l", lty=2, col="red")
}
legend("topright", legend=c("Observed", "Predicted"), col=c("black", "red"), lty=c(1,2))
title("Observed and Predicted Density Dependences Theta Model MAF", outer=TRUE)       

#recalculation using nls
nls <- nls(RepRateFall_mean[cond] ~ rd*MDperKMsqFall_mean[cond]*((1-MDperKMsqFall_mean[cond]/K)^theta), data=Popdata, lower=c(0,0,0), start=c(rd=8.66, K=1.49, theta = 0.023))
#error "missing values or infitiy produced when evaluating model"
#NAs?
count_nas(Popdata$MDperKMsqFall_mean)#0
count_nas(Popdata$RepRateFall_mean)#0


### -----------Including Harvest in Discrete logistic theta-model for each macrounit = on AllMeans
Popdata <- data.frame(MDperKMsqFall_mean=AllMeans$MDperKMsqFall_mean, year=AllMeans$year, macrounit=AllMeans$macrounit, HuntDen_All_mean=AllMeans$HuntDen_All_mean)#remove NAs
Popdata <- Popdata[which(complete.cases(Popdata)==TRUE),]
##huntdenAll no NAs

log_tdh <- function(par, suppar){
  output <- numeric(suppar[2])#steps
  output[1] <- suppar[1]#n0
  for (t in 1:(length(output)-1)){
    output[t+1] <- output[t]+(output[t]*par[1]*(1-((output[t]/par[2])^par[3])))-(suppar[t+2]/2)#from suppar 3 on: harvest data
  }
  return(output)
}

maf_tdh <- function(par, data, suppar){
  model <- log_tdh(par, suppar)
  error <- mean(abs(model-data))
  return (error)
}

mss_tdh <- function(par, data, suppar){
  model <- log_tdh(par, suppar)
  error <- sum((abs(model-data))^2, na.rm=TRUE)/length(data)
  return (error)
}


par(mfrow=c(2,4),oma=c(2,0,2,0))
macrounits <- levels(Popdata$macrounit)
result_out <- list()
rd_out <- list()
parinfo <- data.frame("0-1" = numeric(4), "0-2" = numeric(4), "0-3" = numeric(4), "0-4" = numeric(4), row.names=c("rd", "K", "theta", "MSS"))


for (i in 1:length(macrounits)){
  cond = which(Popdata$macrounit==macrounits[i])  
  
  N_0 <- Popdata$MDperKMsqFall_mean[cond[1]]
  tsteps <- length(Popdata$MDperKMsqFall_mean[cond])
  opt <- optim(par= c(rd=0.3, K=2, theta = 1), fn=mss, method="L-BFGS-B", lower=c(NA,0 , 0), suppar=c(N_0, tsteps, Popdata$HuntDen_All_mean[cond]), data=Popdata$MDperKMsqFall_mean[cond])
  result <- log_tdh(par=c(opt$par[1], opt$par[2], opt$par[3]),suppar=c(N0=N_0,  steps=tsteps,  Popdata$HuntDen_All_mean[cond]))
  names(result) <- Popdata$year[cond]
  result_out <- append(result_out, list(c(result))) 
  names(result) <- Popdata$year[cond]
  rd <- c(((result[-1]-result[-tsteps])/result[-tsteps]),NA)
  rd_out <- append(rd_out, list(c(rd)))
  parinfo[,i] <- c(opt$par, opt$value)
  
  
  plot(Popdata$MDperKMsqFall_mean[cond]~Popdata$year[cond], cex=0.5, main=macrounits[i], xlab="Year",ylab="Population Density")
  lines(result~Popdata$year[cond], col="red")
  plot(rd~result, type="l", lty=2, main=macrounits[i], xlab="Predicted Population Density",ylab="Predicted Reproduction rate")
  remove(result)
  remove(rd)
}
title("Discrete Logistic Theta-Model including Harvest- Population Densities and Reproduction Rates (MAF)", outer=TRUE)       
parinfo  

result_out
