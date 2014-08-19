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



### Discrete logistic model for each macrounit = on AllMeans
Popdata <- data.frame(MDperKMsqFall_mean=AllMeans$MDperKMsqFall_mean, year=AllMeans$year, macrounit=AllMeans$macrounit, RepRateFall_mean =AllMeans$RepRateFall_mean)#remove NAs
Popdata <- Popdata[which(complete.cases(Popdata)==TRUE),]


parinfolist <- list()
for (i in 1:50){

#png("Basic1_MSS.png", width=2200, height=1500)
par(mfrow=c(2,2),oma=c(10,10,10,10),mar=c(10, 10, 10, 10))
macrounits <- levels(Popdata$macrounit)
result_out <- list()
rd_out <- list()
parinfo <- data.frame("0-1" = numeric(3), "0-2" = numeric(3), "0-3" = numeric(3), "0-4" = numeric(3), row.names=c("rd", "K", "MSS"))


for (i in 1:length(macrounits)){
  cond = which(Popdata$macrounit==macrounits[i])  
  
  N_0 <- Popdata$MDperKMsqFall_mean[cond[1]]
  tsteps <- length(Popdata$MDperKMsqFall_mean[cond])
  opt <- optim(par= c(rd=0.3, K=2), fn=mss, method="L-BFGS-B", lower=c(0 , 0), suppar=c(N_0, tsteps), data=Popdata$MDperKMsqFall_mean[cond])
  result <- log_d(c(opt$par[1], opt$par[2]),c(N0=N_0,  steps=tsteps))
  names(result) <- Popdata$year[cond]
  result_out <- append(result_out, list(c(result))) 
  names(result) <- Popdata$year[cond]
  rd <- c(((result[-1]-result[-tsteps])/result[-tsteps]),NA)
  rd_out <- append(rd_out, list(c(rd)))
  parinfo[,i] <- c(opt$par, opt$value)
  
  
  #plot(Popdata$MDperKMsqFall_mean[cond]~Popdata$year[cond], xlab="",ylab="", cex.axis=2.5, cex.main = 3, main=macrounits[i])
  #lines(result~Popdata$year[cond], col="red")
  #plot(rd~result, type="l", lty=2, main=macrounits[i], xlab="Predicted Population Density",ylab="Predicted Reproduction rate")
  remove(result)
  remove(rd)
}
#mtext("Year", 1, 1, outer=TRUE, cex=3)
#mtext("Population Density", 2, 1, outer=TRUE, las=0, cex=3)
#mtext("Basic1", 3, 1, outer=TRUE, cex=3.5)
#parinfo  
#dev.off()
parinfolist <- append(parinfolist, as.list(parinfo))
}
w <- which(names(parinfolist) == "X0.4")
parinfolist[w]
hist(parinfolist[[which(names(parinfolist) == "X0.4")]])


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





# Including Harvest in basic logistic model for each macrounit = on AllMeans-----


Popdata <- data.frame(MDperKMsqFall_mean=AllMeans$MDperKMsqFall_mean, year=AllMeans$year, macrounit=AllMeans$macrounit, HuntDen_All_mean=AllMeans$HuntDen_All_mean)#remove NAs
Popdata <- Popdata[which(complete.cases(Popdata)==TRUE),]
##huntdenAll no NAs

log_dh <- function(par, suppar){
  output <- numeric(suppar[2])#steps
  output[1] <- suppar[1]#n0
  for (t in 1:(length(output)-1)){
    output[t+1] <- output[t]+(output[t]*par[1]*(1-(output[t]/par[2])))-(par[3]*suppar[t+2])
    #from suppar 3 on: harvest data, weight term is fitted to avoid negative values if HuntDen<PopDen  
  }
  return(output)
}

maf_dh <- function(par, data, suppar){
  model <- log_dh(par, suppar)
  error <- mean(abs(model-data))
  if (length(which(model < 0)) != 0){error <- length(which(model < 0))}# penalty: the more negative values, the worse MSS
  return (error)
}

mss_dh <- function(par, data, suppar){
  model <- log_dh(par, suppar)
  error <- sum((abs(model-data))^2, na.rm=TRUE)/length(data)
  if (length(which(model < 0)) != 0){error <- length(which(model < 0))}# penalty: the more negative values, the worse MSS
  return (error)
}

png("Basic2_MSS_weight0,15.png", width=2200, height=1500)
par(mfrow=c(2,2),oma=c(10,10,10,10),mar=c(10, 10, 10, 10))
macrounits <- levels(Popdata$macrounit)
result_out <- list()
rd_out <- list()
parinfo <- data.frame("0-1" = numeric(4), "0-2" = numeric(4), "0-3" = numeric(4), "0-4" = numeric(4), row.names=c("rd", "K", "weight of harvest density", "MSS"))


for (i in 1:length(macrounits)){
  cond = which(Popdata$macrounit==macrounits[i])  
  
  N_0 <- Popdata$MDperKMsqFall_mean[cond[1]]
  tsteps <- length(Popdata$MDperKMsqFall_mean[cond])
  opt <- optim(par= c(rd=0.3, K=2, weight=0.15), fn=mss_dh, method="L-BFGS-B", lower=c(0, 0 , 0, 0), suppar=c(N_0, tsteps, Popdata$HuntDen_All_mean[cond]), data=Popdata$MDperKMsqFall_mean[cond])
  result <- log_dh(par=c(opt$par[1], opt$par[2], opt$par[3]),suppar=c(N0=N_0,  steps=tsteps,  Popdata$HuntDen_All_mean[cond]))
  names(result) <- Popdata$year[cond]
  result_out <- append(result_out, list(c(result))) 
  names(result) <- Popdata$year[cond]
  rd <- c(((result[-1]-result[-tsteps])/result[-tsteps]),NA)
  rd_out <- append(rd_out, list(c(rd)))
  parinfo[,i] <- c(opt$par, opt$value)
  
  
  plot(Popdata$MDperKMsqFall_mean[cond]~Popdata$year[cond], xlab="",ylab="", cex.axis=2.5, cex.main = 3, main=macrounits[i])
  lines(result~Popdata$year[cond], col="red")
    #plot(rd~result, type="l", lty=2, main=macrounits[i], xlab="Predicted Population Density",ylab="Predicted Reproduction rate")
  remove(result)
  remove(rd)
}
mtext("Year", 1, 1, outer=TRUE, cex=3)
mtext("Population Density", 2, 1, outer=TRUE, las=0, cex=3)
mtext("Basic2", 3, 1, outer=TRUE, cex=3.5)
parinfo  
dev.off()


xyplot(MDperKMsqFall_mean + HuntDen_All_mean ~ year |macrounit, data=Popdata, type="l",auto.key = list(space = "top", text = c("PopDen", "HuntDen"), points = FALSE, lines = TRUE))
length(which(Popdata$HuntDen_All_mean[cond]>Popdata$MDperKMsqFall_mean[cond]))




# On WholeAreaMeans
# opt_maf <- optim(par= c(rd=0.3, K=2), fn=maf, suppar=c(N0=N_0, steps=tsteps), data=WholeAreaMeans$MDperKMsqFall_mean)
# result_maf <- log_d(c(opt_maf$par[1], opt_maf$par[2]),c(N0=N_0,  steps=tsteps)) #MAF: 0.3056459, rd:0.05669221, K:2.12359422
# rd_maf <- c(((result_maf[-1]-result_maf[-tsteps])/result_maf[-tsteps]),NA)
# opt_mss <- optim(par= c(rd=0.3, K=2), fn=mss, suppar=c(N0=N_0, steps=tsteps), data=WholeAreaMeans$MDperKMsqFall_mean)
# result_mss <- log_d(c(opt_mss$par[1], opt_mss$par[2]),c(N0=N_0,  steps=tsteps)) #MSS:8.654571, rd:0.03246672, K:3.33934080
# rd_mss <- c(((result_mss[-1]-result_mss[-tsteps])/result_mss[-tsteps]),NA)
# 
# plot(rd_maf~result_maf, type="l", col="red")#correct:linear relationship of rd and N
# cor(rd_maf,result_maf, use="pairwise.complete.obs")#--1
# lines(rd_mss~result_mss, type="l", col="blue")#correct:linear
# cor(rd_mss,result_mss, use="pairwise.complete.obs") #-1
# 
# plot(WholeAreaMeans$MDperKMsqFall_mean, cex=0.5, main="Discrete logistic model")
# lines(result_maf, col="red") 
# lines(result_mss, col="blue")
# legend("topleft", legend=c("MAF", "MSS"), col=c("red", "blue"), lty=1)
