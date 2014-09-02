
log_d <- function(par, suppar){
  output <- numeric(suppar[2])
  output[1] <- suppar[1]
  for (t in 1:(length(output)-1)){
    output[t+1] <- output[t]+output[t]*par[1]*(1-output[t]/par[2])
  }
  return(output)
}

#comparing basic 1 and statespace 1 in MU2

Popdata <- data.frame(MDperKMsqFall_mean=AllMeans$MDperKMsqFall_mean, year=AllMeans$year, macrounit=AllMeans$macrounit, RepRateFall_mean =AllMeans$RepRateFall_mean)#remove NAs
Popdata <- Popdata[which(complete.cases(Popdata)==TRUE),]

macrounits <- levels(Popdata$macrounit)
i=2
cond = which(Popdata$macrounit==macrounits[i]) 
N_0 <- Popdata$MDperKMsqFall_mean[cond[1]]
tsteps <- length(Popdata$MDperKMsqFall_mean[cond])
suppar=c(N_0, tsteps)
plot(Popdata$MDperKMsqFall_mean[cond]~Popdata$year[cond], xlab="",ylab="", cex.axis=2.5, cex.main = 3, main=macrounits[i], type="l", ylim=c(0,ylimmax[i]))

par1= c(rd=0.06, K=2.3)
basic1 <- log_d(par1,c(N0=N_0,  steps=tsteps))
lines(basic1~Popdata$year[cond], col="red",lwd=2)
par2 = c(0.19,20.1)
statespace1 <- log_d(par2,c(N0=N_0,  steps=tsteps))
lines(statespace1~Popdata$year[cond], col="red",lwd=2)

#STATESPACE: MEAN vs MEDIAN
Popdata <- data.frame(MDperKMsqFall_mean=AllMeans$MDperKMsqFall_mean, year=AllMeans$year, macrounit=AllMeans$macrounit, RepRateFall_mean =AllMeans$RepRateFall_mean)#remove NAs
Popdata <- Popdata[which(complete.cases(Popdata)==TRUE),]
macrounits <- levels(Popdata$macrounit)
ylimmax <- c(3.5,3.5,3.5,6.5)

png("StateSpace1_median_vs_max-mitmedian K green.png", width=2200, height=1500)

par(mfrow=c(2,2),oma=c(10,10,10,10),mar=c(10, 10, 10, 10))
for (i in 1:length(macrounits)){
  cond = which(Popdata$macrounit==macrounits[i]) 
  N_0 <- Popdata$MDperKMsqFall_mean[cond[1]]
  tsteps <- length(Popdata$MDperKMsqFall_mean[cond])
  suppar=c(N_0, tsteps)
  plot(Popdata$MDperKMsqFall_mean[cond]~Popdata$year[cond], xlab="",ylab="", cex.axis=2.5, cex.main = 3, main=macrounits[i], type="l", ylim=c(0,ylimmax[i]))
  
  par1= as.numeric(c(parinfo[1,i], parinfo[5,i]))
  mean <- log_d(par1,c(N0=N_0,  steps=tsteps))
  lines(mean~Popdata$year[cond], col="red",lwd=2)
  par2 = as.numeric(c(parinfo[4,i], parinfo[5,i]))
  median <- log_d(par2,c(N0=N_0,  steps=tsteps))
  lines(median~Popdata$year[cond], col="blue",lwd=2)
  par3= as.numeric(c(parinfo[1,i], parinfo[8,i]))
  bothmedian <- log_d(par3,c(N0=N_0,  steps=tsteps))
  lines(bothmedian~Popdata$year[cond], col="green",lwd=2)
  }
mtext("Year", 1, 1, outer=TRUE, cex=3)
mtext("Population Density", 2, 1, outer=TRUE, las=0, cex=3)
mtext("StateSpace1 - mean vs. median", 3, 1, outer=TRUE, cex=3.5)
dev.off()

###STATESPACE2----
Popdata <- data.frame(MDperKMsqFall_mean=AllMeans$MDperKMsqFall_mean, year=AllMeans$year, macrounit=AllMeans$macrounit, RepRateFall_mean =AllMeans$RepRateFall_mean)#remove NAs
Popdata <- Popdata[which(complete.cases(Popdata)==TRUE),]
macrounits <- levels(Popdata$macrounit)
ylimmax <- c(3.5,3.5,3.5,6.5)
index_mean <-  c(1,3,5,7)
index_median <- c(2,4,6,8)
png("StateSpace1_median_vs_max-mitmedian K green.png", width=2200, height=1500)
par(mfrow=c(2,2),oma=c(10,10,10,10),mar=c(10, 10, 10, 10))
for (i in 1:length(macrounits)){
  cond = which(Popdata$macrounit==macrounits[i]) 
  N_0 <- Popdata$MDperKMsqFall_mean[cond[1]]
  tsteps <- length(Popdata$MDperKMsqFall_mean[cond])
  suppar=c(N_0, tsteps)
  plot(Popdata$MDperKMsqFall_mean[cond]~Popdata$year[cond], xlab="",ylab="", cex.axis=2.5, cex.main = 3, main=macrounits[i], type="l", ylim=c(0,ylimmax[i]))
#   lines(N_list[[index_mean[i]]]~Popdata$year[cond], col="red", lwd=2)  
#   lines(N_list[[index_median[i]]]~Popdata$year[cond],lty=2,col="red",lwd=2)
  par1= as.numeric(c(parinfo[1,i], parinfo[5,i]))
  mean <- log_d(par1,c(N0=N_0,  steps=tsteps))
  lines(mean~Popdata$year[cond], col="red",lwd=2)
  par2 = as.numeric(c(parinfo[4,i], parinfo[5,i]))
  median <- log_d(par2,c(N0=N_0,  steps=tsteps))
  lines(median~Popdata$year[cond], col="blue",lwd=2)
  par3= as.numeric(c(parinfo[4,i], parinfo[8,i]))
  bothmedian <- log_d(par3,c(N0=N_0,  steps=tsteps))
  lines(bothmedian~Popdata$year[cond], col="green",lwd=2)
}
mtext("Year", 1, 1, outer=TRUE, cex=3)
mtext("Population Density", 2, 1, outer=TRUE, las=0, cex=3)
mtext("StateSpace2 - mean vs. median", 3, 1, outer=TRUE, cex=3.5)
dev.off()


