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