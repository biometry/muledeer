
#AllMeans <- read.csv("AllMeans.csv")
Popdata <- data.frame(MDperKMsqFall_mean=AllMeans$MDperKMsqFall_mean, year=AllMeans$year, macrounit=AllMeans$macrounit, HuntDen_All_mean=AllMeans$HuntDen_All_mean)#remove NAs
Popdata <- Popdata[which(complete.cases(Popdata)==TRUE),]
library(rjags)




# 5.2 A simple model ----------------------------------------


# Specify the model in JAGS:
sink("model_ss.jags")
cat("
    model{
    #Priors    
    N.est[1] ~ dunif(0.00001, 10)             # Prior for initial population density (max in dataset:6.6)
    mean.rmax ~ dunif(0.00001, 4)            # Prior for mean rmax (geommean in dataset 0.2, max 3.4)
    K ~ dunif(0.00001, 100)                     # prior for carrying capacity
    sigma.proc ~ dunif(0.00001, 50)           # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)         # Variance = sd^2
    tau.proc <- pow(sigma.proc, -2)           # Tau = sqrt(sd)
    sigma.obs ~ dunif(0.00001, 50)            # Prior for sd of observation process
    sigma2.obs <- pow(sigma.obs, 2)
    tau.obs <- pow(sigma.obs, -2)
    pH ~ dunif(0,1)    
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){
    rmax[t] ~ dnorm(mean.rmax, tau.proc)T(0, 10)        # Truncated normal distribution to prevent values < 0
    rd[t] <- rmax[t]*(1-N.est[t]/K)
    N.est[t+1] <- N.est[t]+N.est[t]*rd[t]-pH*HuntDen[t]
    }
    # Observation process
    for (t in 1:T){
    y[t] ~ dnorm(N.est[t], tau.obs) 
    }
    }
    ", fill=TRUE)
sink()



png("StateSpace2_45000_90000.png", width=2200, height=1500)
par(mfrow=c(2,2),oma=c(5,5,5,5),mar=c(10, 10, 10, 10))
parinfo <- data.frame("0-1" = numeric(19), "0-2" = numeric(19), "0-3" = numeric(19), "0-4" = numeric(19), row.names=c("mean of mean.rmax","sd. mean.rmax","se.mean.rmax", "mean of K", "sd.K", "se.K", "mean of sigma2.obs", "mean of sigma2.proc", "mean of sd.N.est","mean of se.N.est","mean of pH","sd.pH", "se.pH","2.5% quantile pH", "25% quantile pH","50% quantile pH","75% quantile pH","97.5% quantile pH","MSS"))
macrounits <- levels(Popdata$macrounit)
ylimmax <- c(3.5,3.5,3.5,6.5)
result_final <- list()
GelmanDiag <- list()


for (i in 1:length(macrounits)){
  
  cond = which(Popdata$macrounit==macrounits[i]) 
  tsteps <- length(Popdata$MDperKMsqFall_mean[cond])
  # Bundle the data:
  model.data <- list(y = Popdata$MDperKMsqFall_mean[cond], T = tsteps, HuntDen = Popdata$HuntDen_All_mean)
  
  # Initial values:
  # inits <- function(){list(sigma.proc = runif(1, 0, 10), mean.rmax = runif(1, 0.01, 10), K = runif(1, 0.1, 10), sigma.obs = runif(1, 0, 10), N.est = c(runif(1, 0.1, 10), rep(NA, (tsteps-1))))}
  
  inits1 <- function(){list(sigma.proc = 1, mean.rmax = 0.3, K = 2, sigma.obs = 1, N.est = c(1, rep(NA, (tsteps-1))), pH=0.1);
  list(sigma.proc = 1, mean.rmax = 5, K = 50, sigma.obs = 1, N.est = c(1, rep(NA, (tsteps-1))), pH=0.1)
  list(sigma.proc = 1, mean.rmax = 3, K = 25, sigma.obs = 1, N.est = c(1, rep(NA, (tsteps-1))), pH=0.1)}
  
  
  # Compile the model and run MCMC for burn-in phase:
  model_fit <- jags.model(file = "model_ss.jags", data = model.data, inits = inits, n.chains = 3, n.adapt = 45000)
  
  # Specify parameters whose posterior values are to be saved:
  parameters <- c("mean.rmax", "sigma2.obs", "sigma2.proc", "K", "N.est","pH")
  
  # Continue running the MCMC to produce posterior distributions:
  result_sss <- coda.samples(model = model_fit, variable.names = parameters, n.iter = 90000)
  
  # Extract relevant data for Convergence Diagnostics 
  b <- list(length(result_sss))
  for (g in 1:3){#all 3 chains
    incl <- which (!grepl("N.est", varnames(result_sss[[g]])))#exclude N.est for Gelman Diag
    b[[g]] <- result_sss[[g]][,incl] 
    varnames(b[[g]]) <- varnames(result_sss[[g]][,incl])
  }
  result_final <- append(result_final,b)
  GelmanDiag <- append(GelmanDiag, gelman.diag(result_sss))
  
  
  #extract data
  
  summ <- summary(result_sss)
  res_means <- as.data.frame(summ$statistics)
  res_quant <- as.data.frame(summ$quantiles)
  
  N.est_mean <- res_means[which(rownames(res_means) == "N.est[1]"):(which(rownames(res_means) == "N.est[1]")+tsteps-1),] # plus tstep-1 in order to keep flexible for different lengths of data sets due to NAs
  mean.rmax_mean <- res_means[which(rownames(res_means) == "mean.rmax"),]
  K_mean <- res_means[which(rownames(res_means) == "K"),]
  sigmas_mean <- res_means[c(which(rownames(res_means) == "sigma2.obs"), which(rownames(res_means) == "sigma2.proc")),]
  pH_mean <-   res_means[which(rownames(res_means) == "pH"),]
  
  N.est_quant <- res_quant[which(rownames(res_quant) == "N.est[1]"):(which(rownames(res_means) == "N.est[1]")+tsteps-1),]
  pH_quant <- res_quant[which(rownames(res_quant) == "pH"),]
  
  MSE <-   sum((abs(N.est_mean[,1]-Popdata$MDperKMsqFall_mean[cond]))^2, na.rm=TRUE)/length(Popdata$MDperKMsqFall_mean[cond])
  
  # Make a nice graph with data and model result:
  
  plot(Popdata$MDperKMsqFall_mean[cond]~Popdata$year[cond], cex=0.5,  cex.axis=2.5, cex.main = 3, ylab = "", xlab = "", las = 1, col = "black", type = "p", main=macrounits[i], ylim=c(0,ylimmax[i]))
  polygon(x = c(Popdata$year[cond], rev(Popdata$year[cond])), y = c(N.est_quant[,5], rev(N.est_quant[,1])), col = "gray90", border = "gray90")
  
  lines (Popdata$MDperKMsqFall_mean[cond]~Popdata$year[cond], col="black", lwd=1)
  lines (N.est_mean[,1]~Popdata$year[cond], col="red", lwd=2)
  
  #legend("topleft", legend = c("Observed", "Estimated", "95% Quantile Estimated"), lty = c(1, 1, 1), lwd = c(1,1,1), col = c("black", "red", "grey"), bty = "n", cex = 1)
  # display sd of estimate
  #polygon(x = c(Popdata$year[cond], rev(Popdata$year[cond])), y = c((N.est_mean[,1]+2*N.est_mean[,2]), rev(N.est_mean[,1]-2*N.est_mean[,2])), col = "gray90", border = "gray90")
  #lines (N.est_mean[,1]+2*N.est_mean[,2]~Popdata$year[cond], col="red", lty=2)
  #lines (N.est_mean[,1]-2*N.est_mean[,2]~Popdata$year[cond], col="red", lty=2)
  #legend("topleft", legend = c("Observed", "Estimated", "sd Estimated", "95% CRI"), lty = c(1, 1, 2, 1), lwd = c(2, 2, 1, 2), col = c("black", "red", "red", "grey"), bty = "n", cex = 1)
  
  
  
  parinfo[,i] <- format(c(mean.rmax_mean[,1], mean.rmax_mean[,2],mean.rmax_mean[,4],K_mean[,1],K_mean[,2],K_mean[,4], sigmas_mean[,1], mean(N.est_mean[,2]),mean(N.est_mean[,4]),pH_mean[,1], pH_mean[,2],pH_mean[,4],pH_quant[1:5],MSE),scientific=FALSE)
}
mtext("Year", 1, 1, outer=TRUE, cex=3)
mtext("Population Density", 2, 1, outer=TRUE, las=0, cex=3)
mtext("StateSpace2", 3, 1, outer=TRUE, cex=3.5)       
parinfo  <- format(parinfo, scientific=FALSE)
parinfo
dev.off()

#save Gelman plots
m=0
for (k in c(1,4,7,10)){
  m<- m+1
  name <- paste(c("StateSpace2_45000_90000_Diag_",m,".png"),collapse="")
  png(name, width=2200, height=1500)
  par(mfrow=c(2,3),oma=c(5,5,5,5),mar=c(10, 10, 10, 10)) 
  gelman.plot(result_final[k:(k+2)], auto.layout=FALSE, cex.axis=2.5,cex.main=3,cex=5)
  mtext(paste(c("Gelman Diagnostics: StateSpace2 - MU", as.character(m)), collapse = " "), 3, 1, outer=TRUE, cex=3) 
  dev.off()
}

#save histogramm of pH (always 1st chain)
m=0
png("StateSpace2_45000_90000_pH_histo.png", width=2200, height=1500)
par(mfrow=c(2,2),oma=c(5,5,5,5))
for (k in c(1,4,7,10)){
  m<- m+1
  xmin <- as.numeric(parinfo[14,m])-0.1
  xmax <- as.numeric(parinfo[18,m])+0.1
  hist(result_final[[k]][,"pH"],breaks=200, main=m,xlab="",cex.axis=2.5,cex.main=3,ylab="",xlim=c(xmin,xmax))
  abline(v=as.numeric(parinfo[14,m]), col="red", lwd=3)
  abline(v=as.numeric(parinfo[18,m]), col="red", lwd=3)
}
mtext("pH in StateSpace2", 3, 1, outer=TRUE, cex=3) 
dev.off()



#save density plot of pH (always 1st chain)
m=0
png("StateSpace2_45000_90000_pH_1stchain.png", width=4000, height=4000)
par(mfrow=c(2,2),oma=c(10,10,10,10),mar=c(10,10,10,10))
for (k in c(1,4,7,10)){
  m<- m+1
  d <- density(result_final[[k]][,"pH"])
  xmin <- as.numeric(parinfo[14,m])-0.1
  xmax <- as.numeric(parinfo[18,m])+0.1
  plot(d, main=m,,cex.axis=6,cex.main=5,lwd=3,ylab="",xlab="",xlim=c(xmin,xmax))
  abline(v=as.numeric(parinfo[14,m]), col="red", lwd=3)
  abline(v=as.numeric(parinfo[18,m]), col="red", lwd=3)
}
mtext("pH in StateSpace2", 3, 1, outer=TRUE, cex=6) 
mtext("Value", 1, 1, outer=TRUE, cex=6)
mtext("Density", 2, 1, outer=TRUE, las=0, cex=6)
dev.off()



#save Gelman Diagnostics + parameter values
name2 <-paste(c("StateSpace2_30000_60000_GelmanDiag",".txt"),collapse="")
sink(name2, type=c("output","message"))
print(parinfo)
print(GelmanDiag)
sink()


print("DO YOU WANT MCMCPLOT HTML-FILES?")
print("DO YOU?")
print("DO YOU?")
print("DO YOU?")
save posterior distr+traces in HTML (messy with all 4 MUs, keps overwriting files)
library(mcmcplots)
# m=0
# for (k in c(1,4,7,10)){
#   m<- m+1
#   name3 <- paste(c("StateSpace2_15000_30000_Posterior_",m),collapse="")
#   dire <-   paste(c("C:/Users/lschmidt/Bachelorarbeit/Mule_deer/Popmodels/StateSpace2/30000_60000_3 chains dif initial/HTML",m),collapse="")
#   mcmcplot(result_final[k:(k+2)],filename=name3,style="plain",title=name3,dir=dire)
# }




