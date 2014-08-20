
AllMeans <- read.csv("AllMeans.csv")
Popdata <- data.frame(MDperKMsqFall_mean=AllMeans$MDperKMsqFall_mean, year=AllMeans$year, macrounit=AllMeans$macrounit, HuntDen_All_mean=AllMeans$HuntDen_All_mean)#remove NAs
Popdata <- Popdata[which(complete.cases(Popdata)==TRUE),]
library(rjags)




# 5.2 A simple model ----------------------------------------


# Specify the model in JAGS:
sink("model_ss.jags")
cat("
    model{
    # Priors
    N.est[1] ~ dunif(0.00001, 10)             # Prior for initial population density (max in dataset:6.6)
    mean.rmax ~ dunif(0.00001, 4)            # Prior for mean rmax (geommean in dataset 0.2, max 3.4)
    K ~ dunif(0.00001, 100)                     # prior for carrying capacity
    sigma.proc ~ dunif(0.00001, 50)           # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)         # Variance = sd^2
    tau.proc <- pow(sigma.proc, -2)           # Tau = sqrt(sd)
    sigma.obs ~ dunif(0.00001, 50)            # Prior for sd of observation process
    sigma2.obs <- pow(sigma.obs, 2)
    tau.obs <- pow(sigma.obs, -2)
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){
    rmax[t] ~ dnorm(mean.rmax, tau.proc)T(0, 10)        # Truncated normal distribution to prevent values < 0
    rd[t] <- rmax[t]*(1-N.est[t]/K)
    N.est[t+1] <- N.est[t]+N.est[t]*rd[t]
    }
    # Observation process
    for (t in 1:T){
    y[t] ~ dnorm(N.est[t], tau.obs)
    }
    }
    ", fill=TRUE)
sink()



png("StateSpace1_15000_5000_final.png", width=2200, height=1500)
par(mfrow=c(2,2),oma=c(5,5,5,5),mar=c(10, 10, 10, 10))
parinfo <- data.frame("0-1" = numeric(11), "0-2" = numeric(11), "0-3" = numeric(11), "0-4" = numeric(11), row.names=c("mean of rmax","mean of sd.rmax","mean of se.rmax","mean of K", "sd.K", "se.K", "mean of sigma2.obs", "mean of sigma2.proc", "mean of sd.N.est","mean of se.N.est","MSS"))
macrounits <- levels(Popdata$macrounit)
ylimmax <- c(3.5,3.5,3.5,6.5)
for (i in 1:length(macrounits)){
  
  cond = which(Popdata$macrounit==macrounits[i]) 
  tsteps <- length(Popdata$MDperKMsqFall_mean[cond])
  # Bundle the data:
  model.data <- list(y = Popdata$MDperKMsqFall_mean[cond], T = tsteps)
  
  # Initial values:
  #no idea why these dont work as the follow the priors' conditions
  #inits <- function(){list(sigma.proc = runif(1, 0.00001, 50), mean.rmax = runif(1, 0.00001, 4), K = runif(1, 0.00001, 100), sigma.obs = runif(1, 0.00001, 50), N.est = c(runif(1, 0.00001, 10), rep(NA, (tsteps-1))))}
  inits <- function(){list(sigma.proc = 1, mean.rmax = 0.3, K = 2, sigma.obs = 1, N.est = c(1, rep(NA, (tsteps-1))))}
  
  
  
  # Compile the model and run MCMC for burn-in phase:
  model_fit <- jags.model(file = "model_ss.jags", data = model.data, inits = inits, n.chains = 1, n.adapt = 15000)
  
  # Specify parameters whose posterior values are to be saved:
  parameters <- c("rmax","mean.rmax", "rd", "sigma2.obs", "sigma2.proc", "K", "N.est", "y")
  
  # Continue running the MCMC to produce posterior distributions:
  result_sss <- coda.samples(model = model_fit, variable.names = parameters, n.iter = 5000)
  
  # plot(result_sss) # Produces many plots: 25 for N.est, 24 for rmax, and one for mean.lambda, sigma2.proc and sigma2.obs, respectively
  # summary(result_sss)
  
  
  #extract data
  
  summ <- summary(result_sss)
  res_means <- as.data.frame(summ$statistics)
  res_quant <- as.data.frame(summ$quantiles)
  
  N.est_mean <- res_means[which(rownames(res_means) == "N.est[1]"):(which(rownames(res_means) == "N.est[1]")+tsteps-1),] # plus tstep-1 in order to keep flexible for different lengths of data sets due to NAs
  y_mean <- res_means[which(rownames(res_means) == "y[1]"):(which(rownames(res_means) == "y[1]")+tsteps-1),]
  rd_mean <- res_means[which(rownames(res_means) == "rd[1]"):(which(rownames(res_means) == "rd[1]")+tsteps-2),]
  rd_calc <- c((N.est_mean[-1,1]-N.est_mean[-tsteps,1])/N.est_mean[-tsteps,1],NA) #(N(t+1)-N(t))/N(t)
  rmax_mean <- res_means[which(rownames(res_means) == "rmax[1]"):(which(rownames(res_means) == "rmax[1]")+tsteps-2),]
  K_mean <- res_means[which(rownames(res_means) == "K"),]
  sigmas_mean <- res_means[c(which(rownames(res_means) == "sigma2.obs"), which(rownames(res_means) == "sigma2.proc")),]
  
  
  N.est_quant <- res_quant[which(rownames(res_quant) == "N.est[1]"):(which(rownames(res_means) == "N.est[1]")+tsteps-1),]
  
  MSE <-   sum((abs(N.est_mean[,1]-Popdata$MDperKMsqFall_mean[cond]))^2, na.rm=TRUE)/length(Popdata$MDperKMsqFall_mean[cond])
  
  # Make a nice graph with data and model result:
  
  plot(Popdata$MDperKMsqFall_mean[cond]~Popdata$year[cond], cex=0.5,  cex.axis=2.5, cex.main = 3, ylab = "", xlab = "", las = 1, col = "black", type = "p", main=macrounits[i], ylim=c(0,ylimmax[i]))
  polygon(x = c(Popdata$year[cond], rev(Popdata$year[cond])), y = c(N.est_quant[,5], rev(N.est_quant[,1])), col = "gray90", border = "gray90")
  
  lines (y_mean[,1]~Popdata$year[cond], col="black", lwd=1)
  lines (N.est_mean[,1]~Popdata$year[cond], col="red", lwd=2)
  #legend("topleft", legend = c("Observed", "Estimated", "95% Quantile Estimated"), lty = c(1, 1, 1), lwd = c(1,1,1), col = c("black", "red", "grey"), bty = "n", cex = 1)
  
  # display sd of estimate
  #polygon(x = c(Popdata$year[cond], rev(Popdata$year[cond])), y = c((N.est_mean[,1]+2*N.est_mean[,2]), rev(N.est_mean[,1]-2*N.est_mean[,2])), col = "gray90", border = "gray90")
  #lines (N.est_mean[,1]+2*N.est_mean[,2]~Popdata$year[cond], col="red", lty=2)
  #lines (N.est_mean[,1]-2*N.est_mean[,2]~Popdata$year[cond], col="red", lty=2)
  #legend("topleft", legend = c("Observed", "Estimated", "sd Estimated", "95% CRI"), lty = c(1, 1, 2, 1), lwd = c(2, 2, 1, 2), col = c("black", "red", "red", "grey"), bty = "n", cex = 1)
  
  
  # rd_mean_sd1 <- c(rd_mean[,1]+2*rd_mean[,2],NA)
  # rd_mean_sd2 <- c(rd_mean[,1]-2*rd_mean[,2],NA)
  # m1 <- min(c(rd_mean_sd1,rd_mean_sd2), na.rm=TRUE)
  # m2 <- max(c(rd_mean_sd1,rd_mean_sd2), na.rm=TRUE)
  # plot(c(rd_mean[,1],NA)~N.est_mean[,1], type="l", lty=1, ylim=c(m1,m2), xlab="Mean of Predicted Population Density", ylab="Mean of Predicted rd", main=macrounits[i]) #logical that shape is not linear because mean of rd is based on different estimates of rmax in each knot
  # lines (rd_mean_sd1~N.est_mean[,1], col="black", lty=2)
  # lines (rd_mean_sd2~N.est_mean[,1], col="black", lty=2)
  # lines (rd_calc~N.est_mean[,1], col="red")
  parinfo[,i] <- c(mean(rmax_mean[,1]), mean(rmax_mean[,2]),mean(rmax_mean[,4]),K_mean[,1],K_mean[,2],K_mean[,4],sigmas_mean[,1], mean(N.est_mean[,2]),mean(N.est_mean[,4]),MSE)
  }
mtext("Year", 1, 1, outer=TRUE, cex=3)
mtext("Population Density", 2, 1, outer=TRUE, las=0, cex=3)
mtext("StateSpace1", 3, 1, outer=TRUE, cex=3.5)       
parinfo  
dev.off()
# #exclude 2 outliers of rmax
# max1 <- which.max(Popdata$RepRateFall_mean)
# max2 <- which.max(Popdata$RepRateFall_mean[-max1])
# hist(Popdata$RepRateFall_mean[-c(max1,max2)])#normally distributed?
# rmax <- max(Popdata$RepRateFall_mean[-c(max1,max2)])
# rmax.sd <- sd(Popdata$RepRateFall_mean[-c(max1,max2)])
# #problem: max(r) != rmax 
# 
# # observation error
# obs.sd <- sd(Popdata$MDperKMsqFall_mean)#0.9217026
# hist(Popdata$MDperKMsqFall_mean, breaks=100)# somwhere near normal distribution
# 
# 
# par(mfrow=c(2,4),oma=c(2,0,2,0))
# macrounits <- levels(Popdata$macrounit)
# result_ss_out <- list()
# rd_ss_out <- list()
# parinfo_ss <- data.frame("0-1" = numeric(3), "0-2" = numeric(3), "0-3" = numeric(3), "0-4" = numeric(3), row.names=c("rd_ss", "K", "MSS"))
# 
# 
# for (i in 1:length(macrounits)){
#   cond = which(Popdata$macrounit==macrounits[i])  
#   
#   N_0 <- Popdata$MDperKMsqFall_mean[cond[1]]
#   tsteps <- length(Popdata$MDperKMsqFall_mean[cond])
#   opt_ss <- optim(par= c(rd_ss=0.3, K=2), fn=mss_ss, method="L-BFGS-B", lower=c(0 , 0), suppar=c(N_0, tsteps, 0.1), data=Popdata$MDperKMsqFall_mean[cond])
#   result_ss <- log_d_ss(c(opt_ss$par[1], opt_ss$par[2]),c(N0=N_0,  steps=tsteps, 0.1))
#   names(result_ss) <- Popdata$year[cond]
#   result_ss_out <- append(result_ss_out, list(c(result_ss))) 
#   names(result_ss) <- Popdata$year[cond]
#   rd_ss <- c(((result_ss[-1]-result_ss[-tsteps])/result_ss[-tsteps]),NA)
#   rd_ss_out <- append(rd_ss_out, list(c(rd_ss)))
#   parinfo_ss[,i] <- c(opt$par, opt$value)
#   
#   
#   plot(Popdata$MDperKMsqFall_mean[cond]~Popdata$year[cond], cex=0.5, main=macrounits[i], xlab="Year",ylab="Population Density")
#   lines(result_ss~Popdata$year[cond], col="red")
#   #lines(result_out[cond]~Popdata$year[cond], col="blue")
#   plot(rd_ss~result_ss, type="l", lty=2, main=macrounits[i], xlab="Predicted Population Density",ylab="Predicted Reproduction rate")
#   #lines(rd_out[cond]~result_out[cond], col="blue",lty=2)
#   remove(result_ss)
#   remove(rd_ss)
# }
# title("Basic Discrete Logistic Model including Observation Error - Population Densities and Reproduction Rates (MSS_SS)", outer=TRUE)       
# parinfo_ss  
