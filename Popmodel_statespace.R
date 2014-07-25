AllMeans <- read.csv("AllMeans.csv")
Popdata <- data.frame(MDperKMsqFall_mean=AllMeans$MDperKMsqFall_mean, year=AllMeans$year, macrounit=AllMeans$macrounit, HuntDen_All_mean=AllMeans$HuntDen_All_mean)#remove NAs
Popdata <- Popdata[which(complete.cases(Popdata)==TRUE),]
library(rjags)


# 5.2 A simple model ----------------------------------------

# ich habe schon versucht, die Startwerte und priors zu erweitern/verkleinern
# in der rjags-manual (auch auf github) auf seite 6 steht etwas zu den "nodes"
# nämlich, dass jeder prior eine "node" wäre. Die Fehlermeldungen kommen jedoch 
# auch z.B. bei node[13], und es gibt ja keine 13 priors



# Specify the model in JAGS:
sink("model_ss.jags")
cat("
    model{
    # Priors
    N.est[1] ~ dunif(0.1, 10)             # Prior for initial population density (max in dataset:6.6)
    mean.rmax ~ dunif(0.1, 10)            # Prior for mean rmax
    K ~ dunif(0.1, 100)                     # prior for carrying capacity
    sigma.proc ~ dunif(0.1, 10)           # Prior for sd of state process
    sigma2.proc <- pow(sigma.proc, 2)   # Variance = sd^2
    tau.proc <- pow(sigma.proc, -2)     # Tau = sqrt(sd)
    sigma.obs ~ dunif(0.1, 10)           # Prior for sd of observation process
    sigma2.obs <- pow(sigma.obs, 2)
    tau.obs <- pow(sigma.obs, -2)
    
    # Likelihood
    # State process
    for (t in 1:(T-1)){
    rmax[t] ~ dnorm(mean.rmax, tau.proc)        # Normal distribution in JAGS is specified with mean and tau (*not* sd!)
    N.est[t+1] <- N.est[t]+(N.est[t]*rmax[t]*(1-(N.est[t]/K)))
    }
    # Observation process
    for (t in 1:T){
    y[t] ~ dnorm(N.est[t], tau.obs)
    }
    }
    ", fill=TRUE)
sink()

i=1
macrounits <- levels(Popdata$macrounit)
cond = which(Popdata$macrounit==macrounits[i]) 
tsteps <- length(Popdata$MDperKMsqFall_mean[cond])
# Bundle the data:
model.data <- list(y = Popdata$MDperKMsqFall_mean[cond], T = tsteps)

# Initial values:
inits <- function(){list(sigma.proc = runif(1, 0.1, 10), mean.rmax = runif(1, 0.1, 10), K = runif(1,0.1,10), sigma.obs = runif(1, 0.1, 10), N.est = c(runif(1, 0.1, 10), rep(NA, (tsteps-1))))}

# Compile the model and run MCMC for burn-in phase:
model_fit <- jags.model(file = "model_ss.jags", data = model.data, inits = inits, n.chains = 3, n.adapt = 5000)

# Specify parameters whose posterior values are to be saved:
parameters <- c("rmax","mean.rmax", "sigma2.obs", "sigma2.proc", "K", "N.est", "y")

# Continue running the MCMC to produce posterior distributions:
result_sss <- coda.samples(model = model_fit, variable.names = parameters, n.iter = 5000)

plot(result_sss) # Produces many plots: 25 for N.est, 24 for rmax, and one for mean.lambda, sigma2.proc and sigma2.obs, respectively
summary(result_sss)

# Make a nice graph with data and model result_sss:
res_table <- as.data.frame(summary(result_sss)$quantiles)
names(res_table) <- c("q2.5", "q25", "q50", "q75", "q97.5")
N.est_median <- res_table$q50[which(rownames(res_table) == "N.est[1]"):which(rownames(res_table) == "N.est[25]")]
N.est_lower <- res_table$q2.5[which(rownames(res_table) == "N.est[1]"):which(rownames(res_table) == "N.est[25]")]
N.est_upper <- res_table$q97.5[which(rownames(res_table) == "N.est[1]"):which(rownames(res_table) == "N.est[25]")]
m1 <- min(c(y, N.est_median, N, N.est_lower))
m2 <- max(c(y, N.est_median, N, N.est_upper))
par(mar = c(4.5, 4, 1, 1), cex = 1.2)
plot(0, 0, ylim = c(m1, m2), xlim = c(0.5, n.years), ylab = "Population size", xlab = "Year", las = 1, col = "black", type = "l", lwd = 2, frame = FALSE, axes = FALSE)
axis(2, las = 1)
axis(1, at = seq(0, n.years, 5), labels = seq(0, n.years, 5))
axis(1, at = 0:n.years, labels = rep("", n.years + 1), tcl = -0.25)
polygon(x = c(1:n.years, n.years:1), y = c(N.est_lower, N.est_upper[n.years:1]), col = "gray90", border = "gray90")
lines(1:n.years, N, col = "red", lwd = 2)
lines(1:n.years, y, col = "black", lwd = 2)
lines(1:n.years, N.est_median, col = "blue", lwd = 2)
legend(x = 1, y = m2, legend = c("True", "Observed", "Estimated"), lty = c(1, 1, 1), lwd = c(2, 2, 2), col = c("red", "black", "blue"), bty = "n", cex = 1)

savePlot("ibex_model_data.png")






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

