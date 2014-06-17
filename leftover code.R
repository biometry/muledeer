###leftovers for plotting, not used because replaced by functions


###draw 4 plots per macrounit and gam-regression lines

par(mfrow=c(2,2),oma=c(2,0,2,0))
macrounits <- levels(AllMeans$macrounit)
for (i in 1:length(macrounits)){
  cond = which(AllMeans$macrounit==macrounits[i])  
  plot(AllMeans$MDperKMsqSpring_mean[cond]~AllMeans$year[cond], type="p", main=macrounits[i], xlab="Year", ylab="Density per km²")
  
  lines(x=AllMeans$year[cond], gam_allpred$fit[cond], col="orange")
  lines(x=AllMeans$year[cond], gam_allpred$fit[cond]  + 2 * gam_allpred$se.fit[cond], col="orange", lty=2)
  lines(x=AllMeans$year[cond], gam_allpred$fit[cond]  - 2 * gam_allpred$se.fit[cond], col="orange", lty=2)
}
title("GAM of whole Area", outer=TRUE)

### GAMs of spring MD Population for each macrounit 
# GAM without Autocorrelation
gamlist <- list()
gampredlist <- list()
for (i in 1:length(macrounits)){
  macrounits <- levels(AllMeans$macrounit)
  gamlist[[i]] <- gamm(MDperKMsqSpring_mean ~ s(year, fx = FALSE, k = -1,bs = "cr"), data=AllMeans, subset=(macrounit == macrounits[i]))
  gampredlist[[i]] <- predict(gamlist[[i]]$gam, se.fit=T, type="response")
}
#GAm with Autocorrelation
gam_corlist <- list()
gam_corpredlist <- list()
for (i in 1:length(macrounits)){
  macrounits <- levels(AllMeans$macrounit)
  gam_corlist[[i]] <- gamm(MDperKMsqSpring_mean ~ s(year, fx = FALSE, k = -1,bs = "cr"), data=AllMeans, correlation = corCAR1(form=~year), subset=(macrounit == macrounits[i]))
  gam_corpredlist[[i]] <- predict(gam_corlist[[i]]$gam, se.fit=T, type="response")
}

# Plot
par(mfrow=c(2,2),oma=c(2,0,2,0))
for (i in 1:length(levels(AllMeans$macrounit))){
  cond = which(AllMeans$macrounit==macrounits[i])  
  plot(AllMeans$MDperKMsqSpring_mean[cond]~AllMeans$year[cond], type="p", main=macrounits[i], xlab="Year", ylab="Density per km²")
  
  lines(x=AllMeans$year[cond], gampredlist[[i]]$fit, col="red") #no correlation
  lines(x=AllMeans$year[cond], gampredlist[[i]]$fit  + 2 * gampredlist[[i]]$se.fit, col="red", lty=2)
  lines(x=AllMeans$year[cond], gampredlist[[i]]$fit  - 2 * gampredlist[[i]]$se.fit, col="red", lty=2)
  
  lines(x=AllMeans$year[cond], gam_corpredlist[[i]]$fit, col="blue") #correlation 
  lines(x=AllMeans$year[cond], gam_corpredlist[[i]]$fit  + 2 * gam_corpredlist[[i]]$se.fit, col="blue", lty=2)
  lines(x=AllMeans$year[cond], gam_corpredlist[[i]]$fit  - 2 * gam_corpredlist[[i]]$se.fit, col="blue", lty=2)  
}
title("GAM per macrounit with/without autocorrelation", outer=TRUE)