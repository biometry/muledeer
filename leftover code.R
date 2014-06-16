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