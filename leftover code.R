

### Check GAMs for Winter Mortality without outliers (09.07)
plot(AllMeans$WinterMort_mean~AllMeans$year) # 2 outliers
outlier <- which.max(AllMeans$WinterMort_mean)
outlier <- c(outlier, which.max(AllMeans$WinterMort_mean[-outlier]))#117,24
AllMeans2 <- AllMeans
AllMeans2$WinterMort_mean[outlier] <- NA
plot(AllMeans2$WinterMort_mean~AllMeans2$year)
gam_temp <- gam(WinterMort_mean ~ s(AvrgWinterMinTemp,  bs="cs"), data=AllMeans2)#slightly better
gam_hunt <- gam(WinterMort_mean ~ s(HuntDen_All_mean,  bs="cs"), data=AllMeans2)#not better
summary(gam_hunt)



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


###Calculating hunting density per area of badlands within each macrounit from badland area size per hunting unit
#load "total_harves_MDgunHarvest" into AllMeans
#unique(mules1962[which(mules1962$year == 1990),]$total_harvest_MDgunHarvest) verify that raw hunting data was produced on macrounit-scale
badlands <- data.frame(hunting_unit = c("4A","4B","4C","4D","4E","4F"), badlands_area_sqkm = c(411,545,568,495,243,127)*2.58998811, macrounit=c("0-4","0-3","0-3","0-2","0-2","0-1")) #from sq miles to sq km
badlands <- melt(t(tapply(badlands$badlands_area_sqkm, INDEX=list(badlands$macrounit), FUN=sum)))
AllMeans <- cbind(AllMeans, HuntDen_All = numeric(length(AllMeans$year)))
for (i in 1:4){
  cond <- which(AllMeans$macrounit == badlands[i,2])
  AllMeans$HuntDen_All[cond] = AllMeans$TotalHarvest_mean[cond] / badlands[i,3]
}
x <- AllMeans$HuntDen_All-AllMeans$HuntDen_All
# turned out to be completely unnecessary 
