### GAMs applied on netto-reproduction-rate, including explanatory variables one by one

library(mgcv)

sink("GAM_results_RepRate.txt", type=c("output","message"))
###effect of population density (=Density-Dependence)
gam_dens <- gam(RepRateFall_mean ~ s(MDperKMsqFall_mean, bs="cs") , data=AllMeans) 
gam_dens <- gam(RepRateFall_mean ~ s(MDperKMsqFall_mean, by=macrounit, bs="cs") , data=AllMeans) 

gam_denspred <- data.frame(MDperKMsqFall_mean=AllMeans$MDperKMsqFall_mean, macrounit=AllMeans$macrounit)
gam_denspred <- order(gam_denspred)
gam_denspred <- cbind(gam_denspred, predict(gam_dens, se.fit=T, newdata=data.frame("MDperKMsqFall_mean"=AllMeans$MDperKMsqFall_mean, "macrounit"=AllMeans$macrounit), type="response"))
macrounitplots(glmobject = gam_denspred,title="gam_dens fall - Relationship Population Density -- Reproduction Rate",xcol="MDperKMsqFall_mean", ycol= "RepRateFall_mean", colour="red", shape="l")

summary(gam_dens)
AIC(gam_dens)
plot(gam_dens)
par(oma=c(2,0,2,0))
gam.check(gam_dens)


###Effect of Average Minimum Winter Temperature

gam_temp <- gam(RepRateFall_mean ~ s(AvrgWinterMinTemp, bs="cs"), data=AllMeans)#1
gam_temp <- gam(RepRateFall_mean ~ s(MDperKMsqFall_mean, bs="cs") + s(AvrgWinterMinTemp, by=macrounit, bs="cs") + macrounit, data=AllMeans)#2
gam_temp <- gam(RepRateFall_mean ~ s(MDperKMsqFall_mean, by=macrounit, bs="cs") +s(AvrgWinterMinTemp, bs="cs") + macrounit, data=AllMeans)#3


g#am_temppred <- data.frame(MDperKMsqFall_mean=AllMeans$MDperKMsqFall_mean)
#gam_temppred <- cbind(gam_temppred, predict(gam_temp, se.fit=T, newdata=data.frame("MDperKMsqFall_mean"=AllMeans$MDperKMsqFall_mean, "macrounit"=AllMeans$macrounit), type="response"))
macrounitplots(glmobject = gam_temppred,xcol="AvrgWinterMinTemp",title="gam_temp fall - effect of Mean Winter Minimum Temperature",colour="red")
macrounitplots(glmobject = gam_temppred,title="gam_temp fall - effect of Mean Winter Minimum Temperature",colour="red")
summary(gam_temp)
AIC(gam_temp)
par(oma=c(2,0,2,0))
gam.check(gam_temp)


# Check Residuals for temporal autocorrelation
gam_tempres <- residuals(gam_temp, type = "deviance")
plot(gam_tempres ~AllMeans$MDperKMsqFall_mean[which(!is.na(AllMeans$RepRateFall_mean))]) #
acf(gam_tempres, na.action = na.pass,main = "Auto-correlation plot for residuals gam_temp fall")



###Effect of hunting Density on each macrounit
gam_hunt <- gam(RepRateFall_mean ~ s(HuntDen_All_mean, bs="cs"), data=AllMeans)#1
gam_hunt <- gam(RepRateFall_mean ~ s(MDperKMsqFall_mean, bs="cs") + s(HuntDen_All_mean, by=macrounit, bs="cs") + macrounit, data=AllMeans)#2
gam_hunt <- gam(RepRateFall_mean ~ s(MDperKMsqFall_mean, by=macrounit, bs="cs") + s(HuntDen_All_mean, bs="cs") + macrounit, data=AllMeans)#3

gam_huntpred <- data.frame(MDperKMsqFall_mean=AllMeans$MDperKMsqFall_mean, macrounit=AllMeans$macrounit, HuntDen_All_mean=AllMeans$HuntDen_All_mean)
gam_huntpred <- cbind(gam_huntpred, predict(gam_hunt, se.fit=T, newdata=data.frame("MDperKMsqFall_mean"=AllMeans$MDperKMsqFall_mean, "macrounit"=AllMeans$macrounit, "HuntDen_All_mean"=AllMeans$HuntDen_All_mean), type="response"))
macrounitplots(glmobject = gam_huntpred,xcol="HuntDen_All_mean",title="gam_hunt fall - effect of Hunting Density",colour="red")
macrounitplots(glmobject = gam_huntpred,title="gam_hunt fall - effect of Hunting Density",colour="red")
summary(gam_hunt)
AIC(gam_hunt)

gam_huntres <- residuals(gam_hunt, type = "deviance")
plot(gam_huntres ~AllMeans$MDperKMsqFall_mean[which(!is.na(AllMeans$RepRateFall_mean))]) #
acf(gam_huntres, na.action = na.pass,main = "Auto-correlation plot for residuals gam_hunt fall")

gam.check(gam_hunt)


###Effect of Oil Well Density on each macrounit


gam_oil <- gam(RepRateFall_mean ~ s(WellDen_mean, bs="cs"), data=AllMeans)#1
gam_oil <- gam(RepRateFall_mean ~ s(MDperKMsqFall_mean, bs="cs") + s(WellDen_mean, by=macrounit, bs="cs") + macrounit, data=AllMeans)#2
gam_oil <- gam(RepRateFall_mean ~ s(MDperKMsqFall_mean, by=macrounit, bs="cs") + s(WellDen_mean, bs="cs") + macrounit, data=AllMeans)#3

gam_oilpred <- data.frame(MDperKMsqFall_mean=AllMeans$MDperKMsqFall_mean, macrounit=AllMeans$macrounit, WellDen_meann=AllMeans$WellDen_mean)
gam_oilpred <- cbind(gam_oilpred, predict(gam_oil, se.fit=T, newdata=data.frame("MDperKMsqFall_mean"=AllMeans$MDperKMsqFall_mean, "macrounit"=AllMeans$macrounit, "WellDen_mean"=AllMeans$WellDen_mean), type="response"))
macrounitplots(glmobject = gam_oilpred,xcol="WellDen_mean",title="gam_oil fall - effect of Oil Well Density",colour="red")
macrounitplots(glmobject = gam_oilpred,title="gam_oil fall - effect of Oil Well Density",colour="red")
summary(gam_oil)
AIC(gam_oil)

gam_oilres <- residuals(gam_oil, type = "deviance")
plot(gam_oilres ~AllMeans$MDperKMsqFall_mean[which(!is.na(AllMeans$RepRateFall_mean))]) #
acf(gam_oilres, na.action = na.pass,main = "Auto-correlation plot for residuals Oil Well Density")

gam.check(gam_oil)


###Effect of Coyote Density on each macrounit

gam_coyote <- gam(RepRateFall_mean ~ s(CoyoteDen_mean, bs="cs"), data=AllMeans)
gam_coyote <- gam(RepRateFall_mean ~ s(MDperKMsqFall_mean, bs="cs") + s(CoyoteDen_mean, by=macrounit, bs="cs") + macrounit, data=AllMeans)
gam_coyote <- gam(RepRateFall_mean ~ s(MDperKMsqFall_mean, by=macrounit, bs="cs") + s(CoyoteDen_mean, bs="cs") + macrounit, data=AllMeans)

gam_coyotepred <- data.frame(MDperKMsqFall_mean=AllMeans$MDperKMsqFall_mean, macrounit=AllMeans$macrounit, CoyoteDen_meann=AllMeans$CoyoteDen_mean)
gam_coyotepred <- cbind(gam_coyotepred, predict(gam_coyote, se.fit=T, newdata=data.frame("MDperKMsqFall_mean"=AllMeans$MDperKMsqFall_mean, "macrounit"=AllMeans$macrounit, "CoyoteDen_mean"=AllMeans$CoyoteDen_mean), type="response"))
macrounitplots(glmobject = gam_coyotepred,xcol="CoyoteDen_mean",title="gam_coyote fall - effect of Coyote Density",colour="red")
macrounitplots(glmobject = gam_coyotepred,title="gam_coyote fall - effect of Coyote Density",colour="red")
summary(gam_coyote)
AIC(gam_coyote)

gam_coyoteres <- residuals(gam_coyote, type = "deviance")
plot(gam_coyoteres ~AllMeans$MDperKMsqFall_mean[which(!is.na(AllMeans$RepRateFall_mean))]) #
acf(gam_coyoteres, na.action = na.pass,main = "Auto-correlation plot for residuals Coyote Density")

gam.check(gam_coyote)



###Effect of Woody Vegetation on each macrounit


gam_woodyveg <- gam(RepRateFall_mean ~ s(WoodyVeg_mean, bs="cs"), data=AllMeans)
gam_woodyveg <- gam(RepRateFall_mean ~ s(MDperKMsqFall_mean, bs="cs") + s(WoodyVeg_mean, by=macrounit, bs="cs") + macrounit, data=AllMeans)
gam_woodyveg <- gam(RepRateFall_mean ~ s(MDperKMsqFall_mean, by=macrounit, bs="cs") + s(WoodyVeg_mean, bs="cs") + macrounit, data=AllMeans)

gam_woodyvegpred <- data.frame(MDperKMsqFall_mean=AllMeans$MDperKMsqFall_mean, macrounit=AllMeans$macrounit, WoodyVeg_meann=AllMeans$WoodyVeg_mean)
gam_woodyvegpred <- cbind(gam_woodyvegpred, predict(gam_woodyveg, se.fit=T, newdata=data.frame("MDperKMsqFall_mean"=AllMeans$MDperKMsqFall_mean, "macrounit"=AllMeans$macrounit, "WoodyVeg_mean"=AllMeans$WoodyVeg_mean), type="response"))
macrounitplots(glmobject = gam_woodyvegpred,xcol="WoodyVeg_mean",title="gam_woodyveg fall - effect of Woody Vegetation",colour="red")
macrounitplots(glmobject = gam_woodyvegpred,title="gam_woodyveg fall - effect of Woody Vegetation",colour="red")
summary(gam_woodyveg)
AIC(gam_woodyveg)
plot(gam_woodyveg)
gam_woodyvegres <- residuals(gam_woodyveg, type = "deviance")
plot(gam_woodyvegres ~AllMeans$MDperKMsqFall_mean[which(!is.na(AllMeans$RepRateFall_mean))]) #
acf(gam_woodyvegres, na.action = na.pass,main = "Auto-correlation plot for residuals Woody Vegetation")

gam.check(gam_woodyveg)



#Combined model
gam_combine <- gam(RepRateFall_mean ~ s(MDperKMsqFall_mean, bs="cs") + s(AvrgWinterMinTemp, bs="cs") + s(HuntDen_All_mean, bs="cs") + s(WellDen_mean, bs="cs") +  s(CoyoteDen_mean, bs="cs")+s(WoodyVeg_mean, bs="cs"), data=AllMeans)
summary(gam_combine)
AIC(gam_combine)
