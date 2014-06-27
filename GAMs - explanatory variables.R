###------GAMs of fall MD population density including explanatory variables, no autocorrelation


# Data Analysis of data that is added for this step
count_nas(AllMeans$AvrgWinterMinTemp)#0, ok
count_nas(AllMeans$FawnFall_mean)#7
count_nas(AllMeans$FemaleFall_mean)#7
count_nas(AllMeans$MaleFall_mean)#7
count_nas(AllMeans$RepRateFall_mean)#7


####Check for outliers of data that is added for this step
plot(AllMeans$AvrgWinterMinTemp)


#Check for collinearity between predictors
source("HighstatLibV6_correlation_functions.R") # Replacement for AEV-package of Zuur et al
z <- AllMeans[,!(names(AllMeans) %in% c("macrounit", "MDperKMsqSpring_mean","WTailDen_mean","FawnFall_mean","MaleFall_mean","FemaleFall_mean"))] 

par(oma=c(2,0,2,0))
pairs(z, lower.panel = panel.smooth2, upper.panel = panel.cor, diag.panel = panel.hist) #not too informative
title("Pairwise Pearson Correlation", outer=TRUE)

cormatrix <- cor(z, use="pairwise.complete.obs") # no correlation >(-)0.4 so ok
write(cormatrix, "correlation matrix predictors.txt", sep="\t")
vif <- corvif(z)# all values <10, so ok



###Effect of Average Minimum Winter Temperature on each macrounit
gam_temp <- gam(MDperKMsqFall_mean ~ s(year, bs="cs") +s(AvrgWinterMinTemp, bs="cs"), data=AllMeans)#1
gam_temp <- gam(MDperKMsqFall_mean ~ s(year, bs="cs") + s(AvrgWinterMinTemp, by=macrounit, bs="cs") + macrounit, data=AllMeans)#2
gam_temp <- gam(MDperKMsqFall_mean ~ s(year, by=macrounit, bs="cs") +s(AvrgWinterMinTemp, bs="cs") + macrounit, data=AllMeans)#3


gam_temppred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit, AvrgWinterMinTemp=AllMeans$AvrgWinterMinTemp)
gam_temppred <- cbind(gam_temppred, predict(gam_temp, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit, "AvrgWinterMinTemp"=AllMeans$AvrgWinterMinTemp), type="response"))
macrounitplots(glmobject = gam_temppred,xcol="AvrgWinterMinTemp",title="gam_temp fall - effect of Mean Winter Minimum Temperature",colour="red")
macrounitplots(glmobject = gam_temppred,title="gam_temp fall - effect of Mean Winter Minimum Temperature",colour="red")
summary(gam_temp)
AIC(gam_temp)
par(oma=c(2,0,2,0))
gam.check(gam_temp)
# title("Gam_all2 fall residual check", outer=TRUE)
# #comparison nullmodell gam_all09 and gam_all2
# anova(gam_all0, gam_all2, test="F")

# Check Residuals for temporal autocorrelation
gam_tempres <- residuals(gam_temp, type = "deviance")
plot(gam_tempres ~AllMeans$year[which(!is.na(AllMeans$MDperKMsqFall_mean))]) #
acf(gam_tempres, na.action = na.pass,main = "Auto-correlation plot for residuals gam_temp fall")



###Effect of hunting Density on each macrounit
gam_hunt <- gam(MDperKMsqFall_mean ~ s(year, bs="cs") + s(HuntDen_All_mean, bs="cs"), data=AllMeans)#1
gam_hunt <- gam(MDperKMsqFall_mean ~ s(year, bs="cs") + s(HuntDen_All_mean, by=macrounit, bs="cs") + macrounit, data=AllMeans)#2
gam_hunt <- gam(MDperKMsqFall_mean ~ s(year, by=macrounit, bs="cs") + s(HuntDen_All_mean, bs="cs") + macrounit, data=AllMeans)#3

gam_huntpred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit, HuntDen_All_mean=AllMeans$HuntDen_All_mean)
gam_huntpred <- cbind(gam_huntpred, predict(gam_hunt, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit, "HuntDen_All_mean"=AllMeans$HuntDen_All_mean), type="response"))
macrounitplots(glmobject = gam_huntpred,xcol="HuntDen_All_mean",title="gam_hunt fall - effect of Hunting Density",colour="red")
macrounitplots(glmobject = gam_huntpred,title="gam_hunt fall - effect of Hunting Density",colour="red")
summary(gam_hunt)
AIC(gam_hunt)

gam_huntres <- residuals(gam_hunt, type = "deviance")
plot(gam_huntres ~AllMeans$year[which(!is.na(AllMeans$MDperKMsqFall_mean))]) #
acf(gam_huntres, na.action = na.pass,main = "Auto-correlation plot for residuals gam_hunt fall")

gam.check(gam_hunt)


###Effect of Oil Well Density on each macrounit


gam_oil <- gam(MDperKMsqFall_mean ~ s(year, bs="cs") + s(WellDen_mean, bs="cs"), data=AllMeans)#1
gam_oil <- gam(MDperKMsqFall_mean ~ s(year, bs="cs") + s(WellDen_mean, by=macrounit, bs="cs") + macrounit, data=AllMeans)#2
gam_oil <- gam(MDperKMsqFall_mean ~ s(year, by=macrounit, bs="cs") + s(WellDen_mean, bs="cs") + macrounit, data=AllMeans)#3

gam_oilpred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit, WellDen_meann=AllMeans$WellDen_mean)
gam_oilpred <- cbind(gam_oilpred, predict(gam_oil, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit, "WellDen_mean"=AllMeans$WellDen_mean), type="response"))
macrounitplots(glmobject = gam_oilpred,xcol="WellDen_mean",title="gam_oil fall - effect of Oil Well Density",colour="red")
macrounitplots(glmobject = gam_oilpred,title="gam_oil fall - effect of Oil Well Density",colour="red")
summary(gam_oil)
AIC(gam_oil)

gam_oilres <- residuals(gam_oil, type = "deviance")
plot(gam_oilres ~AllMeans$year[which(!is.na(AllMeans$MDperKMsqFall_mean))]) #
acf(gam_oilres, na.action = na.pass,main = "Auto-correlation plot for residuals Oil Well Density")

gam.check(gam_oil)


###Effect of Coyote Density on each macrounit

gam_coyote <- gam(MDperKMsqFall_mean ~ s(year, bs="cs") + s(CoyoteDen_mean, bs="cs"), data=AllMeans)
gam_coyote <- gam(MDperKMsqFall_mean ~ s(year, bs="cs") + s(CoyoteDen_mean, by=macrounit, bs="cs") + macrounit, data=AllMeans)
gam_coyote <- gam(MDperKMsqFall_mean ~ s(year, by=macrounit, bs="cs") + s(CoyoteDen_mean, bs="cs") + macrounit, data=AllMeans)

gam_coyotepred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit, CoyoteDen_meann=AllMeans$CoyoteDen_mean)
gam_coyotepred <- cbind(gam_coyotepred, predict(gam_coyote, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit, "CoyoteDen_mean"=AllMeans$CoyoteDen_mean), type="response"))
macrounitplots(glmobject = gam_coyotepred,xcol="CoyoteDen_mean",title="gam_coyote fall - effect of Coyote Density",colour="red")
macrounitplots(glmobject = gam_coyotepred,title="gam_coyote fall - effect of Coyote Density",colour="red")
summary(gam_coyote)
AIC(gam_coyote)

gam_coyoteres <- residuals(gam_coyote, type = "deviance")
plot(gam_coyoteres ~AllMeans$year[which(!is.na(AllMeans$MDperKMsqFall_mean))]) #
acf(gam_coyoteres, na.action = na.pass,main = "Auto-correlation plot for residuals Coyote Density")

gam.check(gam_coyote)



###Effect of Woody Vegetation on each macrounit

gam_woodyveg <- gam(MDperKMsqFall_mean ~ s(year, bs="cs") + s(WoodyVeg_mean, bs="cs"), data=AllMeans)
gam_woodyveg <- gam(MDperKMsqFall_mean ~ s(year, bs="cs") + s(WoodyVeg_mean, by=macrounit, bs="cs") + macrounit, data=AllMeans)
gam_woodyveg <- gam(MDperKMsqFall_mean ~ s(year, by=macrounit, bs="cs") + s(WoodyVeg_mean, bs="cs") + macrounit, data=AllMeans)

gam_woodyvegpred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit, WoodyVeg_meann=AllMeans$WoodyVeg_mean)
gam_woodyvegpred <- cbind(gam_woodyvegpred, predict(gam_woodyveg, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit, "WoodyVeg_mean"=AllMeans$WoodyVeg_mean), type="response"))
macrounitplots(glmobject = gam_woodyvegpred,xcol="WoodyVeg_mean",title="gam_woodyveg fall - effect of Woody Vegetation",colour="red")
macrounitplots(glmobject = gam_woodyvegpred,title="gam_woodyveg fall - effect of Woody Vegetation",colour="red")
summary(gam_woodyveg)
AIC(gam_woodyveg)

gam_woodyvegres <- residuals(gam_woodyveg, type = "deviance")
plot(gam_woodyvegres ~AllMeans$year[which(!is.na(AllMeans$MDperKMsqFall_mean))]) #
acf(gam_woodyvegres, na.action = na.pass,main = "Auto-correlation plot for residuals Woody Vegetation")

gam.check(gam_woodyveg)

##### All explanatory variables, Whole Area Means (equivalent to gam_all3)
gam_combine <- gam(MDperKMsqFall_mean ~ s(year, bs="cs") + s(AvrgWinterMinTemp, bs="cs") + s(HuntDen_All_mean, bs="cs") + s(WellDen_mean, bs="cs") + s(CoyoteDen_mean, bs="cs") + s(WoodyVeg_mean, bs="cs"), data=AllMeans)
gam_combine <- gam(MDperKMsqFall_mean ~ s(year, bs="cs") + s(WellDen_mean, bs="cs") + s(CoyoteDen_mean, bs="cs") + s(WoodyVeg_mean, bs="cs"), data=AllMeans)


#gam_combinepred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit, CoyoteDen_meann=AllMeans$CoyoteDen_mean)
#gam_combinepred <- cbind(gam_combinepred, predict(gam_combine, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit, "CoyoteDen_mean"=AllMeans$CoyoteDen_mean), type="response"))
#macrounitplots(glmobject = gam_combinepred,xcol="CoyoteDen_mean",title="gam_combine fall - effect of Coyote Density",colour="red")
#macrounitplots(glmobject = gam_combinepred,title="gam_combine fall - effect of Coyote Density",colour="red")
summary(gam_combine)
AIC(gam_combine)

gam_combineres <- residuals(gam_combine, type = "deviance")
#plot(gam_combineres ~AllMeans$year[which(!is.na(AllMeans$MDperKMsqFall_mean))]) #
acf(gam_combineres, na.action = na.pass,main = "Auto-correlation plot for residuals Coyote Density")

gam.check(gam_combine)


#### combined gam based on Whole Area Means (equivalent go gam3)

gam_combine3 <- gam(MDperKMsqFall_mean ~ s(year, bs="cs") + s(AvrgWinterMinTemp, bs="cs") + s(HuntDen_All_mean, bs="cs") + s(WellDen_mean, bs="cs") +  s(WoodyVeg_mean, bs="cs"), data=WholeAreaMeans)

#gam_combinepred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit, CoyoteDen_meann=AllMeans$CoyoteDen_mean)
#gam_combinepred <- cbind(gam_combinepred, predict(gam_combine, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit, "CoyoteDen_mean"=AllMeans$CoyoteDen_mean), type="response"))
#macrounitplots(glmobject = gam_combinepred,xcol="CoyoteDen_mean",title="gam_combine fall - effect of Coyote Density",colour="red")
#macrounitplots(glmobject = gam_combinepred,title="gam_combine fall - effect of Coyote Density",colour="red")
summary(gam_combine3)
AIC(gam_combine3)

gam_combine3res <- residuals(gam_combine3, type = "deviance")
#plot(gam_combine3res ~AllMeans$year[which(!is.na(AllMeans$MDperKMsqFall_mean))]) #
acf(gam_combine3res, na.action = na.pass,main = "Auto-correlation plot for residuals Coyote Density")

gam.check(gam_combine)