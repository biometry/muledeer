###------GAMs of fall MD population density including explanatory variables, no autocorrelation

# As this analysis is supposed to result in information how to start with a population density model, some of the data were excluded/modified/recalculated

library(mgcv)
library(reshape2)
library(car)
source("Supplementary_Functions.R")

allmules <- read.csv("muledeer_final_dataset.csv")

# Only use the time series starting in 1962 and drop 3 study sites that were scarcely-surveyed:
mules1962 <- subset(allmules, year>1961)
del <- which(mules1962$StudyArea == "Bible_Camp" | mules1962$StudyArea == "NUTRNP"| mules1962$StudyArea == "SUTRNP")
mules1962 <- mules1962[-del,]

#Data Extraction (more data included than in other R-Scripts!)

# Data Analysis of data that is added for this step
count_nas(AllMeans$AvrgWinterMinTemp)#0, ok
count_nas(AllMeans$FawnFall_mean)#7
count_nas(AllMeans$FemaleFall_mean)#7
count_nas(AllMeans$MaleFall_mean)#7
count_nas(AllMeans$RepRateFall_mean)#7

#Check for outliers of data that is added for this step
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
gam_temp <- gam(MDperKMsqFall_mean ~ s(AvrgWinterMinTemp, by=macrounit) + macrounit, data=AllMeans)
#second try (including year:temp interaction as well)
gam_temp <- gam(MDperKMsqFall_mean ~ s(year, by=AvrgWinterMinTemp) + s(AvrgWinterMinTemp, by=macrounit) + macrounit, data=AllMeans)

gam_temppred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit, AvrgWinterMinTemp=AllMeans$AvrgWinterMinTemp)
gam_temppred <- cbind(gam_temppred, predict(gam_temp, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit, "AvrgWinterMinTemp"=AllMeans$AvrgWinterMinTemp), type="response"))
macrounitplots(glmobject = gam_temppred,xcol="AvrgWinterMinTemp",title="gam_temp fall - effect of Mean Winter Minimum Temperature",colour="red")
macrounitplots(glmobject = gam_temppred,title="gam_temp fall - effect of Mean Winter Minimum Temperature",colour="red")
summary(gam_temp)
AIC(gam_temp)
# par(oma=c(2,0,2,0))
# gam.check(gam_temp)
# title("Gam_all2 fall residual check", outer=TRUE)
# #comparison nullmodell gam_all09 and gam_all2
# anova(gam_all0, gam_all2, test="F")

# Check Residuals for temporal autocorrelation
gam_tempres <- residuals(gam_temp, type = "deviance")
plot(gam_tempres ~AllMeans$year[which(!is.na(AllMeans$MDperKMsqFall_mean))]) #
acf(gam_tempres, na.action = na.pass,main = "Auto-correlation plot for residuals gam_temp fall")


macrounitplots(glmobject = gam_all2pred,title="GAM2 fall - interaction year:MDdensity",colour="red")

summary(gam_all2)
par(oma=c(2,0,2,0))
gam.check(gam_all2)#residual show clear patterns -> normal distribution?variance homogenity?
title("Gam_all2 fall residual check", outer=TRUE)

gam_all3 <- gam(MDperKMsqFall_mean ~ s(year, by=macrounit) + te(AvrgWinterMinTemp,CoyoteDen_mean, bs="cs"), data=AllMeans)
gam_all3pred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit)
gam_all3pred <- cbind(gam_all3pred, predict(gam_all3, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit, "CoyoteDen_mean"=AllMeans$CoyoteDen_mean,"AvrgWinterMinTemp"= AllMeans$AvrgWinterMinTemp), type="response"))
macrounitplots(glmobject = gam_all3pred,title="GAM3 fall - interaction Coyote+Temp",colour="red")

summary(gam_all4$gam)

gam_all4 <- gamm(MDperKMsqFall_mean ~ s(year, by=macrounit) + te(AvrgWinterMinTemp,CoyoteDen_mean, bs="cs"), correlation = corAR1(form=~year|macrounit), data=AllMeans)
vis.gam(gam_all4$gam)
AIC(gam_all4$lme)

gam_all4b <- gamm(MDperKMsqFall_mean ~ s(year, by=macrounit) + s(CoyoteDen_mean, bs="cs") + s(AvrgWinterMinTemp), correlation = corAR1(form=~year|macrounit), data=AllMeans)
summary(gam_all4b$lme)
summary(gam_all4b$gam)
plot(gam_all4b$gam)
AIC(gam_all4b$lme)


gam_all5 <- gamm(MDperKMsqFall_mean ~ s(AvrgWinterMinTemp, bs="cs"), correlation = corAR1(form=~year|macrounit), data=AllMeans)
summary(gam_all5$lme)
plot(gam_all5$gam)

gam_all5 <- gam(MDperKMsqFall_mean ~ s(year, by=macrounit)+ s(AvrgWinterMinTemp, bs="cs"), data=AllMeans)
summary(gam_all5)
plot(gam_all5)
AIC(gam_all5)


plot(gam_all4$gam)
par(mfrow=c(1,1))
vis.gam(gam_all2)
vis.gam(gam_all3)
###Effect of hunting Density on each macrounit
#gam_hunt <- gam(MDperKMsqFall_mean ~ s(HuntDen_All_mean, by=macrounit) + macrounit, data=AllMeans)
#second try(interaction year:hunt):
gam_hunt <- gam(MDperKMsqFall_mean ~ s(year, by=HuntDen_All_mean) + s(HuntDen_All_mean, by=macrounit) + macrounit, data=AllMeans)
gam_huntpred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit, HuntDen_All_mean=AllMeans$HuntDen_All_mean)
gam_huntpred <- cbind(gam_huntpred, predict(gam_hunt, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit, "HuntDen_All_mean"=AllMeans$HuntDen_All_mean), type="response"))
macrounitplots(glmobject = gam_huntpred,xcol="HuntDen_All_mean",title="gam_hunt fall - effect of Hunting Density",colour="red")
macrounitplots(glmobject = gam_huntpred,title="gam_hunt fall - effect of Hunting Density",colour="red")
summary(gam_hunt)
AIC(gam_hunt)

gam_huntres <- residuals(gam_hunt, type = "deviance")
plot(gam_huntres ~AllMeans$year[which(!is.na(AllMeans$MDperKMsqFall_mean))]) #
acf(gam_huntres, na.action = na.pass,main = "Auto-correlation plot for residuals gam_hunt fall")

###Effect of Oil Well Density on each macrounit
#gam_hunt <- gam(MDperKMsqFall_mean ~ s(HuntDen_All_mean, by=macrounit) + macrounit, data=AllMeans)
#second try(interaction year:oil):
gam_oil <- gam(MDperKMsqFall_mean ~ s(year, by=WellDen_mean) + s(WellDen_mean, by=macrounit) + macrounit, data=AllMeans)
gam_oilpred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit, WellDen_meann=AllMeans$WellDen_mean)
gam_oilpred <- cbind(gam_oilpred, predict(gam_oil, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit, "WellDen_mean"=AllMeans$WellDen_mean), type="response"))
macrounitplots(glmobject = gam_oilpred,xcol="WellDen_mean",title="gam_oil fall - effect of Oil Well Density",colour="red")
macrounitplots(glmobject = gam_oilpred,title="gam_oil fall - effect of Oil Well Density",colour="red")
summary(gam_oil)
AIC(gam_oil)

gam_oilres <- residuals(gam_oil, type = "deviance")
plot(gam_oilres ~AllMeans$year[which(!is.na(AllMeans$MDperKMsqFall_mean))]) #
acf(gam_oilres, na.action = na.pass,main = "Auto-correlation plot for residuals Oil Well Density")


### Stepwise AIC selection of the best model? not possible...
gam <- gam(MDperKMsqFall_mean ~ s(year, by=macrounit) + s(AvrgWinterMinTemp, by=macrounit) + s(HuntDen_All_mean, by=macrounit) + s(WellDen_mean, by=macrounit) + macrounit, data=AllMeans)
gam_pred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit, WellDen_meann=AllMeans$WellDen_mean,AvrgWinterMinTemp=AllMeans$AvrgWinterMinTemp, HuntDen_All_mean=AllMeans$HuntDen_All_mean)
gam_pred <- cbind(gam_pred, predict(gam, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit, "WellDen_mean"=AllMeans$WellDen_mean,"AvrgWinterMinTemp"=AllMeans$AvrgWinterMinTemp, "HuntDen_All_mean"=AllMeans$HuntDen_All_mean), type="response"))
macrounitplots(glmobject = gam_pred,title="gam_stepwise",colour="red")
summary(gam)
library(MASS)
gam_step <- stepAIC(gam)#does not work
AIC(gam)
