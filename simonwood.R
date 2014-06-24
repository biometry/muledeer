library(mgcv)
AllMeans <- read.csv("testdata.csv")


###QUESTIONS:
# - Why do the models 1+2 fit a straight line when autocorrelation-term is included
# - Why does autocorrelation-term not improve the AIC (all models)
# - Why does correlation structure change to "ARMA(1,0)"
# - correlation term: form = ~1|year vs form = ~year|macrounit vs. form = ~year


### Models including the effect of macrounits
# no autocorrelation
gam1 <- gam(MDperKMsqFall_mean ~ s(year, by=macrounit) + macrounit, data=AllMeans)
summary(gam1)
plot(gam1) # fits data ok
acf(residuals(gam1,type = "deviance")) # autocorrelation of residuals at year-3
AIC(gam1)#338

# autocorrelation
gam_ac1 <- gamm(MDperKMsqFall_mean ~ s(year, by=macrounit) + macrounit, correlation = corAR1(form=~year|macrounit), data=AllMeans)
summary(gam_ac1$gam)
summary(gam_ac1$lme) # Correlation Structure: ARMA(1,0)
plot(gam_ac1$gam) # straight lines
acf(residuals(gam_ac1$gam,type = "deviance")) # more autocorrelation in residuals than before?
AIC(gam_ac1$lme)#391

### Models for one macrounit only
gam2 <- gam(MDperKMsqFall_mean ~ s(year), data=AllMeans[which(AllMeans$macrounit == "0-1"),])
summary(gam2)
plot(gam2) # fit is still ok
acf(residuals(gam2,type = "deviance")) #autocorrelation of residuals at year-4
AIC(gam2)#74

gam_ac2 <- gamm(MDperKMsqFall_mean ~ s(year), correlation = corAR1(form=~year), data=AllMeans[which(AllMeans$macrounit == "0-1"),])
summary(gam_ac2$gam)
summary(gam_ac2$lme) # Correlation Structure: ARMA(1,0)
plot(gam_ac2$gam) # straight line again
acf(residuals(gam_ac2$gam,type = "deviance")) # autocorrelation in residuals at year-2
AIC(gam_ac2$lme)#81


### Models for all macrounits
gam3 <- gam(MDperKMsqFall_mean ~ s(year), data=AllMeans)
summary(gam3)
plot(gam3) # fit is ok
acf(residuals(gam3,type = "deviance")) #autocorrelation throughout 20 years
AIC(gam3)#457

gam_ac3 <- gamm(MDperKMsqFall_mean ~ s(year), correlation = corAR1(form=~1|year), data=AllMeans)
summary(gam_ac3$gam)
summary(gam_ac3$lme) # Correlation Structure: AR1
plot(gam_ac3$gam) # smoother that fits data
acf(residuals(gam_ac3$gam,type = "deviance")) # exactly the the same autocorrelation
AIC(gam_ac3$lme)#468
