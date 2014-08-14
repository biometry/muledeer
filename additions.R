AllMeans <- read.csv("AllMeans.csv")
source("Supplementary_Functions.R")

#+++data preparation wrong names (plus insted of minus) ---------------

MDperKMsqSpring_mean_tminus1 <- c(NA, AllMeans$MDperKMsqSpring_mean[-length(AllMeans$MDperKMsqSpring_mean)])
AllMeans <- cbind(AllMeans,MDperKMsqSpring_mean_tminus1)

AllMeans <- cbind(AllMeans,AvrgWinterMinTemp_tminus1=c(NA, AllMeans$AvrgWinterMinTemp[-length(AllMeans$AvrgWinterMinTemp)]))

# timelags hunting
AllMeans <- cbind(AllMeans,HuntDen_All_mean_tminus1=c(NA, AllMeans$HuntDen_All_mean[-length(AllMeans$HuntDen_All_mean)]))
AllMeans <- cbind(AllMeans,HuntDen_Aless_mean_tminus1=c(NA, AllMeans$HuntDen_Aless_mean[-length(AllMeans$HuntDen_Aless_mean)]))



#FawnFemaleRatio time lag 
AllMeans <- cbind(AllMeans,FawnFemaleRatio_mean_tminus1 = c(NA, AllMeans$FawnFemaleRatio_mean[-length(AllMeans$FawnFemaleRatio_mean)]))



#+++secondary data exploration:Correlation------------

# Collinearity between PopDen in MUs (as to see wether they differ much)

z <- data.frame(AllMeans[which(AllMeans$macrounit == "0-1"),"MDperKMsqSpring_mean"])
z <- cbind(z, AllMeans[which(AllMeans$macrounit == "0-2"),"MDperKMsqSpring_mean"],AllMeans[which(AllMeans$macrounit == "0-3"),"MDperKMsqSpring_mean"],AllMeans[which(AllMeans$macrounit == "0-4"),"MDperKMsqSpring_mean"])
names(z) <- c("0-1","0-2","0-3","0-4")
pairs(z, lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearities between macrounits")


#### Hunting licences correlation Hunt <-> RepRate(t+1) and Fawn:Female Ratio as information on distribution practice of licenses---------
z <- AllMeans[,names(AllMeans) %in% c("HuntDen_All_mean", "RepRateFall_mean_tminus1", "FawnFemaleRatio_mean_tminus1", "WoodyVeg_mean", "WellDen_mean", "AvrgWinterMinTemp_tminus1", "CoyoteDen_mean")] 

pairs(z, lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Whole Area")
pairs(z[which(AllMeans$macrounit == "0-1"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 1")
pairs(z[which(AllMeans$macrounit == "0-2"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 2")
pairs(z[which(AllMeans$macrounit == "0-3"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 3")
pairs(z[which(AllMeans$macrounit == "0-4"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 4")

#Estimates of Population Densities and growth------------
z <- AllMeans[,names(AllMeans) %in% c("MDperKMsqFall_mean", "MDperKMsqSpring_mean", "RepRateFall_mean_tminus1", "FawnFemaleRatio_mean_tminus1")] 

pairs(z, lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Whole Area")
pairs(z[which(AllMeans$macrounit == "0-1"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 1")
pairs(z[which(AllMeans$macrounit == "0-2"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 2")
pairs(z[which(AllMeans$macrounit == "0-3"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 3")
pairs(z[which(AllMeans$macrounit == "0-4"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 4")


# +++Winter effects ----
### Winter effects on Fawn Female Ratio (coyote,Temp and Hunt and interactions) - same as in Simones model?---------------

#without Coyote
gam_temp <- gam(FawnFemaleRatio_mean ~ s(AvrgWinterMinTemp, bs="cs"), data=AllMeans)#yes
gam_temp <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +s(AvrgWinterMinTemp, bs="cs"), data=AllMeans)#yes
gam_temp <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +s(AvrgWinterMinTemp, by=macrounit, bs="cs")+macrounit, data=AllMeans)#only MU1
gam_temp <- gam(FawnFemaleRatio_mean ~ s(AvrgWinterMinTemp, by=macrounit, bs="cs")+macrounit, data=AllMeans)#only MU1

#Temp*Coyote
gam_coyotettemp <- gam(FawnFemaleRatio_mean ~ te(CoyoteDen_mean,AvrgWinterMinTemp,bs="cs"), data=AllMeans)#yes
gam_coyotettemp <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +te(CoyoteDen_mean,AvrgWinterMinTemp,bs="cs"), data=AllMeans)#yes
gam_coyotettemp <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +te(CoyoteDen_mean,AvrgWinterMinTemp,bs="cs", by=macrounit)+macrounit, data=AllMeans)#None
gam_coyotettemp <- gam(FawnFemaleRatio_mean ~ te(CoyoteDen_mean,AvrgWinterMinTemp,bs="cs", by=macrounit)+macrounit, data=AllMeans)#only MU1

# Hunt only
gam_hunt <- gam(FawnFemaleRatio_mean ~ s(HuntDen_All_mean, bs="cs"), data=AllMeans)#nope
gam_hunt <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +s(HuntDen_All_mean, bs="cs"), data=AllMeans)#no
gam_hunt <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +s(HuntDen_All_mean, by=macrounit, bs="cs")+macrounit, data=AllMeans)#no
gam_hunt <- gam(FawnFemaleRatio_mean ~ s(HuntDen_All_mean, by=macrounit, bs="cs")+macrounit, data=AllMeans)#only MU1

#Temp*Hunt
gam_hunttemp <- gam(FawnFemaleRatio_mean ~ te(HuntDen_All_mean,AvrgWinterMinTemp,bs="cs"), data=AllMeans)#yes!!
gam_hunttemp <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +te(HuntDen_All_mean,AvrgWinterMinTemp,bs="cs"), data=AllMeans)#somewhat:0.0354 *
gam_hunttemp <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +te(HuntDen_All_mean,AvrgWinterMinTemp,bs="cs", by=macrounit)+macrounit, data=AllMeans)#MU1+MU4
gam_hunttemp <- gam(FawnFemaleRatio_mean ~ te(HuntDen_All_mean,AvrgWinterMinTemp,bs="cs", by=macrounit)+macrounit, data=AllMeans)#only MU1:0.00456 **

summary(gam_temp)
summary(gam_hunt)
summary(gam_hunttemp)
AIC(gam_temp)
plot(gam_temp)
gam_temppred<-predict(gam_temp, se.fit=T, newdata=data.frame(year=AllMeans$year, AvrgWinterMinTemp=AllMeans$AvrgWinterMinTemp), type="response")
plot(AllMeans$FawnFemaleRatio_mean~AllMeans$year, cex=0.5, main="GAM Fawn Female Ratio ~ Year and Winter Temp")
lines(gam_temppred$fit[1:51]~AllMeans$year[1:51],type="l", col="red", )

### Winter effects on Winter mortality

gam_temp <- gam(WinterMort_mean ~ s(AvrgWinterMinTemp_tplus1,  bs="cs"), data=AllMeans)#NO EFFECT
gam_hunt <- gam(WinterMort_mean ~ s(HuntDen_All_mean,  bs="cs"), data=AllMeans)#NO EFFECT
gam_hunttemp <- gam(WinterMort_mean ~ te(AvrgWinterMinTemp_tplus1, HuntDen_All_mean,  bs="cs"), data=AllMeans)#NO EFFECT
gam_hunttemp <- gam(WinterMort_mean ~ te(AvrgWinterMinTemp_tplus1, HuntDen_All_mean,  bs="cs",by=macrounit)+macrounit, data=AllMeans)# MU4!


# +++ GAMS PopDen ----
###Effect of hunting Density on each macrounit (timelag t-1 as harvest starts after fall survey!)
gam_hunt <- gam(MDperKMsqFall_mean ~ s(year, bs="cs") + s(HuntDen_All_mean_tminus1, bs="cs"), data=AllMeans)#1
gam_hunt <- gam(MDperKMsqFall_mean ~ s(year, bs="cs") + s(HuntDen_All_mean_tminus1, by=macrounit, bs="cs") + macrounit, data=AllMeans)#2
gam_hunt <- gam(MDperKMsqFall_mean ~ s(year, by=macrounit, bs="cs") + s(HuntDen_All_mean_tminus1, bs="cs") + macrounit, data=AllMeans)#3

summary(gam_hunt)

# To Do: same in interactions

# +++ GAMS Reproduction Rate-----

#interaction temp+hunt (= the two sole winter effects) -> useless
gam_hunttemp <- gam(RepRateFall_mean ~ te(HuntDen_All_mean,AvrgWinterMinTemp,bs="cs"), data=AllMeans)#no
gam_hunttemp <- gam(RepRateFall_mean ~ s(MDperKMsqFall_mean, bs="cs") + te(HuntDen_All_mean,AvrgWinterMinTemp, bs="cs", by=macrounit) + macrounit, data=AllMeans)#no
gam_hunttemp <- gam(RepRateFall_mean ~ s(MDperKMsqFall_mean, by=macrounit, bs="cs") +te(HuntDen_All_mean,AvrgWinterMinTemp, bs="cs") + macrounit, data=AllMeans)#no

