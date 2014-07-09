library(lattice)
# effect of WinterTemp on Fawn Female Ratio as in Simones model?

gam_temp <- gam(FawnFemaleRatio_mean ~ s(AvrgWinterMinTemp, bs="cs"), data=AllMeans)
gam_temp <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +s(AvrgWinterMinTemp, bs="cs"), data=AllMeans)
summary(gam_temp)
plot(gam_temp)
AIC(gam_temp) #YES
gam_temppred<-predict(gam_temp, se.fit=T, newdata=data.frame(year=AllMeans$year, AvrgWinterMinTemp=AllMeans$AvrgWinterMinTemp), type="response")
plot(AllMeans$FawnFemaleRatio_mean~AllMeans$year, cex=0.5, main="GAM Fawn Female Ratio ~ Year and Winter Temp")
lines(gam_temppred$fit[1:51]~AllMeans$year[1:51],type="l", col="red", )


# Winter Mortality: Effects of Hunt and Temp
plot(AllMeans$WinterMort_mean~AllMeans$year) # 2 outliers but removing them doesnt change a lot (see "leftover code.R")
plot(AllMeans$WinterMort_mean~AllMeans$HuntDen_All_mean, log="y")
plot(AllMeans$WinterMort_mean~AllMeans$AvrgWinterMinTemp, log="y")
# both dont show a trend -> probably no significant effect

gam_temp <- gam(WinterMort_mean ~ s(AvrgWinterMinTemp,  bs="cs"), data=AllMeans)#NO EFFECT
gam_hunt <- gam(WinterMort_mean ~ s(HuntDen_All_mean,  bs="cs"), data=AllMeans)#NO EFFECT
plot(gam_hunt)
summary(gam_temp)
summary(gam_hunt)
plot(gam_hunt)
AIC(gam_temp)


# Relationship Spring+Fall Data
cor(AllMeans$MDperKMsqSpring_mean, AllMeans$MDperKMsqFall_mean, use="pairwise.complete.obs")#0.71
cor(Spring_tplus1, AllMeans$MDperKMsqFall_mean, use="pairwise.complete.obs")#0.64 -> does not make sense

cor(AllMeans$FawnFemaleRatio_mean, AllMeans$MDperKMsqSpring_mean_tplus1,use="pairwise.complete.obs")#-0.4

#Hunting licenses issued by reproduction rate/population trend
cor(AllMeans$FawnFemaleRatio_mean,AllMeans$HuntDen_All_mean_tplus1,use="pairwise.complete.obs") #-0.05
cor(AllMeans$RepRateFall_mean,AllMeans$HuntDen_All_mean_tplus1, use="pairwise.complete.obs") # -0.11
# -> does not agree with infos on license issuance according to poptrends
# however, this only affects issues for antlerless deer
cor(HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-1")],HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-2")], use="pairwise.complete.obs")
cor(HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-1")],HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-3")], use="pairwise.complete.obs")
cor(HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-1")],HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-4")], use="pairwise.complete.obs")
cor(HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-2")],HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-3")], use="pairwise.complete.obs")
cor(HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-3")],HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-4")], use="pairwise.complete.obs")

# antlerless deer?
count_nas(AllMeans$HuntDen_Aless_mean) #32
count_nas(mules1962$d3) # 0
count_nas(mules1962$d2) # 184 out of 1173 !

plot(WholeAreaMeans$HuntDen_All_mean~WholeAreaMeans$year, type="l",lwd=1, main = "Hunting Densities Whole Area", ylab="Hunting Density per km²", xlab="Year")
lines(WholeAreaMeans$HuntDen_Aless_mean~WholeAreaMeans$year, col="red")#lots of NAs in raw data
legend("topleft", legend=c("Total", "Antlerless"), col=c(1, "red"), lty=1)

# total_harvest.MDgunHarvest.
# 
# antlerless <- mules1962$Antleress_Harvest.MDgunHarvest.
# all <- mules1962$total_harvest.MDgunHarvest.
# mu <- mules1962
# year <- 
# xyplot(Antleress_Harvest.MDgunHarvest.|macrounit,data=mules1962)
# lines(WholeAreaMeans$HuntDen_Aless_mean~WholeAreaMeans$year, col="red")#lots of NAs in raw data
# legend("topleft", legend=c("Total", "Antlerless"), col=c(1, "red"), lty=1)
# 
# #cor(WholeAreaMeans$HuntDen_All_mean,WholeAreaMeans$HuntDen_Aless_mean,use="pairwise.complete.obs")# 0.95
# #cor(AllMeans$FawnFemaleRatio_mean,AllMeans$HuntDen_Aless_mean_tplus1,use="pairwise.complete.obs") #-0.05
#cor(AllMeans$RepRateFall_mean,AllMeans$HuntDen_Aless_mean_tplus1, use="pairwise.complete.obs") # -0.11
