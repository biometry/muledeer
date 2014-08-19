library(lattice)

#recalculate population densities (see "Summary_Data.doc")
sup <- mules1962[c("year","StudyArea", "macrounit", "AreaSqKm_true_fall", "TotalMDFall")]
count_nas(sup$TotalMDFall) #350 NAs
length(which(sup$AreaSqKm_true_fall == 0)) #350 so goes along with NAs in popcounts = correct
sup_pop <- tapply(X=sup$TotalMDFall, INDEX = list(sup$year,sup$macrounit), FUN=sum, na.rm = T) #no of individuals in all surveyed studysites within each macrounit
sup_area <- tapply(X=sup$AreaSqKm_true_fall, INDEX = list(sup$year,sup$macrounit), FUN=sum, na.rm = T) # area of SURVEYED studysites in that year
sup <- sup_pop/sup_area
sup <- melt(sup)
count_nas(sup)#7
count_nas(AllMeans$MDperKMsqFall_mean)#7 so ok
AllMeans <- cbind(AllMeans, MDperKMsqFall_newmean = sup[,3])

plot(AllMeans$MDperKMsqFall_mean~AllMeans$MDperKMsqFall_newmean,cex=0.5, main="Difference in Density Calculation Methods")
lines(abs(AllMeans$MDperKMsqFall_mean-AllMeans$MDperKMsqFall_newmean)~ AllMeans$MDperKMsqFall_newmean, col="red", type="h")
legend("topleft", legend=c("Values", "Sum of Deviance"), col=c(1, "red"), lty=1)
cor(AllMeans$MDperKMsqFall_mean,AllMeans$MDperKMsqFall_newmean, use="pairwise.complete.obs")#0.9903044
xyplot(MDperKMsqFall_mean +MDperKMsqFall_newmean ~ year | macrounit, data=AllMeans, type="l", auto.key=TRUE)
#-> No need to modify densities in all models as difference not that great

#compare GAMs only to verify:
gam_all2 <- gam(MDperKMsqFall_mean ~ s(year, by=macrounit, bs="cs") + macrounit, data=AllMeans)
summary(gam_all2)
gam_all2_newmean <- gam(MDperKMsqFall_newmean ~ s(year, by=macrounit, bs="cs") + macrounit, data=AllMeans)
summary(gam_all2_newmean)
par(oma=c(2,0,2,0))
gam.check(gam_all2)
gam.check(gam_all2_newmean)
title("Gam_all2_newmean fall residual check", outer=TRUE)
# same

#calculate geometric mean of r:
geom <- prod(Popdata$RepRateFall_mean)^(1/length(Popdata$RepRateFall_mean))
mean(Popdata$RepRateFall_mean)

# plot PopDen Reprate(t-1) and (Fawn:Female Ratio)



# effect of WinterTemp on Fawn Female Ratio same as in Simones model?

gam_temp <- gam(FawnFemaleRatio_mean ~ s(AvrgWinterMinTemp, bs="cs"), data=AllMeans)
gam_temp <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +s(AvrgWinterMinTemp, bs="cs"), data=AllMeans)
summary(gam_temp)
plot(gam_temp)
AIC(gam_temp) #YES
gam_temppred<-predict(gam_temp, se.fit=T, newdata=data.frame(year=AllMeans$year, AvrgWinterMinTemp=AllMeans$AvrgWinterMinTemp), type="response")
plot(AllMeans$FawnFemaleRatio_mean~AllMeans$year, cex=0.5, main="GAM Fawn Female Ratio ~ Year and Winter Temp")
lines(gam_temppred$fit[1:51]~AllMeans$year[1:51],type="l", col="red", )
#NO

# Winter Mortality: Effects of Hunt and Temp_tplus1 (tplus1 because data refers to the previous winter)
par(mfrow=c(1,2))
plot(AllMeans$RatioSprAut_mean~AllMeans$year) # 2 outliers but removing them doesnt change a lot (see "leftover code.R")

plot(AllMeans$WinterMort_mean~AllMeans$year) # 2 outliers but removing them doesnt change a lot (see "leftover code.R")
plot(AllMeans$WinterMort_mean~AllMeans$HuntDen_All_mean, log="y")
plot(AllMeans$WinterMort_mean~AllMeans$AvrgWinterMinTemp, log="y")
# both dont show a trend -> probably no significant effect

gam_temp <- gam(WinterMort_mean ~ s(AvrgWinterMinTemp_tplus1,  bs="cs"), data=AllMeans)#NO EFFECT
gam_hunt <- gam(WinterMort_mean ~ s(HuntDen_All_mean,  bs="cs"), data=AllMeans)#NO EFFECT
gam_hunttemp <- gam(WinterMort_mean ~ te(AvrgWinterMinTemp_tplus1, HuntDen_All_mean,  bs="cs"), data=AllMeans)
plot(gam_hunt)
summary(gam_temp)
summary(gam_hunt)
summary(gam_hunttemp)
plot(gam_temp)
AIC(gam_temp)

# glm?
hist(AllMeans$WinterMort_mean)
glm_hunttemp <- glm(WinterMort_mean ~ AvrgWinterMinTemp_tplus1*HuntDen_All_mean, family=gaussian,data=AllMeans)#NO EFFECT

#seperate macrounits
glm_hunttemp <- glm(WinterMort_mean ~ AvrgWinterMinTemp_tplus1*HuntDen_All_mean, family=gaussian,data=AllMeans[which(AllMeans$macrounit == "0-2"),])#NO EFFECT
summary(glm_hunttemp)#NO EFFECT


# Relationship Spring+Fall Data
cor(AllMeans$MDperKMsqSpring_mean, AllMeans$MDperKMsqFall_mean, use="pairwise.complete.obs")#0.71
cor(Spring_tplus1, AllMeans$MDperKMsqFall_mean, use="pairwise.complete.obs")#0.64 -> does not make sense as less effects during winter (?)

plot(AllMeans$MDperKMsqSpring_mean_tplus1~ AllMeans$MDperKMsqFall_mean, main="Spring Population vs. Previous Fall Population", cex=0.7, ylab="Population Density Spring(t+1)", xlab="Population Density Fall(t)")
lines(1*unique(AllMeans$MDperKMsqFall_mean)~ unique(AllMeans$MDperKMsqFall_mean), type="l", col="red")
legend("topleft", legend=c("PopDen Spring(t+1)", "PopDen Fall(t)"), col=c("black", "red"), lty=1,2)

# at least as many datapoints above line than below --> bias in spring observations as less leaf cover?
# logical would be: lower values of spring population due to winter mortality+hunt as only factors betweens spring (t+1) and fall 

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

# GAMS using gamma-distribution

hist(AllMeans$MDperKMsqFall_mean)
gam_all2g <- gam(MDperKMsqFall_mean ~ s(year, by=macrounit, bs="cs") + macrounit, data=AllMeans, family=Gamma)
gam_all2gpred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit)
gam_all2gpred <- cbind(gam_all2gpred, predict(gam_all2g, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit), type="response"))
macrounitplots(glmobject = gam_all2gpred,title="GAM2g fall - interaction",colour="red")
summary(gam_all2g)
AIC(gam_all2g)
par(oma=c(2,0,2,0))
gam.check(gam_all2g)
title("gam_all2g fall residual check", outer=TRUE)
