###------Generalized models for evaluation of effects of different input factors

library(mgcv)
library(reshape2)
source("Supplementary_Functions.R")
allmules <- read.csv("muledeer_final_dataset.csv")

# Only use the time series starting in 1962 and drop 3 study sites that were scarcely-surveyed:
mules1962 <- subset(allmules, year>1961)
del <- which(mules1962$StudyArea == "Bible_Camp" | mules1962$StudyArea == "NUTRNP"| mules1962$StudyArea == "SUTRNP")
mules1962 <- mules1962[-del,]

#Data Extraction
AllMeans <- extract(data=mules1962, fun=mean, xvar=c("MDperKMsqSpring", "MDperKMsqFall","d3","spring_density_coyote_by_macrounit_100km2", "WT_DEER_springsurveysD", "OIL_GAS_insideD", "woody_coverage"), listvar=c("macrounit","year"))
names(AllMeans) <- c("year", "macrounit","MDperKMsqSpring_mean", "MDperKMsqFall_mean", "HuntDen_All_mean", "CoyoteDen_mean", "WTailDen_mean", "WellDen_mean", "WoodyVeg_mean")
WholeAreaMeans <- extract(data=mules1962, fun=mean, xvar=c("MDperKMsqSpring","MDperKMsqFall", "d3","spring_density_coyote_by_macrounit_100km2", "WT_DEER_springsurveysD", "OIL_GAS_insideD", "woody_coverage"), listvar=c("year"))
WholeAreaMeans <- WholeAreaMeans[,-1]
names(WholeAreaMeans) <- c("year","MDperKMsqSpring_mean", "MDperKMsqFall_mean", "HuntDen_All_mean", "CoyoteDen_mean", "WTailDen_mean", "WellDen_mean", "WoodyVeg_mean")

######################SPRING MD DENSITIES##########################################

par(mfrow=c(1,2))
plot(AllMeans$MDperKMsqSpring_mean~AllMeans$year, type="p",pch=16, main="Spring")#comparison datapoints vs. overall mean
lines(WholeAreaMeans$MDperKMsqSpring~WholeAreaMeans$year, col="red")
plot(AllMeans$MDperKMsqFall_mean~AllMeans$year, type="p",pch=16, main="Fall")#comparison datapoints vs. overall mean
lines(WholeAreaMeans$MDperKMsqFall~WholeAreaMeans$year,col="red")
###------1st try: GLM

# linear regression for comparison
glm_all <- glm(MDperKMsqSpring_mean~year,family=gaussian,data=AllMeans)
glm_allpred <- predict(glm_all,se.fit=T, type="response")
plot(AllMeans$MDperKMsqSpring_mean~AllMeans$year, type="p",pch=16)#comparison datapoints vs. mean
lines(glm_allpred$fit~AllMeans$year, col="red")
# glm including macrounits as predictor
glm_all2 <- glm(MDperKMsqSpring_mean~year+macrounit,family=gaussian,data=AllMeans)
glm_all2pred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit)
glm_all2pred <- cbind(glm_all2pred,as.data.frame(predict(glm_all2, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit),se.fit=T, type="response")))
macrounitplots(glmobject = glm_all2pred, title="GLM2 quadratic term - macrounit as predictor",colour="blue")
gam.check(glm_all2)
#-> does not fit the data well,residuals look a lot worse than the ones of gam:
#- especially assumption of variance homogenity (Resids vs.linear pred plot)
#- bias in response vs. fitted values plot, quadratic shape insted of linear
#- bias in histogram towards the right, not symmetrical

###-------- GAM for whole area, no autocorrelation, macrounit as an interaction factor

# One smoother without effect of macrounits
gam_all0 <- gam(MDperKMsqSpring_mean ~ s(year), data=AllMeans)
gam_all0pred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit)
gam_all0pred <- cbind(gam_all0pred, predict(gam_all0, se.fit=T, type="response"))
# One smoother for all macrounits
gam_all1 <- gam(MDperKMsqSpring_mean ~ s(year) + factor(macrounit), data=AllMeans)
gam_all1pred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit)
gam_all1pred <- cbind(gam_all1pred, predict(gam_all1, se.fit=T, type="response"))
# Seperate smoothers for each macrounit
gam_all2 <- gam(MDperKMsqSpring_mean ~ s(year, by=as.numeric(macrounit == "0-1"))+s(year, by=as.numeric(macrounit == "0-2"))+s(year, by=as.numeric(macrounit == "0-3"))+s(year, by=as.numeric(macrounit == "0-4"))+factor(macrounit), data=AllMeans)
gam_all2pred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit)
gam_all2pred <- cbind(gam_all2pred, predict(gam_all2, se.fit=T, type="response"))
macrounitplots(glmobject = gam_all2pred,title="GAM2 - interaction",colour="red")


######################FALL MD DENSITIES##########################################

###----------GLMs

glm_all <- glm(MDperKMsqFall_mean~year,family=gaussian,data=AllMeans)
glm_allpred <- predict(glm_all,se.fit=T, type="response")
lines(glm_allpred$fit~AllMeans$year, col="red")
# glm including macrounits as predictor

glm_all2 <- glm(MDperKMsqFall_mean~year+macrounit,family=gaussian,data=AllMeans)
glm_all2pred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit)
glm_all2pred <- cbind(glm_all2pred,as.data.frame(predict(glm_all2, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit), se.fit=T, type="response")))
macrounitplots(glmobject = glm_all2pred, title="GLM2 quadratic term - macrounit as predictor",colour="blue")
summary(glm_all2)
gam.check(glm_all2)
#-> does not fit the data well,residuals look a lot worse than the ones of spring and all gams:
#- especially assumption of variance homogenity (Resids vs.linear pred plot)
#- bias in response vs. fitted values plot, quadratic shape insted of linear
#- bias in histogram towards the right, not symmetrical

###-------- GAM for whole area, no autocorrelation, macrounit as an interaction factor

plot(AllMeans$MDperKMsqFall_mean~AllMeans$year, type="p",pch=16)#comparison datapoints vs. mean
lines(WholeAreaMeans$MDperKMsqFall_mean~WholeAreaMeans$year)

# One smoother without effect of macrounits
gam_all0 <- gam(MDperKMsqFall_mean ~ s(year), data=AllMeans)
gam_all0pred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit)
gam_all0pred <- cbind(gam_all0pred, predict(gam_all0, se.fit=T, type="response"))
# One smoother for all macrounits
gam_all1 <- gam(MDperKMsqFall_mean ~ s(year) + factor(macrounit), data=AllMeans)
gam_all1pred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit)
gam_all1pred <- cbind(gam_all1pred, predict(gam_all1, se.fit=T, type="response"))
# Seperate smoothers for each macrounit
gam_all2 <- gam(MDperKMsqFall_mean ~ s(year, by=as.numeric(macrounit == "0-1"))+s(year, by=as.numeric(macrounit == "0-2"))+s(year, by=as.numeric(macrounit == "0-3"))+s(year, by=as.numeric(macrounit == "0-4"))+factor(macrounit), data=AllMeans)
gam_all2pred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit)
gam_all2pred <- cbind(gam_all2pred, predict(gam_all2, se.fit=T, type="response"))
macrounitplots(glmobject = gam_all2pred,title="GAM2 - interaction",colour="red")