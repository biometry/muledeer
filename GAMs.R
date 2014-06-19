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
gam_all0pred <- cbind(gam_all0pred, predict(gam_all0, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit), type="response"))
macrounitplots(gam_all0pred, col=2)
# One smoother for all macrounits
gam_all1 <- gam(MDperKMsqFall_mean ~ s(year) + factor(macrounit), data=AllMeans)
gam_all1pred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit)
gam_all1pred <- cbind(gam_all1pred, predict(gam_all1, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit), type="response"))
# Seperate smoothers for each macrounit
gam_all2 <- gam(MDperKMsqFall_mean ~ s(year, by=macrounit) + macrounit, data=AllMeans)
gam_all2 <- gam(MDperKMsqFall_mean ~ s(year, by=macrounit) + macrounit,  data=AllMeans)


gam_all2pred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit)
gam_all2pred <- cbind(gam_all2pred, predict(gam_all2, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit), type="response"))
macrounitplots(glmobject = gam_all2pred,title="GAM2 fall - interaction",colour="red")
summary(gam_all2)
par(oma=c(2,0,2,0))
gam.check(gam_all2)#residual show clear patterns -> normal distribution?variance homogenity?
title("Gam_all2 fall residual check", outer=TRUE)
#comparison nullmodell gam_all09 and gam_all2
anova(gam_all0, gam_all2, test="F")

# Check Residuals for temporal autocorrelation
gam_all2res <- residuals(gam_all2, type = "deviance")
plot(gam_all2res ~AllMeans$year[which(!is.na(AllMeans$MDperKMsqFall_mean))]) #pattern not too bad?
acf(gam_all2res, na.action = na.pass,main = "Auto-correlation plot for residuals Gam_all2 fall")


###-------- GAM for whole area including autocorrelation, macrounit as an interaction factor

# Seperate smoothers for each macrounit
gam_all2c <- gamm(MDperKMsqFall_mean ~ s(year, by=macrounit) + macrounit, data=AllMeans, correlation=corAR1(form= ~ year | macrounit))
gam_all2c <- gamm(MDperKMsqFall_mean ~ s(year, by=macrounit) + macrounit, control=list(niterEM=0,optimMethod="GCV.Cp"), data=AllMeans, correlation=corAR1(form= ~ year| macrounit))

# why is smooth so different?
#GCP.Cp is default method of gam AND gamm

gam_all2cpred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit)
gam_all2cpred <- cbind(gam_all2cpred, predict(gam_all2c$gam, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit), type="response"))
macrounitplots(glmobject = gam_all2cpred,title="Gam_all2c fall - interaction",colour="red")
summary(gam_all2c$gam)
par(oma=c(2,0,2,0))
gam.check(gam_all2c$gam)
title("gam_all2c fall residual check", outer=TRUE)
AIC(gam_all2c$lme)
#comparison nullmodell gam_all09 and gam_all2c
anova(gam_all0, gam_all2c, test="F")

gamm:
control=list(niterEM=0,optimMethod="GCV.Cp")
gamm(formula,random=NULL,correlation=NULL,family=gaussian(),
     data=list(),weights=NULL,subset=NULL,na.action,knots=NULL,
     control=list(niterEM=0,optimMethod="L-BFGS-B"),
     niterPQL=20,verbosePQL=TRUE,method="ML",...)

control=list(niterEM=0,optimMethod="GCV.Cp")






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
macrounitplots(glmobject = glm_all2pred, title="GLM2 - macrounit as predictor",colour="blue")
gam.check(glm_all2)
#-> does not fit the data well,residuals look a lot worse than the ones of gam:
#- especially assumption of variance homogenity (Resids vs.linear pred plot)
#- bias in response vs. fitted values plot, quadratic shape insted of linear
#- bias in histogram towards the right, not symmetrical

###-------- GAM for whole area, no autocorrelation, macrounit as an interaction factor

# One smoother without effect of macrounits
gam_all0 <- gam(MDperKMsqSpring_mean ~ s(year), data=AllMeans)
gam_all0pred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit)
gam_all0pred <- cbind(gam_all0pred, predict(gam_all0, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit), type="response"))
# One smoother for all macrounits
gam_all1 <- gam(MDperKMsqSpring_mean ~ s(year) + factor(macrounit), data=AllMeans)
gam_all1pred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit)
gam_all1pred <- cbind(gam_all1pred, predict(gam_all1, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit), type="response"))
# Seperate smoothers for each macrounit
gam_all2 <- gam(MDperKMsqSpring_mean ~ s(year, by=macrounit) + macrounit, data=AllMeans)
gam_all2pred <- data.frame(year=AllMeans$year, macrounit=AllMeans$macrounit)
gam_all2pred <- cbind(gam_all2pred, predict(gam_all2, se.fit=T, newdata=data.frame("year"=AllMeans$year, "macrounit"=AllMeans$macrounit), type="response"))
macrounitplots(glmobject = gam_all2pred,title="GAM2 spring - interaction",colour="red")
par(oma=c(2,0,2,0))
gam.check(gam_all2)
title("Gam_all2 spring residual check", outer=TRUE)

