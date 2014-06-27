library(mgcv)
library(reshape2)
source("Supplementary_Functions.R")


allmules <- read.csv("muledeer_final_dataset.csv")

### ALLMEANS: means per macrounit
# Only use the time series starting in 1962 and drop 3 study sites that were scarcely-surveyed:
mules1962 <- subset(allmules, year>1961)
del <- which(mules1962$StudyArea == "Bible_Camp" | mules1962$StudyArea == "NUTRNP"| mules1962$StudyArea == "SUTRNP")
mules1962 <- mules1962[-del,]

#Data Extraction and Modification
AllMeans <- extract(data=mules1962, fun=mean, xvar=c("MDperKMsqSpring", "MDperKMsqFall", "Average_of_Minimum_temperature_11_4", "d3","fall_density_coyote_by_macrounit_100km2", "WT_DEER_springsurveysD", "OIL_GAS_insideD", "woody_coverage", "MaleFall", "FemaleFall", "FawnFall"), listvar=c("macrounit","year"))
names(AllMeans) <- c("year", "macrounit","MDperKMsqSpring_mean", "MDperKMsqFall_mean",  "AvrgWinterMinTemp", "HuntDen_All_mean", "CoyoteDen_mean", "WTailDen_mean", "WellDen_mean", "WoodyVeg_mean", "MaleFall_mean", "FemaleFall_mean", "FawnFall_mean")

#add t+1/t+2 timelags to AllMeans to find out wether explanatory variables have a lagged effect
AllMeans <- cbind(AllMeans, "MDperKMsqFall_mean_t+1"= c(AllMeans$MDperKMsqFall_mean[-1],NA), "MDperKMsqFall_mean_t+2"= c(AllMeans$MDperKMsqFall_mean[-c(1,2)],NA,NA))

#calculate reproduction rate
AllMeans <- cbind(AllMeans, "RepRateFall_mean" = (AllMeans$FawnFall_mean/(AllMeans$FemaleFall_mean+AllMeans$MaleFall_mean)))


#recalculate densities to sq km²
AllMeans$CoyoteDen_mean <- AllMeans$CoyoteDen_mean / 100 
AllMeans$WTailDen_mean <- AllMeans$WTailDen_mean / 10 
AllMeans$WellDen_mean <- AllMeans$WellDen_mean / 10 


###WHOLEAREAMEANS: Means for whole area
WholeAreaMeans <- extract(data=mules1962, fun=mean, xvar=c("MDperKMsqSpring", "MDperKMsqFall", "Average_of_Minimum_temperature_11_4", "d3","fall_density_coyote_by_macrounit_100km2", "WT_DEER_springsurveysD", "OIL_GAS_insideD", "woody_coverage", "MaleFall", "FemaleFall", "FawnFall"), listvar=c("year"))
WholeAreaMeans <- WholeAreaMeans[,-1]
names(WholeAreaMeans) <- c("year","MDperKMsqSpring_mean", "MDperKMsqFall_mean",  "AvrgWinterMinTemp", "HuntDen_All_mean", "CoyoteDen_mean", "WTailDen_mean", "WellDen_mean", "WoodyVeg_mean", "MaleFall_mean", "FemaleFall_mean", "FawnFall_mean")


# for Gita
write.csv(AllMeans, "AllMeans.csv")





