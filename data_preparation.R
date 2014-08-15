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
AllMeans <- extract(data=mules1962, fun=mean, xvar=c("MDperKMsqSpring", "MDperKMsqFall", "Average_of_Minimum_temperature_11_4", "d3","d2", "fall_density_coyote_by_macrounit_100km2", "WT_DEER_springsurveysD", "OIL_GAS_insideD", "woody_coverage", "MaleFall", "FemaleFall", "FawnFall"), listvar=c("macrounit","year"))
names(AllMeans) <- c("year", "macrounit","MDperKMsqSpring_mean", "MDperKMsqFall_mean",  "AvrgWinterMinTemp", "HuntDen_All_mean", "HuntDen_Aless_mean", "CoyoteDen_mean", "WTailDen_mean", "WellDen_mean", "WoodyVeg_mean", "MaleFall_mean", "FemaleFall_mean", "FawnFall_mean")

#add t+1/t+2 timelags to AllMeans to find out wether explanatory variables have a lagged effect
AllMeans <- cbind(AllMeans, "MDperKMsqFall_mean_tplus1"= c(AllMeans$MDperKMsqFall_mean[-1],NA), "MDperKMsqFall_mean_tplus2"= c(AllMeans$MDperKMsqFall_mean[-c(1,2)],NA,NA))




#add timelags to Spring Population Density
MDperKMsqSpring_mean_tminus1 <- c(NA, AllMeans$MDperKMsqSpring_mean[-length(AllMeans$MDperKMsqSpring_mean)])
AllMeans <- cbind(AllMeans,MDperKMsqSpring_mean_tminus1)

AllMeans <- cbind(AllMeans, "MDperKMsqSpring_mean_tplus1"= c(AllMeans$MDperKMsqSpring_mean[-1],NA))


# timelags hunting
AllMeans <- cbind(AllMeans,HuntDen_All_mean_tminus1=c(NA, AllMeans$HuntDen_All_mean[-length(AllMeans$HuntDen_All_mean)]))
AllMeans <- cbind(AllMeans,HuntDen_Aless_mean_tminus1=c(NA, AllMeans$HuntDen_Aless_mean[-length(AllMeans$HuntDen_Aless_mean)]))

# timelag AvrgminTemp
AllMeans <- cbind(AllMeans,AvrgWinterMinTemp_tminus1=c(NA, AllMeans$AvrgWinterMinTemp[-length(AllMeans$AvrgWinterMinTemp)]))

# Calculate Winter mortality
WinterMort_mean <- AllMeans$MDperKMsqFall_mean - AllMeans$MDperKMsqSpring_mean_tplus1 
AllMeans <- cbind(AllMeans,WinterMort_mean)


#calculate fawn:total number ratio
AllMeans <- cbind(AllMeans, "FawnTotalRatio_mean" = (AllMeans$FawnFall_mean/(AllMeans$FemaleFall_mean+AllMeans$MaleFall_mean)))
AllMeans <- cbind(AllMeans, "FawnTotalRatio_mean_tminus1" = c(NA, AllMeans$FawnTotalRatio_mean[-length(AllMeans$FawnTotalRatio_mean)]))

# calculate Fawn:Female ratio
FawnFemaleRatio_mean = (AllMeans$FawnFall_mean/AllMeans$FemaleFall_mean)
AllMeans <- cbind(AllMeans, FawnFemaleRatio_mean)
#FawnFemaleRatio time lag 
AllMeans <- cbind(AllMeans,FawnFemaleRatio_mean_tminus1 = c(NA, AllMeans$FawnFemaleRatio_mean[-length(AllMeans$FawnFemaleRatio_mean)]))
# calculate reproduction rate per individual
AllMeans <- cbind(AllMeans, "RepRateFall_mean" = ((AllMeans$MDperKMsqFall_mean_tplus1) - (AllMeans$MDperKMsqFall_mean))/AllMeans$MDperKMsqFall_mean)
AllMeans <- cbind(AllMeans, "RepRateFall_mean_tminus1" = c(NA, AllMeans$RepRateFall_mean[-length(AllMeans$RepRateFall_mean)]))
                                    


#recalculate densities to sq km²
AllMeans$CoyoteDen_mean <- AllMeans$CoyoteDen_mean / 100 
AllMeans$WTailDen_mean <- AllMeans$WTailDen_mean / 10 
AllMeans$WellDen_mean <- AllMeans$WellDen_mean / 10 


###WHOLEAREAMEANS: Means for whole area
WholeAreaMeans <- extract(data=mules1962, fun=mean, xvar=c("MDperKMsqSpring", "MDperKMsqFall", "Average_of_Minimum_temperature_11_4", "d3", "d2","fall_density_coyote_by_macrounit_100km2", "WT_DEER_springsurveysD", "OIL_GAS_insideD", "woody_coverage", "MaleFall", "FemaleFall", "FawnFall"), listvar=c("year"))
WholeAreaMeans <- WholeAreaMeans[,-1]
names(WholeAreaMeans) <- c("year","MDperKMsqSpring_mean", "MDperKMsqFall_mean",  "AvrgWinterMinTemp", "HuntDen_All_mean", "HuntDen_Aless_mean", "CoyoteDen_mean", "WTailDen_mean", "WellDen_mean", "WoodyVeg_mean", "MaleFall_mean", "FemaleFall_mean", "FawnFall_mean")


# for Gita
write.csv(AllMeans, "AllMeans.csv")
write.csv(WholeAreaMeans, "WholeAreaMeans.csv")

remove(allmules)

remove(del)
remove(FawnFemaleRatio_mean)
remove(WinterMort_mean)

