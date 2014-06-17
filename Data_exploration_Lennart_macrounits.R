# Load extra packages:
library(lattice) # for xyplot
library(ggplot2)

# load supplementary functions
source("Supplementary_Functions.R")

# Load whole data set:
allmules <- read.csv("muledeer_final_dataset.csv")

# Gaps in Data
NAs_fall <- melt(tapply(X=allmules$TotalMDFall, INDEX=list(allmules$year), FUN=count_nas))
plot(NAs_fall, type="h", main="NAs in dataset")
NAs_spring <- melt(tapply(X=allmules$NumbMDSpring, INDEX=list(allmules$year), FUN=count_nas))
#plot(NAs_spring, type="h")
lines(NAs_spring, col="red")
legend("topright", legend=c("TotalMDFall","NumbMDSpring"), col=c("black","red"),lty=1)
#-> more gaps in fall, best to use data from 1962-2012
#-> but fall data includes sex+fawn data, more useful for future models

# Only use the time series starting in 1962 and drop 3 study sites that were scarcely-surveyed(see report:
mules1962 <- subset(allmules, year>1961)
del <- which(mules1962$StudyArea == "Bible_Camp" | mules1962$StudyArea == "NUTRNP"| mules1962$StudyArea == "SUTRNP")
mules1962 <- mules1962[-del,]

### --- Data Extraction
#raw data per macrounit
AllData <- extract(data=mules1962, columns=c("year", "macrounit","MDperKMsqSpring","d3","spring_density_coyote_by_macrounit_100km2", "WT_DEER_springsurveysD", "OIL_GAS_insideD", "woody_coverage"))
names(AllData) <- c("year", "macrounit","MDperKMsqSpring", "HuntDen_All", "CoyoteDen", "WTailDen", "WellDen", "WoodyVeg")
#means per macrounit
AllMeans <- extract(data=mules1962, fun=mean, xvar=c("MDperKMsqSpring", "MDperKMsqFall","d3","spring_density_coyote_by_macrounit_100km2", "WT_DEER_springsurveysD", "OIL_GAS_insideD", "woody_coverage"), listvar=c("macrounit","year"))
names(AllMeans) <- c("year", "macrounit","MDperKMsqSpring_mean", "MDperKMsqFall_mean", "HuntDen_All_mean", "CoyoteDen_mean", "WTailDen_mean", "WellDen_mean", "WoodyVeg_mean")
#means of whole area 
WholeAreaMeans <- extract(data=mules1962, fun=mean, xvar=c("MDperKMsqSpring","MDperKMsqFall", "d3","spring_density_coyote_by_macrounit_100km2", "WT_DEER_springsurveysD", "OIL_GAS_insideD", "woody_coverage"), listvar=c("year"))
WholeAreaMeans <- WholeAreaMeans[,-1]
names(WholeAreaMeans) <- c("year","MDperKMsqSpring_mean", "MDperKMsqFall_mean", "HuntDen_All_mean", "CoyoteDen_mean", "WTailDen_mean", "WellDen_mean", "WoodyVeg_mean")


### ---Comparison of Spring/Fall MD Densities
par(mfrow=c(1,2))
plot(AllMeans$MDperKMsqSpring_mean~AllMeans$year, type="p",pch=16, main="Spring")#comparison datapoints vs. overall mean
lines(WholeAreaMeans$MDperKMsqSpring~WholeAreaMeans$year, col="red")
plot(AllMeans$MDperKMsqFall_mean~AllMeans$year, type="p",pch=16, main="Fall")#comparison datapoints vs. overall mean
lines(WholeAreaMeans$MDperKMsqFall~WholeAreaMeans$year,col="red")
# show same pattern


### ---Plots of all factors per macrounit

# Spring population density (sq km) (raw data per study area)

PopDenSpring <- data.frame(StudyArea = mules1962$StudyArea, macrounit = mules1962$macrounit, year = mules1962$year, MDperKMsqSpring = mules1962$MDperKMsqSpring)
PopDenSpring_mean <- t(tapply(X=PopDenSpring$MDperKMsqSpring, INDEX=c(list(PopDenSpring$macrounit), list(PopDenSpring$year)), FUN=mean, na.rm=T))
PopDenSpring_mean_m <- melt(PopDenSpring_mean)
names(PopDenSpring_mean_m) <- c("year", "macrounit", "PopDenSpring_mean")

# Harvest density antlered+antlerless (sq km) (raw data per macrounit, density per size of badlands in each huntingunit)
HuntDenAll  <- data.frame(StudyArea = mules1962$StudyArea, macrounit = mules1962$macrounit, year = mules1962$year, HuntDenAll = mules1962$d3) 
HuntDenAll_mean <- t(tapply(X=HuntDenAll$HuntDenAll, INDEX=c(list(HuntDenAll$macrounit), list(HuntDenAll$year)), FUN = mean, na.rm=T)) #not really necessary as hunting data only available per huntin macrounit (but easiest way to extract data)
HuntDenAll_mean_m <- melt(HuntDenAll_mean)
names(HuntDenAll_mean_m) <- c("year", "macrounit", "HuntDen_mean")

#--------
# Coyote Density (raw data per study area, regrouped into macrounits)(100 sq km)
CoyoteDen <- data.frame(StudyArea = mules1962$StudyArea, macrounit = mules1962$macrounit, year = mules1962$year, CoyoteDen = mules1962$spring_density_coyote_by_macrounit_100km2)
CoyoteDen_mean <- t(tapply(X=CoyoteDen$CoyoteDen, INDEX=c(list(CoyoteDen$macrounit), list(CoyoteDen$year)), FUN = mean, na.rm=T))
CoyoteDen_mean_m <- melt(CoyoteDen_mean)
names(CoyoteDen_mean_m) <- c("year", "macrounit", "CoyoteDen_mean")
CoyoteDen_mean_m$CoyoteDen_mean <- CoyoteDen_mean_m$CoyoteDen_mean / 100 #per sq km


# White Tail Deer Spring Density (10 sq km) (raw data per study site)
WTailDen <- data.frame(StudyArea = mules1962$StudyArea, macrounit = mules1962$macrounit, year = mules1962$year, WTailDen = mules1962$WT_DEER_springsurveysD)
WTailDen_mean <- t(tapply(X=WTailDen$WTailDen, INDEX=c(list(WTailDen$macrounit), list(WTailDen$year)), FUN = mean, na.rm=T))
WTailDen_mean_m <- melt(WTailDen_mean)
names(WTailDen_mean_m) <- c("year", "macrounit", "WTailDen_mean")
WTailDen_mean_m$WTailDen_mean <- WTailDen_mean_m$WTailDen_mean / 10 # per sq km

#-----

# Density of active oil wells within study site including 1000m buffer (10sq km) (only group "oil+gas")
WellDen <- data.frame(StudyArea = mules1962$StudyArea, macrounit = mules1962$macrounit, year = mules1962$year, WellDen = mules1962$OIL_GAS_insideD)
WellDen_mean <- t(tapply(X=WellDen$WellDen, INDEX=c(list(WellDen$macrounit), list(WellDen$year)), FUN=mean, na.rm=T))
WellDen_mean_m <- melt(WellDen_mean)
names(WellDen_mean_m) <- c("year", "macrounit", "WellDen_mean")
WellDen_mean_m$WellDen_mean <- WellDen_mean_m$WellDen_mean / 10 # in sq km

# Woody Vegetation from model (percentage per study site)
WoodyVeg <- data.frame(StudyArea = mules1962$StudyArea, macrounit = mules1962$macrounit, year = mules1962$year, WoodyVeg = mules1962$woody_coverage)
WoodyVeg_mean <- t(tapply(X=WoodyVeg$WoodyVeg, INDEX=c(list(WoodyVeg$macrounit), list(WoodyVeg$year)), FUN = mean, na.rm=T))
WoodyVeg_mean_m <- melt(WoodyVeg_mean)
names(WoodyVeg_mean_m) <- c("year", "macrounit", "WoodyVeg_mean")
WoodyVeg_mean_m <- (WoodyVeg_mean_m)

MyPanel <- function(x, y, subscripts) {
  panel.xyplot(x, y, type = "l")
  panel.abline(lm(y ~ x))}

xyplot(PopDenSpring_mean + HuntDenAll_mean ~ year | macrounit, data=c(PopDenSpring_mean_m, HuntDenAll_mean_m), subscripts=TRUE, panel=function(...) panel.superpose(panel.groups=MyPanel, ...), type="l", auto.key = list(space = "top", text = c("Population Density", "Hunting Density")), xlab = "Year", ylab= "Density per km?", main = "Mean Spring Population Density and Mean Harvesting Density")

# so will ich http://stackoverflow.com/questions/16361766/multiple-ablines-in-xyplot

xyplot(CoyoteDen_mean + WTailDen_mean ~ year | macrounit, data=c(CoyoteDen_mean_m, WTailDen_mean_m), type="l", auto.key = list(space = "top", text = c("Coyote Density", "WT Deer Density"), points = FALSE, lines = TRUE), xlab = "Year", ylab= "Density per km?", main = "Mean Coyote and White-Tail Deer Spring Population Density")

xyplot(WellDen_mean + (WoodyVeg_mean/100) ~ year | macrounit, data=c(WellDen_mean_m, WoodyVeg_mean_m), type="l", auto.key = list(space = "top", text = c("Density Oil+Gas Sites per km? (incl. 1km-buffer)", "Percentage of Woody Vegetation (divided by 100)"), points = FALSE, lines = TRUE), xlab = "Year", ylab= "Density/Percentage", main = "Mean Oil+Gas Extraction Density and Percentage of Woody Vegetation")

#------------------------------------------
#lm of popdensity  --> of no use....
sub_PopDenSpring01 <- subset(PopDenSpring_mean_m, macrounit == "0-1")
sub_PopDenSpring02 <- subset(PopDenSpring_mean_m, macrounit == "0-2")
sub_PopDenSpring03<- subset(PopDenSpring_mean_m, macrounit == "0-3")
sub_PopDenSpring04 <- subset(PopDenSpring_mean_m, macrounit == "0-4")

lm_sub_PopDenSpring01 <- lm(PopDenSpring_mean~year, data=sub_PopDenSpring01)
lm_sub_PopDenSpring02 <- lm(PopDenSpring_mean~year, data=sub_PopDenSpring02)
lm_sub_PopDenSpring03 <- lm(PopDenSpring_mean~year, data=sub_PopDenSpring03)
lm_sub_PopDenSpring04 <- lm(PopDenSpring_mean~year, data=sub_PopDenSpring04)

par(mfrow=c(2,2))
plot(sub_PopDenSpring01$PopDenSpring_mean~sub_PopDenSpring01$year, main="0 - 1", ylab="per km?", xlab="Year", cex=0.5, col="blue")
abline(lm_sub_PopDenSpring01, add=TRUE, col=4)
plot(sub_PopDenSpring02$PopDenSpring_mean~sub_PopDenSpring02$year, main="0 - 2", ylab="per km?", xlab="Year", cex=0.5, col="blue")
abline(lm_sub_PopDenSpring01, add=TRUE, col=4)
plot(sub_PopDenSpring03$PopDenSpring_mean~sub_PopDenSpring03$year, main="0 - 3", ylab="per km?", xlab="Year", cex=0.5, col="blue")
abline(lm_sub_PopDenSpring01, add=TRUE, col=4)
plot(sub_PopDenSpring04$PopDenSpring_mean~sub_PopDenSpring04$year, main="0 - 4", ylab="per km?", xlab="Year", cex=0.5, col="blue")
abline(lm_sub_PopDenSpring01, add=TRUE, col=4)




