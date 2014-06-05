#Schritte:
- alles in xyplot mit xyplot(PopDens + HuntDens ~ year | macrounit)

 

setwd("C:/Users/lschmidt/Desktop/Mule_deer/Data_exploration")

# Load whole data set:
allmules <- read.csv("muledeer_final_dataset.csv")

# Only use the time series starting in 1962:
mules1962 <- subset(allmules, year>1961)

# Load extra packages:
library(lattice) # for xyplot
library(ggplot2)
theme_set(theme_bw()) # Black-and-white theme for ggplot instead of the default grey background



# ALL IN ONE PER MACROUNIT----------------------------------------------------------

# Spring population density (sq km) (raw data per study area)
PopDenSpring <- data.frame(StudyArea = mules1962$StudyArea, macrounit = mules1962$macrounit, year = mules1962$year, MDperKMsqSpring = mules1962$MDperKMsqSpring)
PopDenSpring_mean <- tapply(X=PopDenSpring$MDperKMsqSpring, INDEX=c(list(PopDenSpring$macrounit), list(PopDenSpring$year)), FUN=mean, na.rm=T)
PopDenSpring_mean <- as.data.frame(PopDenSpring_mean) #notwendig? -> das resultat ist irgendwie trotzdem eine "list"
xyplot(PopDenSpring_mean ~ names(PopDenSpring_mean) | row.names(PopDenSpring_mean), data=PopDenSpring_mean, type="l")


# Harvest density antlered+antlerless (sq km) (raw data per macrounit, density per size of badlands in each huntingunit)
HuntDenAll  <- data.frame(StudyArea = mules1962$StudyArea, macrounit = mules1962$macrounit, year = mules1962$year, HuntDenAll = mules1962$d3) 
HuntDenAll_mean <- tapply(X=HuntDenAll$HuntDenAll, INDEX=c(list(HuntDenAll$macrounit), list(HuntDenAll$year)), FUN = mean, na.rm=T) #not really necessary as hunting data only available per huntin macrounit (but easiest way to extract data)


# Density of active oil wells within study site including 1000m buffer (10sq km) (only group "oil+gas")
WellDen <- data.frame(StudyArea = mules1962$StudyArea, macrounit = mules1962$macrounit, year = mules1962$year, WellDen = mules1962$OIL_GAS_insideD)
WellDen_mean <- tapply(X=WellDen$WellDen, INDEX=c(list(WellDen$macrounit), list(WellDen$year)), FUN=mean, na.rm=T)
WellDen_mean <- WellDen_mean * 100 # in sq km


# Coyote Density (raw data per study area, regrouped into macrounits)(100 sq km)
CoyoteDen <- data.frame(StudyArea = mules1962$StudyArea, macrounit = mules1962$macrounit, year = mules1962$year, CoyoteDen = mules1962$fall_density_coyote_by_macrounit_100km2)
CoyoteDen_mean <- tapply(X=CoyoteDen$CoyoteDen, INDEX=c(list(CoyoteDen$macrounit), list(CoyoteDen$year)), FUN = mean, na.rm=T)
CoyoteDen_mean <- CoyoteDen_mean * 10000 #per sq km


# White Tail Deer Spring Density (10 sq km) (raw data per study site)
WTailDen <- data.frame(StudyArea = mules1962$StudyArea, macrounit = mules1962$macrounit, year = mules1962$year, WTailDen = mules1962$WT_DEER_springsurveysD)
WTailDen_mean <- tapply(X=WTailDen$WTailDen, INDEX=c(list(WTailDen$macrounit), list(WTailDen$year)), FUN = mean, na.rm=T)
WTailDen_mean <- WTailDen_mean * 100 # in sq km


# Woody Vegetation from model (percentage per study site)
WoodyVeg <- data.frame(StudyArea = mules1962$StudyArea, macrounit = mules1962$macrounit, year = mules1962$year, WoodyVeg = mules1962$woody_coverage)
WoodyVeg_mean <- tapply(X=WoodyVeg$WoodyVeg, INDEX=c(list(WoodyVeg$macrounit), list(WoodyVeg$year)), FUN = mean, na.rm=T)
