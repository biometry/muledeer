setwd("/home/gita/Dokumente/Projekte/Mule_deer_Simone/Data_exploration/")

# Load whole data set:
allmules <- read.csv("muledeer_final_dataset.csv")

# Only use the time series starting in 1962:
mules1962 <- subset(allmules, year>1961)

# How long is the time series?
length(levels(as.factor(allmules$year))) # 57 years
length(levels(as.factor(mules1962$year))) # 51 

# How many sites are included?
length(levels(allmules$StudyArea)) # 26


# Load extra packages:
library(lattice) # for xyplot
library(ggplot2)
theme_set(theme_bw()) # Black-and-white theme for ggplot instead of the default grey background

# 1. Mule deer population ----------------------------------------------------

# Identify gaps in the mule deer time series:

count_nas <- function(x) length(which(is.na(x)))

# NAs in full data set:
NAs_fall <- tapply(X=allmules$TotalMDFall, INDEX=list(allmules$year), FUN=count_nas)
summary(NAs_fall)
NAs_spring <- tapply(X=allmules$NumbMDSpring, INDEX=list(allmules$year), FUN=count_nas)
summary(NAs_spring)

# NAs in reduced data set:
NAs_fall <- tapply(X=mules1962$TotalMDFall, INDEX=list(mules1962$year), FUN=count_nas)
sum(NAs_fall)
NAs_spring <- tapply(X=mules1962$NumbMDSpring, INDEX=list(mules1962$year), FUN=count_nas)
sum(NAs_spring)


# Whole population time series --------------------------------------------

# Mean population size over time:
meanpop_fall <- tapply(X=mules1962$TotalMDFall, INDEX=list(mules1962$year), FUN=mean, na.rm=T)
maxpop_fall <- tapply(X=mules1962$TotalMDFall, INDEX=list(mules1962$year), FUN=max, na.rm=T)
minpop_fall <- tapply(X=mules1962$TotalMDFall, INDEX=list(mules1962$year), FUN=min, na.rm=T)

png("Population_counts_mean_range_fall.png", width=600, height=600)
par(cex=1.3)
plot(meanpop_fall~names(meanpop_fall), type="b", ylim=c(0, max(maxpop_fall)), xlab="Year", ylab="Fall population count")
segments(x0=as.numeric(names(meanpop_fall)), y0=minpop_fall, y1=maxpop_fall)
dev.off()

meanpop_spring <- tapply(X=mules1962$NumbMDSpring, INDEX=list(mules1962$year), FUN=mean, na.rm=T)
maxpop_spring <- tapply(X=mules1962$NumbMDSpring, INDEX=list(mules1962$year), FUN=max, na.rm=T)
minpop_spring <- tapply(X=mules1962$NumbMDSpring, INDEX=list(mules1962$year), FUN=min, na.rm=T)

png("Population_counts_mean_range_spring.png", width=600, height=600)
par(cex=1.3)
plot(meanpop_spring~names(meanpop_spring), type="b", ylim=c(0, max(maxpop_spring)), xlab="Year", ylab="Spring population count")
segments(x0=as.numeric(names(meanpop_spring)), y0=minpop_spring, y1=maxpop_spring)
dev.off()

# Separate plots for each site:
png("Population_counts_single_sites_spring.png", width=1100, height=650)
xyplot(NumbMDSpring ~ year | StudyArea, data=mules1962, type="l")
dev.off()

png("Population_counts_single_sites_fall.png", width=1100, height=650)
xyplot(TotalMDFall ~ year | StudyArea, data=mules1962, type="l")
dev.off()


# Same plots with density/km² instead of counts ---------------------------

# Mean population size over time:
meandens_fall <- tapply(X=mules1962$MDperKMsqFall, INDEX=list(mules1962$year), FUN=mean, na.rm=T)
maxdens_fall <- tapply(X=mules1962$MDperKMsqFall, INDEX=list(mules1962$year), FUN=max, na.rm=T)
mindens_fall <- tapply(X=mules1962$MDperKMsqFall, INDEX=list(mules1962$year), FUN=min, na.rm=T)

png("Population_density_mean_range_fall.png", width=600, height=600)
par(cex=1.3)
plot(meandens_fall~names(meandens_fall), type="b", ylim=c(0, max(maxdens_fall)), xlab="Year", ylab="Fall population density / km²")
segments(x0=as.numeric(names(meandens_fall)), y0=mindens_fall, y1=maxdens_fall)
dev.off()

meandens_spring <- tapply(X=mules1962$MDperKMsqSpring, INDEX=list(mules1962$year), FUN=mean, na.rm=T)
maxdens_spring <- tapply(X=mules1962$MDperKMsqSpring, INDEX=list(mules1962$year), FUN=max, na.rm=T)
mindens_spring <- tapply(X=mules1962$MDperKMsqSpring, INDEX=list(mules1962$year), FUN=min, na.rm=T)

png("Population_density_mean_range_spring.png", width=600, height=600)
par(cex=1.3)
plot(meandens_spring~names(meandens_spring), type="b", ylim=c(0, max(maxdens_spring)), xlab="Year", ylab="Spring population density / km²")
segments(x0=as.numeric(names(meandens_spring)), y0=mindens_spring, y1=maxdens_spring)
dev.off()

# Separate plots for each site:
png("Population_density_single_sites_spring.png", width=1100, height=650)
xyplot(MDperKMsqSpring ~ year | StudyArea, data=mules1962, type="l")
dev.off()

png("Population_density_single_sites_fall.png", width=1100, height=650)
xyplot(MDperKMsqFall ~ year | StudyArea, data=mules1962, type="l")
dev.off()


# Fawn:female ratio -------------------------------------------------------
# (Fall data)

ffratio_mean <- tapply(X=mules1962$Fw.FratioFall, INDEX=list(mules1962$year), FUN=mean, na.rm=T)
ffratio_max <- tapply(X=mules1962$Fw.FratioFall, INDEX=list(mules1962$year), FUN=max, na.rm=T)
ffratio_min <- tapply(X=mules1962$Fw.FratioFall, INDEX=list(mules1962$year), FUN=min, na.rm=T)

png("Fawn_female_ratio_mean_range.png", width=600, height=600)
par(cex=1.3)
plot(ffratio_mean~names(ffratio_mean), type="b", ylim=c(0, max(ffratio_max)), xlab="Year", ylab="Fawn:female ratio")
segments(x0=as.numeric(names(ffratio_mean)), y0=ffratio_min, y1=ffratio_max)
dev.off()

png("Fawn_female_ratio_single_sites.png", width=1100, height=650)
xyplot(Fw.FratioFall ~ year | StudyArea, data=mules1962, type="l")
dev.off()


# Male:female ratio -------------------------------------------------------

mfratio_mean <- tapply(X=mules1962$M.FratioFall, INDEX=list(mules1962$year), FUN=mean, na.rm=T)
mfratio_max <- tapply(X=mules1962$M.FratioFall, INDEX=list(mules1962$year), FUN=max, na.rm=T)
mfratio_min <- tapply(X=mules1962$M.FratioFall, INDEX=list(mules1962$year), FUN=min, na.rm=T)

png("Male_female_ratio_mean_range.png", width=600, height=600)
par(cex=1.3)
plot(mfratio_mean~names(mfratio_mean), type="b", ylim=c(0, max(mfratio_max)), xlab="Year", ylab="Male:female ratio")
segments(x0=as.numeric(names(mfratio_mean)), y0=mfratio_min, y1=mfratio_max)
dev.off()

png("Male_female_ratio_single_sites.png", width=1100, height=650)
xyplot(M.FratioFall ~ year | StudyArea, data=mules1962, type="l")
dev.off()


# 2. Hunting -----------------------------------------------------------------


# Density of animals killed -----------------------------------------------

# Total number of mule deer killed (density / km²) per hunting unit (4A-F)
png("MD_hunted_density_single_units_all.png", width=600, height=500)
ggplot(data=mules1962, aes(x=year, y=d3, colour=unit))+geom_line()+ylab("Mule deer killed / km²")
dev.off()

# Antlered m.d. killed (density / km²) per hunting unit (4A-F)
png("MD_hunted_density_single_units_antlered.png", width=600, height=500)
ggplot(data=mules1962, aes(x=year, y=d1, colour=unit))+geom_line()+ylab("Mule deer killed / km²")
dev.off()

# Antler-less m.d. killed (density / km²) per hunting unit (4A-F)
png("MD_hunted_density_single_units_antlerless.png", width=600, height=500)
ggplot(data=mules1962, aes(x=year, y=d2, colour=unit))+geom_line()+ylab("Mule deer killed / km²")
dev.off()



# Absolute number of animals killed ---------------------------------------

# Total number of mule deer killed per hunting unit (4A-F)
png("MD_hunted_number_single_units_all.png", width=600, height=500)
ggplot(data=mules1962, aes(x=year, y=total_harvest.MDgunHarvest., colour=unit))+geom_line()+ylab("No. mule deer killed")
dev.off()

# Antlered m.d. killed per hunting unit (4A-F)
png("MD_hunted_number_single_units_antlered.png", width=600, height=500)
ggplot(data=mules1962, aes(x=year, y=Antlered_harvest.MDgunHarvest., colour=unit))+geom_line()+ylab("No. mule deer killed")
dev.off()

# Antler-less m.d. killed per hunting unit (4A-F)
png("MD_hunted_number_single_units_antlerless.png", width=600, height=500)
ggplot(data=mules1962, aes(x=year, y=Antleress_Harvest.MDgunHarvest., colour=unit))+geom_line()+ylab("No. mule deer killed")
dev.off()



# Licenses issued ---------------------------------------------------------

# No. hunting licenses issued (for whole area or even the whole of ND?)

# Since these data are the same for all study sites, choose a random site:
data_3V <- subset(mules1962, StudyArea=="3_V")
png("Hunting_licenses_issued.png", width=550, height=500)
plot(Deer_Licenses_Issued~year, data=data_3V, type="b", xlab="Year", ylab="No. licenses issued")
dev.off()






