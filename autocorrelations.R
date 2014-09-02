
### AUTOCORRELATION IN DATA ----
par(mfrow=c(2,3))
# AllMeans
acf(AllMeans$MDperKMsqFall_mean, na.action = na.pass,main = "Data", ylab="", xlab="",cex.lab=4)

# WholeAreaMeans
acf(WholeAreaMeans$MDperKMsqFall_mean, na.action = na.pass,main = "PopDen WholeAreaMeans")

# AllMeans by macorunit
macrounits <- levels(AllMeans$macrounit)
for (i in 1:length(macrounits)){
  cond = which(AllMeans$macrounit==macrounits[i])  
  acf(AllMeans$MDperKMsqFall_mean[cond], na.action = na.pass,main = c("PopDen AllMeans macrounit ", i))
}  


### AUTOCORRELATION IN RESIDUALS ----


#Nullmodels
acf(gam_all0res, na.action = na.pass,main = "AllMeans:PopDen ~ s(year)")

acf(gam_all3res, na.action = na.pass,main = "WholeAreaMeans:PopDen ~ s(year)")

macrounits <- levels(AllMeans$macrounit)
for (i in 1:length(macrounits)){
  cond = which(AllMeans$macrounit==macrounits[i])  
  acf(gam_all2res[cond], na.action = na.pass,main = c("AllMeans:PopDen ~ s(year, by=macrounit) , macrounit ",i))
}

#Univariate models ("see-Gam-explanatory.R")

acf(gam_woodyvegres, na.action = na.pass,main = "PopDen ~ s(year) + s(WoodyVeg)")
acf(gam_tempres, na.action = na.pass,main = "PopDen ~ s(AvrgWinMinTemp)")
acf(gam_tempres, na.action = na.pass,main = "PopDen ~ s(year) + s(AvrgWinMinTemp)")
acf(gam_huntres, na.action = na.pass,main = "PopDen ~ s(HuntDen)")
acf(gam_huntres, na.action = na.pass,main = "PopDen ~ s(year) + s(HuntDen)")
acf(gam_oilres, na.action = na.pass,main = "PopDen ~ s(WellDen)")
acf(gam_oilres, na.action = na.pass,main = "PopDen ~ s(year) + s(WellDen)")
acf(gam_coyoteres, na.action = na.pass,main = "PopDen ~ s(CoyoteDen)")
acf(gam_coyoteres, na.action = na.pass,main = "PopDen ~ s(year) + s(CoyoteDen)")
acf(gam_ffratiores, na.action = na.pass,main = "PopDen ~ s(FFRatio)")
acf(gam_ffratiores, na.action = na.pass,main = "PopDen ~ s(year) + s(FFRatio)")