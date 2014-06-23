ex <- AllMeans[which(AllMeans$macrounit == "0-1"),c(-MDpersqKmSpring,-CoyoteDen,-WellDen, -Fw.FratioFall)]
write.csv(AllMeans, "data.csv")

# no autocorrelation included
gam_all2 <- gam(MDperKMsqFall_mean ~ s(year, by=macrounit) + macrounit, data=AllMeans)

xx <- read.csv("data.csv")
