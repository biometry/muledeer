library(lattice)
xyplot(AllMeans$MDperKMsqFall_mean + AllMeans$HuntDen_All_mean ~ year | macrounit, data=AllMeans, subscripts=TRUE, type="l", auto.key = list(space = "top", text = c("Population Density", "Hunting Density")), xlab = "Year", ylab= "Density per km²", main = "Fall Population Density and Hunting Density")
