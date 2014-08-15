### Data Exploration: Collinearities----

#Check for collinearity between predictors
source("HighstatLibV6_correlation_functions.R") # Replacement for AEV-package of Zuur et al
z <- AllMeans[,!(names(AllMeans) %in% c("macrounit", "MDperKMsqSpring_mean","WTailDen_mean","FawnFall_mean","MaleFall_mean","FemaleFall_mean"))] 

par(oma=c(2,0,2,0))
pairs(z, lower.panel = panel.smooth2, upper.panel = panel.cor, diag.panel = panel.hist) #not too informative
title("Pairwise Pearson Correlation", outer=TRUE)

sink("cormatrices.txt", type="output")
cormatrix <- cor(z, use="pairwise.complete.obs") # no correlation >(-)0.4 so ok
cormatrix1 <- cor(z[which(AllMeans$macrounit == "0-1"),], use="pairwise.complete.obs") # no correlation >(-)0.4 so ok
cormatrix2 <- cor(z[which(AllMeans$macrounit == "0-2"),], use="pairwise.complete.obs") # no correlation >(-)0.4 so ok
cormatrix3 <- cor(z[which(AllMeans$macrounit == "0-3"),], use="pairwise.complete.obs") # no correlation >(-)0.4 so ok
cormatrix4 <- cor(z[which(AllMeans$macrounit == "0-4"),], use="pairwise.complete.obs") # no correlation >(-)0.4 so ok
print("Correlation Matrix All Macrounits")
print(cormatrix)
print("Correlation Matrix Macrounit 1")
print(cormatrix1)
print("Correlation Matrix Macrounit 2")
print(cormatrix2)
print("Correlation Matrix Macrounit 3")
print(cormatrix3)
print("Correlation Matrix Macrounit 4")
print(cormatrix4)
sink(NULL)

pairs(z, lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Whole Area")
pairs(z[which(AllMeans$macrounit == "0-1"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 1")
pairs(z[which(AllMeans$macrounit == "0-2"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 2")
pairs(z[which(AllMeans$macrounit == "0-3"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 3")
pairs(z[which(AllMeans$macrounit == "0-4"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 4")

vif <- corvif(z)# all values <10, so ok


# Collinearity between PopDen in MUs (as to see wether they differ much)

z <- data.frame(AllMeans[which(AllMeans$macrounit == "0-1"),"MDperKMsqSpring_mean"])
z <- cbind(z, AllMeans[which(AllMeans$macrounit == "0-2"),"MDperKMsqSpring_mean"],AllMeans[which(AllMeans$macrounit == "0-3"),"MDperKMsqSpring_mean"],AllMeans[which(AllMeans$macrounit == "0-4"),"MDperKMsqSpring_mean"])
names(z) <- c("0-1","0-2","0-3","0-4")
pairs(z, lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearities between macrounits")


#### Hunting licences correlation Hunt <-> RepRate(t+1) and Fawn:Female Ratio as information on distribution practice of licenses---------
z <- AllMeans[,names(AllMeans) %in% c("HuntDen_All_mean", "RepRateFall_mean_tminus1", "FawnFemaleRatio_mean_tminus1", "WoodyVeg_mean", "WellDen_mean", "AvrgWinterMinTemp_tminus1", "CoyoteDen_mean")] 

pairs(z, lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Whole Area")
pairs(z[which(AllMeans$macrounit == "0-1"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 1")
pairs(z[which(AllMeans$macrounit == "0-2"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 2")
pairs(z[which(AllMeans$macrounit == "0-3"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 3")
pairs(z[which(AllMeans$macrounit == "0-4"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 4")

#Estimates of Population Densities and growth------------
z <- AllMeans[,names(AllMeans) %in% c("MDperKMsqFall_mean", "MDperKMsqSpring_mean", "RepRateFall_mean_tminus1", "FawnFemaleRatio_mean_tminus1")] 

pairs(z, lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Whole Area")
pairs(z[which(AllMeans$macrounit == "0-1"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 1")
pairs(z[which(AllMeans$macrounit == "0-2"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 2")
pairs(z[which(AllMeans$macrounit == "0-3"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 3")
pairs(z[which(AllMeans$macrounit == "0-4"),], lower.panel = panel.smooth, diag.panel = panel.hist, upper.panel = panel.cor, main = "Collinearity Analysis Macrounit 4")
