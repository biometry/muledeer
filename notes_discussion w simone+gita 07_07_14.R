

FawnFemaleRatio_mean = (AllMeans$FawnFall_mean/AllMeans$FemaleFall_mean)
AllMeans <- cbind(AllMeans, FawnFemaleRatio_mean)
gam_temp <- gam(FawnFemaleRatio_mean ~ s(AvrgWinterMinTemp, bs="cs"), data=AllMeans)

gam_temp <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +s(AvrgWinterMinTemp, bs="cs"), data=AllMeans)
summary(gam_temp)
head(AllMeans)
AIC(gam_temp)
Spring_tplus1 <- c(NA, AllMeans$MDperKMsqSpring_mean[-length(AllMeans$MDperKMsqSpring_mean)])

plot(Spring_tplus1~AllMeans$MDperKMsqFall_mean)
ratio <- Spring_tplus1/AllMeans$MDperKMsqFall_mean
plot(ratio~AllMeans$HuntDen_All_mean, log="y")
plot(ratio~AllMeans$AvrgWinterMinTemp, log="y")

gam_temp <- gam(ratio ~ s(AvrgWinterMinTemp,  bs="cs"), data=AllMeans)#2
gam_hunt <- gam(ratio ~ s(HuntDen_All_mean,  bs="cs"), data=AllMeans)#2

summary(gam_hunt)
plot(gam_hunt)
AIC(gam_temp)

HuntDen_all_mean_tplus1 <- c(NA, AllMeans$HuntDen_All_mean[-length(AllMeans$HuntDen_All_mean)])
AllMeans$HuntDen_All_mean

cor(AllMeans$MDperKMsqSpring_mean, AllMeans$MDperKMsqFall_mean, use="pairwise.complete.obs")
cor(Spring_tplus1, AllMeans$MDperKMsqFall_mean, use="pairwise.complete.obs")
cor(FawnFemaleRatio_mean,Spring_tplus1,use="pairwise.complete.obs")
cor(FawnFemaleRatio_mean,HuntDen_all_mean_tplus1,use="pairwise.complete.obs")
cor(ratio,HuntDen_all_mean_tplus1, use="pairwise.complete.obs")

cor(HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-1")],HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-2")], use="pairwise.complete.obs")
cor(HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-1")],HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-3")], use="pairwise.complete.obs")
cor(HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-1")],HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-4")], use="pairwise.complete.obs")
cor(HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-2")],HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-3")], use="pairwise.complete.obs")
cor(HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-3")],HuntDen_all_mean_tplus1[which(AllMeans$macrounit == "0-4")], use="pairwise.complete.obs")


?cor
