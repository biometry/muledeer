#interaction temp+hunt on growth rate (= the two sole winter effects) -> useless as growth rate is affected by all factors over the whole year
gam_hunttemp <- gam(RepRateFall_mean ~ te(HuntDen_All_mean,AvrgWinterMinTemp,bs="cs"), data=AllMeans)#no
gam_hunttemp <- gam(RepRateFall_mean ~ s(MDperKMsqFall_mean, bs="cs") + te(HuntDen_All_mean,AvrgWinterMinTemp, bs="cs", by=macrounit) + macrounit, data=AllMeans)#no
gam_hunttemp <- gam(RepRateFall_mean ~ s(MDperKMsqFall_mean, by=macrounit, bs="cs") +te(HuntDen_All_mean,AvrgWinterMinTemp, bs="cs") + macrounit, data=AllMeans)#no


# +++Winter effects ----
### Winter effects on Fawn Female Ratio (coyote,Temp and Hunt and interactions) - same as in Simones model?---------------

#without Coyote
gam_temp <- gam(FawnFemaleRatio_mean ~ s(AvrgWinterMinTemp, bs="cs"), data=AllMeans)#yes
gam_temp <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +s(AvrgWinterMinTemp, bs="cs"), data=AllMeans)#yes
gam_temp <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +s(AvrgWinterMinTemp, by=macrounit, bs="cs")+macrounit, data=AllMeans)#only MU1
gam_temp <- gam(FawnFemaleRatio_mean ~ s(AvrgWinterMinTemp, by=macrounit, bs="cs")+macrounit, data=AllMeans)#only MU1

#Temp*Coyote
gam_coyotettemp <- gam(FawnFemaleRatio_mean ~ te(CoyoteDen_mean,AvrgWinterMinTemp,bs="cs"), data=AllMeans)#yes
gam_coyotettemp <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +te(CoyoteDen_mean,AvrgWinterMinTemp,bs="cs"), data=AllMeans)#yes
gam_coyotettemp <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +te(CoyoteDen_mean,AvrgWinterMinTemp,bs="cs", by=macrounit)+macrounit, data=AllMeans)#None
gam_coyotettemp <- gam(FawnFemaleRatio_mean ~ te(CoyoteDen_mean,AvrgWinterMinTemp,bs="cs", by=macrounit)+macrounit, data=AllMeans)#only MU1

# Hunt only
gam_hunt <- gam(FawnFemaleRatio_mean ~ s(HuntDen_All_mean, bs="cs"), data=AllMeans)#nope
gam_hunt <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +s(HuntDen_All_mean, bs="cs"), data=AllMeans)#no
gam_hunt <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +s(HuntDen_All_mean, by=macrounit, bs="cs")+macrounit, data=AllMeans)#no
gam_hunt <- gam(FawnFemaleRatio_mean ~ s(HuntDen_All_mean, by=macrounit, bs="cs")+macrounit, data=AllMeans)#only MU1

#Temp*Hunt
gam_hunttemp <- gam(FawnFemaleRatio_mean ~ te(HuntDen_All_mean,AvrgWinterMinTemp,bs="cs"), data=AllMeans)#yes!!
gam_hunttemp <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +te(HuntDen_All_mean,AvrgWinterMinTemp,bs="cs"), data=AllMeans)#somewhat:0.0354 *
gam_hunttemp <- gam(FawnFemaleRatio_mean ~ s(year, bs="cs") +te(HuntDen_All_mean,AvrgWinterMinTemp,bs="cs", by=macrounit)+macrounit, data=AllMeans)#MU1+MU4
gam_hunttemp <- gam(FawnFemaleRatio_mean ~ te(HuntDen_All_mean,AvrgWinterMinTemp,bs="cs", by=macrounit)+macrounit, data=AllMeans)#only MU1:0.00456 **

summary(gam_temp)
summary(gam_hunt)
summary(gam_hunttemp)
AIC(gam_temp)
plot(gam_temp)
gam_temppred<-predict(gam_temp, se.fit=T, newdata=data.frame(year=AllMeans$year, AvrgWinterMinTemp=AllMeans$AvrgWinterMinTemp), type="response")
plot(AllMeans$FawnFemaleRatio_mean~AllMeans$year, cex=0.5, main="GAM Fawn Female Ratio ~ Year and Winter Temp")
lines(gam_temppred$fit[1:51]~AllMeans$year[1:51],type="l", col="red", )

### Winter effects on Winter mortality

gam_temp <- gam(WinterMort_mean ~ s(AvrgWinterMinTemp_tplus1,  bs="cs"), data=AllMeans)#NO EFFECT
gam_hunt <- gam(WinterMort_mean ~ s(HuntDen_All_mean,  bs="cs"), data=AllMeans)#NO EFFECT
gam_hunttemp <- gam(WinterMort_mean ~ te(AvrgWinterMinTemp_tplus1, HuntDen_All_mean,  bs="cs"), data=AllMeans)#NO EFFECT
gam_hunttemp <- gam(WinterMort_mean ~ te(AvrgWinterMinTemp_tplus1, HuntDen_All_mean,  bs="cs",by=macrounit)+macrounit, data=AllMeans)# MU4!

