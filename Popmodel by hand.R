log_tdh <- function(par, suppar){
  output <- numeric(suppar[2])#steps
  output[1] <- suppar[1]#n0
  for (t in 1:(length(output)-1)){
    output[t+1] <- output[t]+(output[t]*par[1]*(1-((output[t]/par[2])^par[3])))-suppar[t+2]#from suppar 3 on: harvest data
  }
  return(output)
}

a <- Popdata$MDperKMsqFall_mean[1]
b <- a+a*2.00813824*(1-((a/1.71814299)^0.03431441))-Popdata$HuntDen_All_mean[1]
c <- b+b*2.00813824*(1-((b/1.71814299)^0.03431441))-Popdata$HuntDen_All_mean[2]#produces NA
# 
(b/1.71814299)# = -0.1030193
(-0.1030193)^0.03431441#=-0.9249726
# but it does not work with b --> ??

plot(Popdata$MDperKMsqFall_mean[cond]~Popdata$year[cond], type="l")
lines(Popdata$HuntDen_All_mean[cond]~Popdata$year[cond], col="red")

?exp
