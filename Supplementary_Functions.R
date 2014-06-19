###---------SUPPLEMENTARY FUNCTIONS

library(reshape2)


###---------Count NAs
count_nas <- function(x) length(which(is.na(x)))

###---------Data Extraction
extract <- function(data, columns=NULL, fun = NULL, xvar, listvar){
  if (is.null(fun)) { #only raw data?
    a <- subset(data, select=columns)
    return(a)}
  else {
    indexlist <- list()
      for (i in 1:length(listvar)){ # tapply needs as a list
        indexlist[[i]] <- data[,listvar[i]]
      }
    b <- t(tapply(X=data[,xvar[1]], INDEX = indexlist, FUN=fun, na.rm = T))
    b <- melt(b)
      if (length(xvar) > 1){ # if several columns should be processed
        for (i in 2:length(xvar)){
          c <- t(tapply(X=data[,xvar[i]], INDEX = indexlist, FUN=fun, na.rm = T))
          c <- melt(c)
          b <- cbind(b, c[,NCOL(c)])
        }
      }
    names(b) <- names(data[,c(rev(listvar),xvar)])# for some reason, tapply always reverses the order of listvar    
    return(b)}    
}

#problems left: 
#1. names might not be assigned correctly
#2. if one of the xvars is shorter than the one processed before because ofna.rm, that cell is still assigned NA
#   also if there are rows missing due to NAs in first column -> values of other columns might be lost
#   SOLUTION: y(including NAs)~AllMeans$year[which(!is.na(AllMeans$y))]

# tested via: mean(mules1962$d3[which(mules1962$year == 1962 & mules1962$macrounit == "0-1")])

###-------- One plot per macrounit displaying gam-object

macrounitplots <- function(glmobject, title, colour){
  par(mfrow=c(2,2),oma=c(2,0,2,0))
  macrounits <- levels(glmobject[,"macrounit"])
  for (i in 1:length(macrounits)){
    cond = which(glmobject[,"macrounit"]==macrounits[i])  
    plot(AllMeans[cond,"MDperKMsqSpring_mean"]~AllMeans[cond,"year"], type="p", main=macrounits[i], xlab="Year", ylab="Density per km²")
    
    lines(x=AllMeans[cond,"year"], glmobject[cond,"fit"], col=colour)
    lines(x=AllMeans[cond,"year"], glmobject[cond,"fit"]  + 2 * glmobject[cond,"se.fit"], col=colour, lty=2)
    lines(x=AllMeans[cond,"year"], glmobject[cond,"fit"]  - 2 * glmobject[cond,"se.fit"], col=colour, lty=2)
  }
  title(as.character(title), outer=TRUE)
}