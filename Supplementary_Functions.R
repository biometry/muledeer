###---------SUPPLEMENTARY FUNCTIONS

library(reshape2)


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

#only problem left: names might not be assigned correctly
# tested via: mean(mules1962$d3[which(mules1962$year == 1962 & mules1962$macrounit == "0-1")])