MIGaussian <- function(myDist1 = NULL, myDist2 = NULL, uis=NULL){
  library(ppcor)
  if( is.null( myDist1 ) || is.null(myDist2) )
  { stop("The input data file is required") }
  if(length(myDist1) != length(myDist2)){
    stop("myDist1 vector must be of the same length of myDist2")
  }
  if(!is.null(uis)){
    if(length(myDist1) != nrow(uis)){
      stop("uis data frame must have the same numbre of samples of myDist1 and myDist2")
    }
  }
  data = cbind(myDist1,myDist2)
  if(!is.null(uis)){
    data = cbind(data,uis)
  }
  cor = ppcor::pcor(data)$estimate[1,2]
  MI  = (-log(1 - cor * cor)/2) * length(myDist1);
  return(MI)
}
