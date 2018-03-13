#' Discretize two distributions with
#' @description This function discretizes two distributions by maximizing their mutual information, with MDL.
#'
#' @param myDist1 [a vector]
#' A vector that contains the observational data of the first variable.
#' @param myDist2 [a vector]
#' A vector that contains the observational data of the second variable.
#' @param maxbins [an int] The maximum number of bins to test for.
#' @return A list with the two vectors containing the cutpoints of the best discretization for both variables.
#' @export
#' @useDynLib miic

discretizeMutual <- function(myDist1 = NULL, myDist2 = NULL, maxbins=50)
{
  result = list()
  #### Check the input arguments
  if( is.null( myDist1 ) || is.null(myDist2) )
  { stop("The input data file is required") }
  if (base::requireNamespace("Rcpp", quietly = TRUE)) {
    rescpp <- .Call('mydiscretizeMutual', myDist1, myDist2, maxbins, PACKAGE = "miic")
  }
  niterations = length(rescpp$cutpoints1)/maxbins

  result$niterations = niterations
  for(i in 0:(niterations-1)){
    clean_cutpoints1 = rescpp$cutpoints1[(maxbins*i)+(1:maxbins)]
    clean_cutpoints1 = clean_cutpoints1[clean_cutpoints1 != -1]
    clean_cutpoints1 = unique(sort(myDist1)[c(1, clean_cutpoints1+1, length(myDist1))])
    clean_cutpoints2 = rescpp$cutpoints2[(maxbins*i)+(1:maxbins)]
    clean_cutpoints2 = clean_cutpoints2[clean_cutpoints2 != -1]
    clean_cutpoints2 = unique(sort(myDist2)[c(1, clean_cutpoints2+1, length(myDist2))])
    result[[paste0("iteration",i+1)]] = list(cutpoints1=clean_cutpoints1, cutpoints2=clean_cutpoints2)
  }
  result$cutpoints1 = result[[paste0("iteration", niterations)]]$cutpoints1
  result$cutpoints2 = result[[paste0("iteration", niterations)]]$cutpoints2

  result
}
