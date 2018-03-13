#' Discretize a real valued distribution
#' @description This function discretizes a distribution with MDL.
#'
#' @param myDist [a vector]
#' A vector that contains the observational data.
#' @param maxbins [an int] The maximum number of bins to test for.
#' @return A vector containing the cutpoints of the best discretization according to the MDL principle.
#' @export
#' @useDynLib miic

discretizeMDL <- function(myDist = NULL, maxbins=20)
{
  result = list()
  #### Check the input arguments
  if( is.null( myDist ) )
  { stop("The input data file is required") }
  if (base::requireNamespace("Rcpp", quietly = TRUE)) {
    result <- .Call('mydiscretizeMDL', myDist, maxbins, PACKAGE = "miic")
  }
  result
}
