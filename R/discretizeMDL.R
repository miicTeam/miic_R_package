#' Discretize a real valued distribution
#' @description This function performs minimum description length (MDL)-optimal histogram density estimation
#' as described in Kontkanen and Myllymäki (2007) and returns the cutpoints found to give the best model
#' according to the MDL principle.
#'
#' @param myDist [a vector]
#' A vector that contains the distribution to be discretized.
#' @param maxbins [an int]
#' The maximum number of bins allowed by the algorithm.
#' @return A list containing the cutpoints of the best discretization.
#'
#' @export
#' @useDynLib miic
#'
#' @examples
#' library(miic)
#' \dontrun{
#' # Bimodal normal distribution
#' N <- 300
#' modes <- sample(1:2, size = N, replace = T)
#' myDist <- as.numeric(modes == 1) * rnorm(N, mean = 0, sd = 1) + as.numeric(modes == 2) * rnorm(N, mean = 5, sd = 2)
#' MDL_disc <- discretizeMDL(myDist)
#' hist(myDist, breaks = MDL_disc$cutpoints)
#'
#' N <- 2000
#' modes <- sample(1:2, size = N, replace = T)
#' myDist <- as.numeric(modes == 1) * rnorm(N, mean = 0, sd = 1) + as.numeric(modes == 2) * rnorm(N, mean = 5, sd = 2)
#' MDL_disc <- discretizeMDL(myDist)
#' hist(myDist, breaks = MDL_disc$cutpoints)
#' }
#'
discretizeMDL <- function(myDist = NULL, maxbins = 20) {
  result <- list()
  #### Check the input arguments
  if (is.null(myDist)) {
    stop("The input data file is required")
  }
  if (base::requireNamespace("Rcpp", quietly = TRUE)) {
    result <- mydiscretizeMDL(myDist, maxbins)
  }
  result
}
