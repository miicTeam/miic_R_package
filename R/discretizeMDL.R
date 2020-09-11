#' Discretize a real valued distribution
#' @description This function performs minimum description length (MDL)-optimal histogram density estimation
#' as described in Kontkanen and Myllymäki (2007) and returns the cutpoints found to give the best model
#' according to the MDL principle.
#' 
#' @references 
#' \itemize{
#' \item Kontkanen P, Myllymäki P. MDL histogram density estimation. Artificial Intelligence and Statistics 2007 Mar 11 (pp. 219-226).
#' }
#'
#' @param x [a vector]
#' A vector that contains the distribution to be discretized.
#' @param max_bins [an int]
#' The maximum number of bins allowed by the algorithm.
#' @return A list containing the cutpoints of the best discretization.
#'
#' @export
#' @useDynLib miic
#'
#' @examples
#' library(miic)
#' # Bimodal normal distribution
#' N <- 300
#' modes <- sample(1:2, size = N, replace = TRUE)
#' x <- as.numeric(modes == 1) * rnorm(N, mean = 0, sd = 1) +
#'      as.numeric(modes == 2) * rnorm(N, mean = 5, sd = 2)
#' MDL_disc <- discretizeMDL(x)
#' hist(x, breaks = MDL_disc$cutpoints)
#'
#' N <- 2000
#' modes <- sample(1:2, size = N, replace = TRUE)
#' x <- as.numeric(modes == 1) * rnorm(N, mean = 0, sd = 1) +
#'      as.numeric(modes == 2) * rnorm(N, mean = 5, sd = 2)
#' MDL_disc <- discretizeMDL(x)
#' hist(x, breaks = MDL_disc$cutpoints)
#'
discretizeMDL <- function(x = NULL, max_bins = 20) {
  result <- list()
  #### Check the input arguments
  if (is.null(x) || !is.vector(x)) {
    stop("The input data is required and must be a vector")
  }
  if (base::requireNamespace("Rcpp", quietly = TRUE)) {
    result <- mydiscretizeMDL(x, max_bins)
  }
  result
}
