#' Evaluate the effective number of samples
#' @description This function evaluates the effective number of samples in a dataset.
#'
#' @param inputData [a data frame]
#' A data frame that contains the observational data. Each
#' column corresponds to one variable and each row is a sample that gives the
#' values for all the observed variables. The column names correspond to the
#' names of the observed variables. Data must be discrete like.
#' @param plot [a boolean value] if the autocorrelation plot has to be done. It will be performed only if all values of the correlation vector are positive.
#' @return A list containing the autocorrelation decay, the effective number of samples, and the result of an exponentiality test with alpha = 0.05
#' @export
#' @useDynLib miic

miic.evaluate.effn <- function(inputData = NULL, plot=T)
{
  result = list()
  #### Check the input arguments
  if( is.null( inputData ) )
  { stop("The input data file is required") }
  inData <- c(colnames(inputData), as.vector(as.character(t(as.matrix(inputData)))))
  if (base::requireNamespace("Rcpp", quietly = TRUE)) {
    result <- .Call('evaluateEffn', inData, ncol(inputData), nrow(inputData),PACKAGE = "miic")
  }
  if(length(which(result$correlation > 0)) == length(result$correlation)){

    fit1 <- MASS::fitdistr(result$correlation, "exponential")
    pval = stats::ks.test(result$correlation, "pexp", fit1$estimate)$p.value
    if(pval < 0.05){
      result$exponential_decay= FALSE
    } else {
      result$exponential_decay= TRUE
    }

    if(plot){
      graphics::par(mar=rep(1.5, 4), oma=c(3,3,3,1), las=1)

      graphics::plot(0:(length(result$correlation)-1),result$correlation, type="l", log="y")
      graphics::title("Autocorrelation between n distant samples",ylab="Autocorrelation with lag", xlab="n")
    }
  }else {
    result$exponential_decay= FALSE
  }
   result
}
