#' Evaluate the effective number of samples
#' @description This function evaluates the effective number of samples in a dataset.
#'
#' @param input_data [a data frame]
#' A data frame that contains the observational data. Each
#' column corresponds to one variable and each row is a sample that gives the
#' values for all the observed variables. The column names correspond to the
#' names of the observed variables. Data must be discrete like.
#' @param plot [a boolean value] if the autocorrelation plot has to be done. It will be performed only if all values of the correlation vector are positive.
#' @return A list containing the autocorrelation decay, the effective number of samples, and the result of an exponentiality test with alpha = 0.05
#' @export
#' @useDynLib miic

miic.evaluate.effn <- function(input_data = NULL, plot = T) {
  result <- list()
  #### Check the input arguments
  if (is.null(input_data)) {
    stop("The input data file is required")
  }
  if (!is.data.frame(input_data)) {
    stop("The input data is not a dataframe")
  }
  # Remove rows with only NAs
  input_data <- input_data[rowSums(is.na(input_data)) != ncol(input_data), ]
  if (length(input_data) == 0) {
    stop("The input data is empty or contains only NAs")
  }

  if (base::requireNamespace("Rcpp", quietly = TRUE)) {
    cpp_input <- data.frame(t(sapply(input_data, as.character)),
                            stringsAsFactors = FALSE)
    # Call C++ function
    result <- evaluateEffn(cpp_input)
  }
  if (length(which(result$correlation > 0)) == length(result$correlation)) {
    fit1 <- MASS::fitdistr(result$correlation, "exponential")
    pval <- stats::ks.test(result$correlation, "pexp", fit1$estimate)$p.value
    if (pval < 0.05) {
      result$exponential_decay <- FALSE
    } else {
      result$exponential_decay <- TRUE
    }

    if (plot) {
      graphics::par(
        mar = rep(1.5, 4),
        oma = c(3, 3, 3, 1),
        las = 1
      )

      graphics::plot(0:(length(result$correlation) - 1),
        result$correlation,
        type = "l",
        log = "y"
      )
      graphics::title("Autocorrelation between n distant samples",
        ylab = "Autocorrelation with lag",
        xlab = "n"
      )
    }
  } else {
    result$exponential_decay <- FALSE
  }
  result
}
