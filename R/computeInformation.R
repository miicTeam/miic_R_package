#' Compute (conditional) mutual information
#' @description For discrete variables, the computation is based on the
#' empirical frequency minus a complexity cost (computed as BIC or with the
#' Normalized Maximum Likelihood). When continuous variables are present, each
#' continuous variable is discretized where the partitioning is chosen by
#' maximizing the mutual information minus the complexity cost. The estimation
#' based on the optimally discretized distributions effectively approaches the
#' mutual information computed on the original continuous variables.
#'
#' @details For a pair of continuous variables \eqn{X} and \eqn{Y}, the mutual
#' information \eqn{I(X;Y)} will be computed iteratively. In each iteration, the
#' algorithm optimizes first the partitioning of \eqn{X} and then that of
#' \eqn{Y}, while maximizing
#' \deqn{Ik(X_{d};Y_{d}) = I(X_{d};Y_{d}) - cplx(X_{d};Y_{d})}
#' where \eqn{cplx(X_{d}; Y_{d})} is the complexity cost of the current
#' partitioning (see Affeldt 2016 and Cabeli 2020). Upon convergence, the
#' information terms \eqn{I(X_{d};Y_{d})} and \eqn{Ik(X_{d};Y_{d})}, as well as
#' the partitioning of \eqn{X_{d}} and \eqn{Y_{d}} in terms of cutpoints, are
#' returned.
#'
#' For conditional mutual information with conditioning set \eqn{U}, the
#' computation is done based on
#' \deqn{
#'   Ik(X;Y|U) = 0.5*(Ik(X_{d};Y_{d},U_{d}) - Ik(X_{d};U_{d})
#'                  + Ik(Y_{d};X_{d},U_{d}) - Ik(Y_{d};U_{d})),
#' }
#' where each of the four summands is estimated independently.
#'
#' @references
#' \itemize{
#' \item Verny et al., \emph{PLoS Comp. Bio. 2017.}
#'   https://doi.org/10.1371/journal.pcbi.1005662
#' \item Cabeli et al., \emph{PLoS Comp. Bio. 2020.}
#'   https://doi.org/10.1371/journal.pcbi.1007866
#' \item Affeldt et al., \emph{Bioinformatics 2016}
#' }
#'
#' @param X [a vector]
#' A vector that contains the observational data of the first variable.
#' @param Y [a vector]
#' A vector that contains the observational data of the second variable.
#' @param df_conditioning [a data frame]
#' The data frame of the observations of the conditioning variables.
#' @param maxbins [an integer]
#' When the data contain continuous variables, the maximum number of bins
#' allowed during the discretization. A smaller number makes the computation
#' faster, a larger number allows finer discretization.
#' @param cplx [a string]
#' The complexity model:
#' \itemize{
#' \item["mdl"] Minimum description Length
#' \item["nml"] Normalized Maximum Likelihood, less costly compared to "mdl" in
#' the finite sample case and will allow for more bins.
#' }
#' @param n_eff [an integer]
#' The number of effective samples. When there is significant autocorrelation in
#' the samples you may want to specify a number of effective samples that is
#' lower than the number of points in the distribution.
#' @param sample_weights [a vector of floats]
#' Individual weights for each sample, used for the same reason as the effective
#' sample number but with individual precision.
#' @param is_continuous [a vector of booleans]
#' Specify if each variable is to be treated as continuous (TRUE) or discrete
#' (FALSE), must be of length `ncol(df_conditioning) + 2`, in the order
#' \eqn{X, Y, U1, U2, ...}. If not specified, factors and character vectors are
#' considered as discrete, and numerical vectors as continuous.
#' @param plot [a boolean]
#' Specify whether the XY joint space with discretization scheme is to be
#' plotted (requires `ggplot2` and `gridExtra`).
#'
#' @return A list that contains :
#' \itemize{
#' \item cutpoints1: Only when \eqn{X} is continuous, a vector containing
#'   the cutpoints for the partitioning of \eqn{X}.
#' \item cutpoints2: Only when \eqn{Y} is continuous, a vector containing
#'   the cutpoints for the partitioning of \eqn{Y}.
#' \item niterations: Only when at least one of the input variables is
#'   continuous, the number of iterations it takes to reach the convergence of
#'   the estimated information.
#' \item iterationN: Only when at least one of the input variables is
#'   continuous, the list of vectors of cutpoints of each iteration.
#' \item info: The estimation of (conditional) mutual information without the
#' complexity cost.
#' \item infok: The estimation of (conditional) mutual information with the
#' complexity cost (\eqn{Ik = I - cplx}).
#' \item plot: Only when `plot == TRUE`, the plot object.
#' }
#' @export
#' @useDynLib miic
#' @importFrom stats density sd
#'
#' @examples
#' library(miic)
#' N <- 1000
#' # Dependence, conditional independence : X <- Z -> Y
#' Z <- runif(N)
#' X <- Z * 2 + rnorm(N, sd = 0.2)
#' Y <- Z * 2 + rnorm(N, sd = 0.2)
#' res <- computeMutualInfo(X, Y, plot = FALSE)
#' message("I(X;Y) = ", res$info)
#' res <- computeMutualInfo(X, Y, df_conditioning = matrix(Z, ncol = 1), plot = FALSE)
#' message("I(X;Y|Z) = ", res$info)
#'
#' \donttest{
#' # Conditional independence with categorical conditioning variable : X <- Z -> Y
#' Z <- sample(1:3, N, replace = TRUE)
#' X <- -as.numeric(Z == 1) + as.numeric(Z == 2) + 0.2 * rnorm(N)
#' Y <- as.numeric(Z == 1) + as.numeric(Z == 2) + 0.2 * rnorm(N)
#' res <- miic::computeMutualInfo(X, Y, cplx = "nml")
#' message("I(X;Y) = ", res$info)
#' res <- miic::computeMutualInfo(X, Y, matrix(Z, ncol = 1), is_continuous = c(TRUE, TRUE, FALSE))
#' message("I(X;Y|Z) = ", res$info)
#'
#'
#' # Independence, conditional dependence : X -> Z <- Y
#' X <- runif(N)
#' Y <- runif(N)
#' Z <- X + Y + rnorm(N, sd = 0.1)
#' res <- computeMutualInfo(X, Y, plot = TRUE)
#' message("I(X;Y) = ", res$info)
#' res <- computeMutualInfo(X, Y, df_conditioning = matrix(Z, ncol = 1), plot = TRUE)
#' message("I(X;Y|Z) = ", res$info)
#' }
#'
computeMutualInfo <- function(X, Y,
                              df_conditioning = NULL,
                              maxbins = NULL,
                              cplx = c("nml", "mdl"),
                              n_eff = -1,
                              sample_weights = NULL,
                              is_continuous = NULL,
                              plot = FALSE) {
  cplx <- tryCatch(
    {match.arg(cplx)},
    error = function(e) {
      if (grepl("object .* not found", e$message)) {
        message(e, "")
        return("")
      }
      return(toString(cplx))
    }
  )
  cplx <- match.arg(cplx)

  input_data = data.frame(X, Y)
  if (!is.null(df_conditioning)) {
    input_data <- data.frame(input_data, df_conditioning)
  }

  if (!is.null(sample_weights) && length(sample_weights) != nrow(input_data)) {
    stop(paste(
      "Differing number of rows between `sample_weights` and input data:",
      length(sample_weights),
      length(X)
    ))
  }

  complete_row <- stats::complete.cases(input_data)
  n_rows_na <- sum(!complete_row)
  if (n_rows_na > 0) {
    input_data <- input_data[complete_row, ]
    warning(paste0(
      "Removed ", n_rows_na, " rows containing at least one NA value."
    ))
  }

  n_samples <- nrow(input_data)
  n_nodes <- ncol(input_data)

  if (n_samples < 3) {
    stop(paste0("Insufficient number of complete rows: ", nrow(input_data)))
  }

  if (is.null(is_continuous)) {
    is_continuous <- sapply(input_data, is.numeric)
  } else if (length(is_continuous) != n_nodes) {
    stop(paste(
      "Length of `is_continuous` does not match number of input variables:",
      length(is_continuous),
      n_nodes
    ))
  }

  # Numeric factor matrix, level starts from 0
  input_factor <- as.matrix(apply(input_data, 2,
    function(x) (as.numeric(factor(x, levels = unique(x))) - 1))
  )
  max_level_list <- as.numeric(apply(input_factor, 2, max)) + 1
  # Data list, numeric for continuous columns, -1 for discrete columns
  input_double <- matrix(nrow = n_samples, ncol = n_nodes)
  # Order list, order(column) for continuous columns (index starting from 0),
  # -1 for discrete columns
  input_order <- matrix(nrow = n_samples, ncol = n_nodes)
  for (i in c(1: n_nodes)) {
    if (is_continuous[i]) {
      input_double[, i] <- as.numeric(input_data[, i])
      input_order[, i] <- order(input_data[, i], na.last=NA) - 1
    } else {
      input_double[, i] <- rep_len(-1, n_samples)
      input_order[, i] <- rep_len(-1, n_samples)
    }
  }

  arg_list <- list(
    "cplx" = cplx,
    "is_continuous" = is_continuous,
    "levels" = max_level_list,
    "n_eff" = n_eff,
    "n_nodes" = n_nodes,
    "n_samples" = n_samples
  )
  # Continuous variables will be discretized during the computation
  if (any(is_continuous)) {
    initbins <- min(30, round(n_samples**(1 / 3)))
    if (is.null(maxbins) || maxbins > n_samples || maxbins < initbins) {
      maxbins <- min(n_samples, 5 * initbins, 50)
    }
    arg_list[["max_bins"]] <- maxbins
  }
  if (!is.null(sample_weights)) {
    arg_list[["sample_weights"]] <- sample_weights[complete_row, ]
  }
  cpp_input <- list(
    "factor" = as.vector(input_factor),
    "double" = as.vector(input_double),
    "order" = as.vector(input_order)
  )
  # Call cpp code
  rescpp <- mydiscretizeMutual(cpp_input, arg_list)

  result <- list()
  result$info <- rescpp$info
  result$infok <- rescpp$infok

  X_num <- if (is_continuous[1]) input_double[, 1] else input_factor[, 1]
  Y_num <- if (is_continuous[2]) input_double[, 2] else input_factor[, 2]

  if (any(is_continuous)) {
    # Parse cutpointsmatrix
    epsilon <- min(c(sd(X_num), sd(Y_num))) / 100
    niterations <- nrow(rescpp$cutpointsmatrix) / maxbins
    result$niterations <- niterations
    for (i in 0:(niterations - 1)) {
      result[[paste0("iteration", i + 1)]] <- list()
      for (l in 1:2) {
        if (!is_continuous[l]) next

        data <- if (l == 1) X_num else Y_num
        clean_cutpoints <- rescpp$cutpointsmatrix[, l][(maxbins*i) + (1:maxbins)]
        clean_cutpoints <- clean_cutpoints[clean_cutpoints != -1]
        clean_cutpoints <- sort(data)[clean_cutpoints + 1]

        uniquedata <- sort(unique(data))
        if (length(clean_cutpoints) > 0) {
          # Take midpoints between two consecutive unique values instead of
          # the values themselves
          clean_cutpoints <- sapply(clean_cutpoints, function(x) {
            if (x < uniquedata[length(uniquedata)]) {
              return((min(uniquedata[uniquedata > x]) +
                max(uniquedata[uniquedata <= x])) / 2)
            } else {
              return(x)
            }
          })
        }
        clean_cutpoints <- c(uniquedata[1] - epsilon, clean_cutpoints)
        if (max(clean_cutpoints) < uniquedata[length(uniquedata)]) {
          clean_cutpoints <- c(
            clean_cutpoints,
            uniquedata[length(uniquedata)] + epsilon
          )
        }
        result[[paste0("iteration", i + 1)]][[paste0("cutpoints", l)]] <-
          clean_cutpoints
      }
    }
    for (l in 1:n_nodes) {
      result[[paste0("cutpoints", l)]] <-
        result[[paste0("iteration", niterations)]][[paste0("cutpoints", l)]]
    }
  }

  if (plot) {
    nameDist1 <- deparse(substitute(X))
    nameDist2 <- deparse(substitute(Y))
    if (base::requireNamespace("ggplot2", quietly = TRUE) &&
        base::requireNamespace("gridExtra", quietly = TRUE)) {
      if (all(is_continuous[1:2])) {
        result$plot <- jointplot_hist(X_num, Y_num, result, nameDist1, nameDist2)
      } else if (any(is_continuous[1:2])) {
        result$plot <- barplot_disc(
          input_data[, 1],
          input_data[, 2],
          result,
          !is_continuous,
          nameDist1,
          nameDist2
        )
      } else {
        result$plot <- grid_plot(
          input_data[, 1],
          input_data[, 2],
          nameDist1,
          nameDist2
        )
      }
    } else {
      warning("Plotting requires ggplot2 and gridExtra.")
    }
  }

  return(result)
}

#' Compute (conditional) three-point information
#' @description Three point information is defined based on mutual information.
#' For discrete variables, the computation is based on the
#' empirical frequency minus a complexity cost (computed as BIC or with the
#' Normalized Maximum Likelihood). When continuous variables are present, each
#' continuous variable is discretized where the partitioning is chosen by
#' maximizing the mutual information minus the complexity cost.
#'
#' @details For variables \eqn{X}, \eqn{Y}, \eqn{Z} and a set of conditioning
#' variables \eqn{U}, the conditional three point information is defined as
#' \deqn{Ik(X;Y;Z|U) = Ik(X;Y|U) - Ik(X;Y|U,Z)}, where \eqn{Ik} is the
#' regularized conditional mutual information.
#' See \code{\link{computeMutualInfo}} for the definition of \eqn{Ik}.
#'
#' @references
#' \itemize{
#' \item Verny et al., \emph{PLoS Comp. Bio. 2017.}
#'   https://doi.org/10.1371/journal.pcbi.1005662
#' \item Cabeli et al., \emph{PLoS Comp. Bio. 2020.}
#'   https://doi.org/10.1371/journal.pcbi.1007866
#' \item Affeldt et al., \emph{Bioinformatics 2016}
#' }
#'
#' @param X [a vector]
#' A vector that contains the observational data of the first variable.
#' @param Y [a vector]
#' A vector that contains the observational data of the second variable.
#' @param Z [a vector]
#' A vector that contains the observational data of the third variable.
#' @param df_conditioning [a data frame]
#' The data frame of the observations of the set of conditioning variables
#' \eqn{U}.
#' @param maxbins [an integer]
#' When the data contain continuous variables, the maximum number of bins
#' allowed during the discretization. A smaller number makes the computation
#' faster, a larger number allows finer discretization.
#' @param cplx [a string]
#' The complexity model:
#' \itemize{
#' \item["mdl"] Minimum description Length
#' \item["nml"] Normalized Maximum Likelihood, less costly compared to "mdl" in
#' the finite sample case and will allow for more bins.
#' }
#' @param n_eff [an integer]
#' The number of effective samples. When there is significant autocorrelation in
#' the samples you may want to specify a number of effective samples that is
#' lower than the number of points in the distribution.
#' @param sample_weights [a vector of floats]
#' Individual weights for each sample, used for the same reason as the effective
#' sample number but with individual precision.
#' @param is_continuous [a vector of booleans]
#' Specify if each variable is to be treated as continuous (TRUE) or discrete
#' (FALSE), must be of length `ncol(df_conditioning) + 2`, in the order
#' \eqn{X, Y, U1, U2, ...}. If not specified, factors and character vectors are
#' considered as discrete, and numerical vectors as continuous.
#'
#' @return A list that contains :
#' \itemize{
#' \item I3: The estimation of (conditional) three-point information without the
#' complexity cost.
#' \item I3k: The estimation of (conditional) three-point information with the
#' complexity cost (\eqn{I3k = I3 - cplx}).
#' \item I2: For reference, the estimation of (conditional) mutual information
#' \eqn{I(X;Y|U)} used in the estimation of \eqn{I3}.
#' \item I2k: For reference, the estimation of regularized (conditional) mutual
#' information \eqn{Ik(X;Y|U)} used in the estimation of \eqn{I3k}.
#' }
#' @export
#' @useDynLib miic
#' @importFrom stats density sd
#'
#' @examples
#' library(miic)
#' N <- 1000
#' # Dependence, conditional independence : X <- Z -> Y
#' Z <- runif(N)
#' X <- Z * 2 + rnorm(N, sd = 0.2)
#' Y <- Z * 2 + rnorm(N, sd = 0.2)
#' res <- computeThreePointInfo(X, Y, Z)
#' message("I(X;Y;Z) = ", res$I3)
#' message("Ik(X;Y;Z) = ", res$I3k)
#'
#' \donttest{
#' # Independence, conditional dependence : X -> Z <- Y
#' X <- runif(N)
#' Y <- runif(N)
#' Z <- X + Y + rnorm(N, sd = 0.1)
#' res <- computeThreePointInfo(X, Y, Z)
#' message("I(X;Y;Z) = ", res$I3)
#' message("Ik(X;Y;Z) = ", res$I3k)
#' }
#'
computeThreePointInfo <- function(X, Y, Z,
                              df_conditioning = NULL,
                              maxbins = NULL,
                              cplx = c("nml", "mdl"),
                              n_eff = -1,
                              sample_weights = NULL,
                              is_continuous = NULL) {
  cplx <- tryCatch(
    {match.arg(cplx)},
    error = function(e) {
      if (grepl("object .* not found", e$message)) {
        message(e, "")
        return("")
      }
      return(toString(cplx))
    }
  )
  cplx <- match.arg(cplx)

  input_data = data.frame(X, Y, Z)
  if (!is.null(df_conditioning)) {
    input_data <- data.frame(input_data, df_conditioning)
  }

  if (!is.null(sample_weights) && length(sample_weights) != nrow(input_data)) {
    stop(paste(
      "Differing number of rows between `sample_weights` and input data:",
      length(sample_weights),
      length(X)
    ))
  }

  complete_row <- stats::complete.cases(input_data)
  n_rows_na <- sum(!complete_row)
  if (n_rows_na > 0) {
    input_data <- input_data[complete_row, ]
    warning(paste0(
      "Removed ", n_rows_na, " rows containing at least one NA value."
    ))
  }

  n_samples <- nrow(input_data)
  n_nodes <- ncol(input_data)

  if (n_samples < 3) {
    stop(paste0("Insufficient number of complete rows: ", nrow(input_data)))
  }

  if (is.null(is_continuous)) {
    is_continuous <- sapply(input_data, is.numeric)
  } else if (length(is_continuous) != n_nodes) {
    stop(paste(
      "Length of `is_continuous` does not match number of input variables:",
      length(is_continuous),
      n_nodes
    ))
  }

  # Numeric factor matrix, level starts from 0
  input_factor <- as.matrix(apply(input_data, 2,
    function(x) (as.numeric(factor(x, levels = unique(x))) - 1))
  )
  max_level_list <- as.numeric(apply(input_factor, 2, max)) + 1
  # Data list, numeric for continuous columns, -1 for discrete columns
  input_double <- matrix(nrow = n_samples, ncol = n_nodes)
  # Order list, order(column) for continuous columns (index starting from 0),
  # -1 for discrete columns
  input_order <- matrix(nrow = n_samples, ncol = n_nodes)
  for (i in c(1: n_nodes)) {
    if (is_continuous[i]) {
      input_double[, i] <- as.numeric(input_data[, i])
      input_order[, i] <- order(input_data[, i], na.last=NA) - 1
    } else {
      input_double[, i] <- rep_len(-1, n_samples)
      input_order[, i] <- rep_len(-1, n_samples)
    }
  }

  arg_list <- list(
    "cplx" = cplx,
    "is_continuous" = is_continuous,
    "levels" = max_level_list,
    "n_eff" = n_eff,
    "n_nodes" = n_nodes,
    "n_samples" = n_samples
  )
  # Continuous variables will be discretized during the computation
  if (any(is_continuous)) {
    initbins <- min(30, round(n_samples**(1 / 3)))
    if (is.null(maxbins) || maxbins > n_samples || maxbins < initbins) {
      maxbins <- min(n_samples, 5 * initbins, 50)
    }
    arg_list[["max_bins"]] <- maxbins
  }
  if (!is.null(sample_weights)) {
    arg_list[["sample_weights"]] <- sample_weights[complete_row, ]
  }
  cpp_input <- list(
    "factor" = as.vector(input_factor),
    "double" = as.vector(input_double),
    "order" = as.vector(input_order)
  )
  # Call cpp code
  rescpp <- miicRGetInfo3Point(cpp_input, arg_list)

  return(rescpp)
}

grid_plot <- function(X, Y, nameDist1, nameDist2) {
  plot_df <- data.frame(table(X, Y), stringAsFactors = TRUE)
  hist2d <- ggplot2::ggplot(plot_df, ggplot2::aes(x=X, y=Y)) +
    ggplot2::geom_tile(
      ggplot2::aes(fill=plot_df$Freq),
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#f4f5fc",
      high = "#0013a3",
      position = "left"
    ) +
    ggplot2::xlab(nameDist1) + ggplot2::ylab(nameDist2) +
    ggplot2::theme_classic()

  g <- ggplot2::ggplot_build(hist2d)
  labels <- g$layout$panel_params[[1]]$y$get_labels()
  labels <- labels[labels != "NA"]

  side_hist_top <- ggplot2::ggplot(data.frame(X), ggplot2::aes(x = X)) +
    ggplot2::geom_bar(color="black", fill="white") +
    theme_side_hist() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(
      5.5, 5.5, -30, 5.5, "pt")) +
    ggplot2::scale_y_continuous(
      labels = labels, # Pass hist2d's labels to align cutpoints on X axis
      breaks = seq(0, 0.1, length.out = length(labels))
    ) +
    ggplot2::ylab("X")

  side_hist_right <- ggplot2::ggplot(data.frame(Y), ggplot2::aes(x = Y)) +
    ggplot2::geom_bar(color="black", fill="white") +
    theme_side_hist() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(
      5.5, 5.5, 5.5, -30, "pt")) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::ylab("Y") +
    ggplot2::coord_flip()

  empty <- ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(1, 1), colour = "white") +
    theme_side_hist()

  return(
    gridExtra::grid.arrange(
      side_hist_top,
      empty,
      hist2d,
      side_hist_right,
      ncol = 2,
      nrow = 2,
      widths = c(4.2, 1),
      heights = c(1, 4.2)
    )
  )
}
