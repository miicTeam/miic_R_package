#' Iterative dynamic programming for (conditional) mutual information through optimized discretization.
#' @description This function chooses cutpoints in the input distributions by maximizing the mutual
#' information minus a complexity cost (computed as BIC or with the Normalized Maximum Likelihood). The
#' (conditional) mutual information computed on the optimized discretized distributions effectively approaches
#' the mutual information computed on the original continuous variables.
#'
#' @details For a pair of variables \eqn{X} and \eqn{Y}, the algorithm will in turn choose cutpoints on \eqn{X}
#' then on \eqn{Y}, maximizing \eqn{I(X_{d};Y_{d}) - cplx(X_{d};Y_{d})} where \eqn{cplx(X_{d};Y_{d})} is the
#' complexity cost of the considered discretizations of \eqn{X} and \eqn{Y} (see Affeldt 2016 and Cabeli 2020).
#' When the value \eqn{I(X_{d};Y_{d})} is stable between two iterations the discretization scheme of
#' \eqn{X_{d}} and \eqn{Y_{d}} is returned as well as \eqn{I(X_{d};Y_{d})} and \eqn{I(X_{d};Y_{d})-cplx(X_{d};Y_{d})}.
#'
#' With a set of conditioning variables \eqn{U}, the discretization scheme maximizes each term of the sum
#' \eqn{I(X;Y|U) \sim 0.5*(I(X_{d};Y_{d}, U_{d}) - I(X_{d};U_{d}) + I(Y_{d};X_{d}, U_{d}) - I(Y_{d};U_{d}))}.
#'
#' Discrete variables can be passed as factors and will be used "as is" to maximize each term.
#'
#'
#' @references
#' \itemize{
#' \item Verny et al., \emph{PLoS Comp. Bio. 2017.}  https://doi.org/10.1371/journal.pcbi.1005662
#' \item Cabeli et al., \emph{PLoS Comp. Bio. 2020.}  https://doi.org/10.1371/journal.pcbi.1007866
#' \item Affeldt et al., \emph{Bioinformatics 2016}
#' }
#'
#' @param x [a vector]
#' The \eqn{X} vector that contains the observational data of the first variable.
#' @param y [a vector]
#' The \eqn{Y} vector that contains the observational data of the second variable.
#' @param matrix_u [a numeric matrix]
#' The matrix with the observations of as many columns as conditioning variables.
#' @param maxbins [an int]
#' The maximum number of bins desired in the discretization. A lower number makes the computation faster, a higher
#' number allows finer discretization (by default : 5 * cubic root of N).
#' @param cplx [a string]
#' The complexity used in the dynamic programming. Either "mdl" for Minimum description Length or
#' "nml" for Normalized Maximum Likelihood, which is less costly in the finite sample case and
#' will allow more bins than mdl.
#' @param n_eff [an int]
#' The number of effective samples. When there is significant autocorrelation in the samples you may
#' want to specify a number of effective samples that is lower than the number of points in the distribution.
#' @param sample_weights [a vector of floats]
#' Individual weights for each sample, used for the same reason as the effective sample number but with individual
#' precision.
#' @param is_discrete [a vector of booleans]
#' Specify if each variable is to be treated as discrete (TRUE) or continuous (FALSE) in a
#' logical vector of length ncol(matrix_u) + 2, in the order [X, Y, U1, U2...]. By default,
#' factors and character vectors are treated as discrete, and numerical vectors as continuous.
#' @param plot [a boolean]
#' Specify if the XY joint space with discretization scheme is to be plotted or not (requires
#' ggplot2 and gridExtra).
#'
#' @return A list that contains :
#' \itemize{
#' \item{two vectors containing the cutpoints for each variable :
#' \emph{cutpoints1} corresponds to \emph{x},
#' \emph{cutpoints2} corresponds to \emph{y}.}
#' \item{\emph{niterations} is the number of iterations performed before convergence of the (C)MI estimation.}
#' \item{\emph{iteration1, iteration2, ...}, lists containing the cutpoint vectors for each iteration.}
#' \item{\emph{info} and \emph{infok}, the estimated (C)MI value and (C)MI minus the complexity cost.}
#' \item{if \emph{plot} == TRUE, a plot object (requires ggplot2 and gridExtra).}
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
#' res <- discretizeMutual(X, Y, plot = FALSE)
#' message("I(X;Y) = ", res$info)
#' res <- discretizeMutual(X, Y, matrix_u = matrix(Z, ncol = 1), plot = FALSE)
#' message("I(X;Y|Z) = ", res$info)
#'
#' \donttest{
#' # Conditional independence with categorical conditioning variable : X <- Z -> Y
#' Z <- sample(1:3, N, replace = TRUE)
#' X <- -as.numeric(Z == 1) + as.numeric(Z == 2) + 0.2 * rnorm(N)
#' Y <- as.numeric(Z == 1) + as.numeric(Z == 2) + 0.2 * rnorm(N)
#' res <- miic::discretizeMutual(X, Y, cplx = "nml")
#' message("I(X;Y) = ", res$info)
#' res <- miic::discretizeMutual(X, Y, matrix(Z, ncol = 1), is_discrete = c(FALSE, FALSE, TRUE))
#' message("I(X;Y|Z) = ", res$info)
#'
#'
#' # Independence, conditional dependence : X -> Z <- Y
#' X <- runif(N)
#' Y <- runif(N)
#' Z <- X + Y + rnorm(N, sd = 0.1)
#' res <- discretizeMutual(X, Y, plot = TRUE)
#' message("I(X;Y) = ", res$info)
#' res <- discretizeMutual(X, Y, matrix_u = matrix(Z, ncol = 1), plot = TRUE)
#' message("I(X;Y|Z) = ", res$info)
#' }
#'
discretizeMutual <- function(x,
                             y,
                             matrix_u = NULL,
                             maxbins = NULL,
                             cplx = "nml",
                             n_eff = NULL,
                             sample_weights = NULL,
                             is_discrete = NULL,
                             plot = TRUE) {
  nameDist1 <- deparse(substitute(x))
  nameDist2 <- deparse(substitute(y))
  # Check the input arguments
  if (is.null(matrix_u)) {
    nbrU <- 0
  } else {
    nbrU <- ncol(matrix_u)
  }

  if (is.null(is_discrete)) {
    is_discrete <- c(
      (is.character(x) || is.factor(x)),
      (is.character(y) || is.factor(y))
    )
    if (nbrU > 0) {
      for (z in 1:nbrU) {
        is_discrete <- c(is_discrete, (is.character(matrix_u[, z]) ||
          is.factor(matrix_u[, z])))
      }
    }
  }

  if (all(is_discrete[1:2])) {
    stop("Either x or y must be continuous to be discretized.")
  }

  if (!(is.vector(x) || is.factor(x)) ||
    !(is.vector(y) || is.factor(y))) {
    stop(
      paste0(
        "Please provide the two samples x and y as numerical vectors ",
        "for continuous variables and factors or character vectors ",
        "for discrete variables."
      )
    )
  }

  if (length(x) != length(y)) {
    stop(
      paste0(
        "The two samples must have the same number of observation ",
        "(found ",
        length(x),
        " and ",
        length(y),
        " )."
      )
    )
  }

  if ((!is.null(sample_weights)) &&
    (length(sample_weights) != length(x))) {
    stop(
      paste0(
        "The sample weight vector must be of the same length as the ",
        "number of observations (found ",
        length(sample_weights),
        " while there are ",
        length(x),
        " observations)."
      )
    )
  }

  if ((!is.null(matrix_u) && !is.matrix(matrix_u)) ||
    (!is.null(matrix_u) && nrow(matrix_u) != length(x))) {
    stop(
      paste0(
        "matrix_u is not a matrix or its number of rows differs from",
        " the number of observations."
      )
    )
  }

  if (!is.null(is_discrete) && (length(is_discrete) != (2 + nbrU))) {
    stop(
      paste0(
        "The vector passed as is_discrete argument must be the same",
        " length as the number of variables, which is ncol(matrix_u) ",
        "+ 2."
      )
    )
  }

  # Remove rows for which any input vector is NA
  matrix_u_NA <- matrix()
  NArows <- logical(length(x))
  NArows <- NArows | is.na(x)
  NArows <- NArows | is.na(y)
  if (!is.null(matrix_u)) {
    for (k in 1:ncol(matrix_u)) {
      NArows <- NArows | is.na(matrix_u[, k])
    }
    matrix_u_NA <- matrix(nrow = length(which(!NArows)), ncol = ncol(matrix_u))
  }

  if (length(which(NArows)) > 0) {
    warning(paste0(
      "Found ", length(which(NArows)),
      " rows with NAs in at least one of the inputs. Running on ",
      length(which(!NArows)), " samples."
    ))
    x <- x[!NArows]
    y <- y[!NArows]
  }
  if (length(which(!NArows)) < 3) {
    stop(paste0(
      "There are only ", length(which(!NArows)),
      " rows with complete data, not enough non-NA samples."
    ))
  }

  if (!is.null(matrix_u)) {
    for (k in 1:ncol(matrix_u)) {
      matrix_u_NA[, k] <- matrix_u[!NArows, k]
    }
  }

  initbins <- min(30, round(length(x)**(1 / 3)))

  if (is.null(maxbins) || maxbins > length(x) || maxbins < initbins) {
    maxbins <- min(length(x), 5 * initbins, 50)
  }

  # Converting factors to discrete numerical variables
  X_orig <- x
  Y_orig <- y
  if (is_discrete[1]) {
    x <- as.factor(x)
    levels(x) <- 1:nlevels(x)
    x <- as.numeric(x)
  }
  if (is_discrete[2]) {
    y <- as.factor(y)
    levels(y) <- 1:nlevels(y)
    y <- as.numeric(y)
  }
  if (nbrU > 0) {
    for (l in 0:(nbrU - 1)) {
      if (is_discrete[l + 3]) {
        matrix_u_NA[, (l + 1)] <- as.factor(matrix_u_NA[, (l + 1)])
        levels(matrix_u_NA[, (l + 1)]) <- 1:nlevels(matrix_u_NA[, (l + 1)])
        matrix_u_NA[, (l + 1)] <- as.numeric(matrix_u_NA[, (l + 1)])
      }
    }
  }
  is_continuous <- !is_discrete

  # Pass complexity parameter as int
  if (cplx == "mdl") {
    intcplx <- 0
  } else if (cplx == "nml") {
    intcplx <- 1
  } else {
    warning(
      paste0(
        "cplx parameter not understood, please specify either ",
        "\'mdl\' or \'nml\'. Running with the default option ",
        "(nml)."
      )
    )
    intcplx <- 1
  }

  if (is.null(n_eff)) {
    n_eff <- length(x)
  }

  if (is.null(sample_weights)) {
    sample_weights <- numeric(0);
  }

  input_data = data.frame(x,y)
  if(!all(is.na(matrix_u_NA))) input_data = cbind(input_data, matrix_u_NA)
  n_samples <- nrow(input_data)
  n_nodes <- ncol(input_data)
  input_factor <- apply(input_data, 2, function(x)
                        (as.numeric(factor(x, levels = unique(x))) - 1))
  input_factor[is.na(input_factor)] <- -1
  max_level_list <- as.numeric(apply(input_factor, 2, max)) + 1
  input_factor <- as.vector(as.matrix(input_factor))
  # Order list, order(column) for continuous columns (index starting from 0, NA
  # mapped to -1), -1 for discrete columns
  input_order <- matrix(nrow = n_samples, ncol = n_nodes)
  for (i in c(1: n_nodes)) {
    if (is_continuous[i]) {
      n_NAs <- sum(is.na(input_data[, i]))
      input_order[, i] <- c(order(input_data[, i], na.last=NA) - 1,
                            rep_len(-1, n_NAs))
    } else {
      input_order[, i] <- rep_len(-1, n_samples)
    }
  }
  input_order <- as.vector(input_order)

  var_names <- colnames(input_data)

  arg_list <- list(
    "cplx" = cplx,
    "is_continuous" = is_continuous,
    "levels" = max_level_list,
    "max_bins" = maxbins,
    "n_eff" = n_eff,
    "n_nodes" = n_nodes,
    "n_samples" = n_samples,
    "sample_weights" = sample_weights
  )
  cpp_input <- list("factor" = input_factor, "order" = input_order)
  # Call cpp code
  if (base::requireNamespace("Rcpp", quietly = TRUE)) {
    rescpp <- mydiscretizeMutual(cpp_input, arg_list)
  }
  niterations <- nrow(rescpp$cutpointsmatrix) / maxbins

  result <- list()
  epsilon <- min(c(sd(x), sd(y))) / 100
  result$niterations <- niterations
  for (i in 0:(niterations - 1)) {
    result[[paste0("iteration", i + 1)]] <- list()
    for (l in 1:2) {
      if (!is_discrete[l]) {
        clean_cutpoints <- rescpp$cutpointsmatrix[, l][(maxbins * i) +
          (1:maxbins)]
        clean_cutpoints <- clean_cutpoints[clean_cutpoints != -1]
        if (l == 1) {
          data <- x
        } else {
          if (l == 2) {
            data <- y
          } else {
            data <- matrix_u[, l - 2]
          }
        }
        uniquedata <- sort(unique(data))

        clean_cutpoints <- sort(data)[clean_cutpoints + 1]
        if (length(clean_cutpoints) > 0) {
          # Take midpoints between two consecutive unique values instead of
          # the values themselves
          clean_cutpoints <- sapply(clean_cutpoints, function(x) {
            if (x < uniquedata[length(uniquedata)]) {
              ((min(uniquedata[uniquedata > x]) +
                max(uniquedata[uniquedata <= x])) / 2)
            } else {
              x
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
  }
  for (l in 1:(nbrU + 2)) {
    result[[paste0("cutpoints", l)]] <-
      result[[paste0("iteration", niterations)]][[paste0("cutpoints", l)]]
  }

  result$info <- rescpp$info
  result$infok <- rescpp$infok

  if (plot) {
    if (base::requireNamespace("ggplot2", quietly = TRUE) & base::requireNamespace("gridExtra", quietly = TRUE)) {
      if (!any(is_discrete[1:2])) {
        result$plot <- jointplot_hist(x, y, result, nameDist1, nameDist2)
      } else if (!all(is_discrete[1:2])) {
        result$plot <- barplot_disc(
          X_orig,
          Y_orig,
          result,
          is_discrete,
          nameDist1,
          nameDist2
        )
      }
    } else {
      warning("Plotting requires ggplot2 and gridExtra.")
    }
  }

  result
}

# Plot functions

axisprint <- function(x) {
  sprintf("%6s", x)
}

theme_side_hist <- function() {
  ggplot2::theme_classic() +
    ggplot2::theme(
      title = ggplot2::element_text(
        family = "",
        face = "plain",
        color = NA,
        size = ggplot2::theme_classic()$text$size
      ),
      axis.text = ggplot2::element_text(
        color = NA,
        size = ggplot2::theme_classic()$axis.text$size
      ),
      axis.line.x = ggplot2::element_line(colour = NA),
      axis.line.y = ggplot2::element_line(colour = NA),
      axis.ticks = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = "white")
    )
}


jointplot_hist <- function(X, Y, result, nameDist1, nameDist2,
                           title = "Joint histogram") {
  cut_points1 <- result$cutpoints1
  cut_points2 <- result$cutpoints2

  # Custom density matrix for colour filling relative to 2D bin area in
  # the 2d histogram
  bin_count <- table(
    cut(X, cut_points1, include.lowest = TRUE),
    cut(Y, cut_points2, include.lowest = TRUE)
  )
  bin_areas <-
    (cut_points1[-1] - cut_points1[1:(length(cut_points1) - 1)]) %*%
    t(cut_points2[-1] - cut_points2[1:(length(cut_points2) - 1)])
  fill_density <- bin_count / bin_areas
  min_density <- 1 / max(bin_areas[bin_count > 0]) / sum(fill_density)
  fill_density <- fill_density / sum(fill_density)
  fill_density_flat <- data.frame(
    xstart = numeric(),
    xend = numeric(),
    ystart = numeric(),
    yend = numeric(),
    density = numeric()
  )
  for (j in 1:(ncol(fill_density))) {
    for (i in 1:(nrow(fill_density))) {
      fill_density_flat[(j - 1) * nrow(fill_density) + i, ] <-
        c(
          cut_points1[i],
          cut_points1[i + 1],
          cut_points2[j],
          cut_points2[j + 1],
          fill_density[i, j]
        )
    }
  }
  fill_density_flat[fill_density_flat$density == 0, "density"] <- NA
  fill_density_flat$logdensity = log(fill_density_flat$density)

  hist2d <- ggplot2::ggplot(fill_density_flat) +
    ggplot2::geom_rect(
      ggplot2::aes_string(
        xmin = "xstart",
        xmax = "xend",
        ymin = "ystart",
        ymax = "yend",
        fill = "logdensity"
      ),
      na.rm = TRUE,
      show.legend = FALSE
    ) +
    ggplot2::scale_fill_gradient(
      low = "#f4f5fc",
      high = "#0013a3",
      position = "left",
      na.value = "white",
      limits = log(c(
        min_density,
        max(fill_density_flat$density,na.rm = TRUE))
      )
    ) +
    ggplot2::geom_vline(
      xintercept = cut_points1,
      linetype = "dashed",
      color = "grey"
    ) +
    ggplot2::geom_hline(
      yintercept = cut_points2,
      linetype = "dashed",
      color = "grey"
    ) +
    ggplot2::geom_point(
      data = data.frame(X, Y),
      ggplot2::aes(x = X, y = Y),
      shape = 21,
      alpha = .7,
      fill = "#ffef77",
      size = 2
    ) +
    ggplot2::xlab(nameDist1) + ggplot2::ylab(nameDist2) +
    ggplot2::theme_classic()

  g <- ggplot2::ggplot_build(hist2d)
  labels <- g$layout$panel_params[[1]]$y$get_labels()
  labels <- labels[labels != "NA"]

  side_hist_top <- ggplot2::ggplot(data.frame(X), ggplot2::aes(x = X)) +
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      breaks = cut_points1,
      colour = "black",
      fill = "white"
    ) +
    # Overlay with transparent density plot
    ggplot2::geom_density(
      adjust = 0.5,
      alpha = .5,
      fill = "#c1c6ee"
    ) +
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
    ggplot2::geom_histogram(
      ggplot2::aes(y = ggplot2::after_stat(density)),
      breaks = cut_points2,
      colour = "black",
      fill = "white"
    ) +
    # Overlay with transparent density plot
    ggplot2::geom_density(
      adjust = 0.5,
      alpha = .5,
      fill = "#c1c6ee"
    ) +
    theme_side_hist() +
    ggplot2::theme(
      plot.margin = ggplot2::margin(
      5.5, 5.5, 5.5, -30, "pt")) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::ylab("Y") +
    ggplot2::coord_flip()

  I2 <- result$info
  I2p <- result$infok

  empty <- ggplot2::ggplot() + ggplot2::geom_point(ggplot2::aes(1, 1), colour = "white") +
    ggplot2::geom_text(ggplot2::aes(
      x = 1,
      y = 0.5,
      label = paste0("I = ", round(I2, 3))
    )) +
    ggplot2::geom_text(ggplot2::aes(
      x = 1,
      y = 0,
      label = paste0("Ik = ", round(I2p, 3))
    )) +
    ggplot2::theme(
      axis.ticks = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank()
    )


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
  # , bottom = title))
}

barplot_disc <- function(X,
                         Y,
                         result,
                         is_discrete,
                         nameDist1,
                         nameDist2,
                         title = "Joint histogram") {
  cut_points <- result$cutpoints2
  Xname <- nameDist1
  Yname <- nameDist2
  if (is_discrete[2]) {
    # Swap X and Y: by convention, X is the discrete variable
    X_ <- X
    X <- Y
    Y <- X_
    Xname <- nameDist2
    Yname <- nameDist1
    cut_points <- result$cutpoints1
  }

  plot_df <- data.frame(
    X = X,
    Y = Y,
    stringAsFactors = TRUE
  )

  barplot <- ggplot2::ggplot(plot_df, ggplot2::aes(x = X, y = Y)) +
    ggplot2::geom_boxplot(outlier.shape = NA, fill = "#c1c6ee") +
    ggplot2::geom_jitter(width = 0.2, height = 0) +
    ggplot2::geom_hline(
      yintercept = cut_points,
      linetype = "dashed",
      color = "grey",
      size = 1
    ) +
    ggplot2::theme_classic() +
    ggplot2::xlab(Xname) +
    ggplot2::ylab(Yname)

  return(barplot)
}
