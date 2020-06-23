#' Iterative dynamic programming for (conditional) mutual information through optimized discretization.
#' @description This function chooses cutpoints in the input distributions by maximizing the mutual
#' information minus a complexity cost (computed as BIC or with the Normalized Maximum Likelihood ). The
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
#' \item Verny et al., \emph{PLoS Comp. Bio. 2017.  https://doi.org/10.1371/journal.pcbi.1005662
#' \item Cabeli et al., \emph{PLoS Comp. Bio. 2020.  https://doi.org/10.1371/journal.pcbi.1007866
#' \item Affeldt et al., \emph{Bioinformatics 2016}
#' }
#'
#' @param X [a vector]
#' A vector that contains the observational data of the first variable.
#' @param Y [a vector]
#' A vector that contains the observational data of the second variable.
#' @param matrixU [a numeric matrix]
#' The matrix with the observations of as many columns as conditioning variables.
#' @param maxbins [an int]
#' The maximum number of bins desired in the discretization. A lower number makes the computation faster, a higher
#' number allows finer discretization (by default : 5 * cubic root of N).
#' @param cplx [a string]
#' The complexity used in the dynamic programming. Either "mdl" for Minimum description Length or
#' "nml" for Normalized Maximum Likelihood, which is less costly in the finite sample case and
#' will allow more bins than mdl.
#' @param Neff [an int]
#' The number of effective samples. When there is significant autocorrelation in the samples you may
#' want to specify a number of effective samples that is lower than the number of points in the distribution.
#' @param sample_weights [a vector of floats]
#' Individual weights for each sample, used for the same reason as the effective sample number but with individual
#' precision.
#' @param is_discrete [a vector of booleans]
#' Specify if each variable is to be treated as discrete (TRUE) or continuous (FALSE) in a
#' logical vector of length ncol(matrixU) + 2, in the order [X, Y, U1, U2...]. By default,
#' factors and character vectors are treated as discrete, and numerical vectors as continuous.
#' @param plot [a boolean]
#' Specify if the XY joint space with discretization scheme is to be plotted or not (requires
#' ggplot2 and gridExtra).
#'
#' @return A list that contains :
#' \itemize{
#' \item{two vectors containing the cutpoints for each variable : \emph{cutpoints1}
#' corresponds to /emph{myDist1}, /emph{cutpoints2} corresponds to /emph{myDist2}.}
#' \item{\emph{niterations} is the number of iterations performed before convergence of the (C)MI estimation.}
#' \item{\emph{iterationN}, lists contatining the cutpoint vectors for each iteration.}
#' \item{\emph{info} and \emph{infok}, the estimated (C)MI value and (C)MI minus the complexity cost.}
#' \item{if $emph{plot} == T, a plot object (requires ggplot2 and gridExtra).}
#' @export
#' @useDynLib miic
#'
#' @examples
#' library(miic)
#' \dontrun{
#' N <- 1000
#' # Dependence, conditional independence
#' Z <- runif(N)
#' X <- Z * 2 + rnorm(N, sd = 0.2)
#' Y <- Z * 2 + rnorm(N, sd = 0.2)
#' res <- discretizeMutual(X, Y, plot = T)
#' cat("I(X;Y) = ", res$info)
#' res <- discretizeMutual(X, Y, matrixU = matrix(Z, ncol = 1), plot = T)
#' cat("I(X;Y|Z) = ", res$info)
#'
#' # Conditional independence with categorical conditioning variable
#' Z <- sample(1:3, N, replace = T)
#' X <- -as.numeric(Z == 1) + as.numeric(Z == 2) + 0.2 * rnorm(N)
#' Y <- as.numeric(Z == 1) + as.numeric(Z == 2) + 0.2 * rnorm(N)
#' res <- miic::discretizeMutual(X, Y, cplx = "nml")
#' cat("I(X;Y) = ", res$info)
#' res <- miic::discretizeMutual(X, Y, matrix(Z, ncol = 1), is_discrete = c(F, F, T))
#' cat("I(X;Y|Z) = ", res$info)
#'
#'
#' # Independence, conditional dependence
#' X <- runif(N)
#' Y <- runif(N)
#' Z <- X + Y + runif(N, sd = 0.1)
#' res <- discretizeMutual(X, Y, plot = T)
#' cat("I(X;Y) = ", res$info)
#' res <- discretizeMutual(X, Y, matrixU = matrix(Z, ncol = 1), plot = T)
#' cat("I(X;Y|Z) = ", res$info)
#' }
#'
discretizeMutual <- function(X,
                             Y,
                             matrixU = NULL,
                             maxbins = NULL,
                             cplx = "nml",
                             Neff = NULL,
                             sample_weights = NULL,
                             is_discrete = NULL,
                             plot = T) {
  nameDist1 <- deparse(substitute(X))
  nameDist2 <- deparse(substitute(Y))
  # Check the input arguments
  if (is.null(matrixU)) {
    nbrU <- 0
  } else {
    nbrU <- ncol(matrixU)
  }

  if (is.null(is_discrete)) {
    is_discrete <- c(
      (is.character(X) || is.factor(X)),
      (is.character(Y) || is.factor(Y))
    )
    if (nbrU > 0) {
      for (z in 1:nbrU) {
        is_discrete <- c(is_discrete, (is.character(matrixU[, z]) ||
          is.factor(matrixU[, z])))
      }
    }
  }

  if (all(is_discrete[1:2])) {
    stop("Either X or Y must be continuous to be discretized.")
  }

  if (!(is.vector(X) || is.factor(X)) ||
    !(is.vector(Y) || is.factor(Y))) {
    stop(
      paste0(
        "Please provide the two samples X and Y as numerical vectors ",
        "for continuous variables and factors or character vectors ",
        "for discrete variables."
      )
    )
  }

  if (length(X) != length(Y)) {
    stop(
      paste0(
        "The two samples must have the same number of observation ",
        "(found ",
        length(X),
        " and ",
        length(Y),
        " )."
      )
    )
  }

  if ((!is.null(sample_weights)) &&
    (length(sample_weights) != length(X))) {
    stop(
      paste0(
        "The sample weight vector must be of the same length as the ",
        "number of observations (found ",
        length(sample_weights),
        " while there are ",
        length(X),
        " observations)."
      )
    )
  }

  if ((!is.null(matrixU) && !is.matrix(matrixU)) ||
    (!is.null(matrixU) && nrow(matrixU) != length(X))) {
    stop(
      paste0(
        "matrixU is not a matrix or its number of rows differs from",
        " the number of observations."
      )
    )
  }

  if (!is.null(is_discrete) && (length(is_discrete) != (2 + nbrU))) {
    stop(
      paste0(
        "The vector passed as is_discrete argument must be the same",
        " length as the number of variables, which is ncol(matrixU) ",
        "+ 2."
      )
    )
  }

  initbins <- NULL
  if ((initbins > length(X)) || is.null(initbins)) {
    initbins <- round(length(X)**(1 / 3))
  }

  if ((maxbins > length(X)) ||
    is.null(maxbins) || (maxbins < initbins)) {
    maxbins <- 5 * initbins
  }

  # Remove rows for which any input vector is NA
  NArows <- logical(length(X))
  NArows <- NArows | is.na(X)
  NArows <- NArows | is.na(Y)
  if (!is.null(matrixU)) {
    for (k in 1:ncol(matrixU)) {
      NArows <- NArows | is.na(matrixU[, k])
    }
    matrixU_NA <- matrix(nrow = length(which(!NArows)), ncol = ncol(matrixU))
  }

  if (length(which(NArows)) > 0) {
    warning(paste0(
      "Found ", length(which(NArows)),
      " rows with NAs in at least one of the inputs. Running on ",
      length(which(!NArows)), " samples."
    ))
    X <- X[!NArows]
    Y <- Y[!NArows]
  }

  if (!is.null(matrixU)) {
    for (k in 1:ncol(matrixU)) {
      matrixU_NA[, k] <- matrixU[!NArows, k]
    }
  }

  # Converting factors to discrete numerical variables
  X_orig <- X
  Y_orig <- Y
  if (is_discrete[1]) {
    X <- as.factor(X)
    levels(X) <- 1:nlevels(X)
    X <- as.numeric(X)
  }
  if (is_discrete[2]) {
    Y <- as.factor(Y)
    levels(Y) <- 1:nlevels(Y)
    Y <- as.numeric(Y)
  }
  if (nbrU > 0) {
    for (l in 0:(nbrU - 1)) {
      if (is_discrete[l + 3]) {
        matrixU_NA[, (l + 1)] <- as.factor(matrixU_NA[, (l + 1)])
        levels(matrixU_NA[, (l + 1)]) <- 1:nlevels(matrixU_NA[, (l + 1)])
        matrixU_NA[, (l + 1)] <- as.numeric(matrixU_NA[, (l + 1)])
      }
    }
  }
  cnt_vec <- as.numeric(!is_discrete)

  # Converting matrix to flat vector to pass to cpp
  if (is.null(matrixU)) {
    flatU <- c(0)
  } else {
    class(matrixU_NA) <- "numeric"
    flatU <- as.vector(matrixU_NA)
  }

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

  if (is.null(Neff)) {
    Neff <- length(X)
  }

  if (is.null(sample_weights)) {
    sample_weights <- rep(1, length(X))
  }

  # Number of unique values for each input
  nlevels <- numeric(nbrU + 2)
  for (i in 1:(nbrU + 2)) {
    if (i == 1) {
      nlevels[i] <- length(unique(X))
    } else if (i == 2) {
      nlevels[i] <- length(unique(Y))
    } else {
      nlevels[i] <- length(unique(matrixU_NA[, i - 2]))
    }
  }

  # Call cpp code
  if (base::requireNamespace("Rcpp", quietly = TRUE)) {
    rescpp <- mydiscretizeMutual(
        X,
        Y,
        flatU,
        nbrU,
        maxbins,
        initbins,
        intcplx,
        cnt_vec,
        nlevels,
        Neff,
        sample_weights
      )
  }
  niterations <- nrow(rescpp$cutpointsmatrix) / maxbins

  result <- list()
  epsilon <- min(c(sd(X), sd(Y))) / 100
  result$niterations <- niterations
  for (i in 0:(niterations - 1)) {
    result[[paste0("iteration", i + 1)]] <- list()
    for (l in 1:2) {
      if (!is_discrete[l]) {
        clean_cutpoints <- rescpp$cutpointsmatrix[, l][(maxbins * i) +
          (1:maxbins)]
        clean_cutpoints <- clean_cutpoints[clean_cutpoints != -1]
        if (l == 1) {
          data <- X
        } else {
          if (l == 2) {
            data <- Y
          } else {
            data <- matrixU[, l - 2]
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
  result$efinfo <- rescpp$efinfo

  if (plot) {
    if (require(ggplot2, quietly = TRUE) & require(gridExtra, quietly = TRUE)) {
      if (!any(is_discrete[1:2])) {
        result$plot <- jointplot_hist(X, Y, result, nameDist1, nameDist2)
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
  theme_classic() %+replace%
    theme(
      title = element_text(
        family = "",
        face = "plain",
        color = NA,
        size = theme_classic()$text$size
      ),
      axis.text = element_text(
        color = NA,
        size = theme_classic()$axis.text$size
      ),
      axis.line.x = element_line(colour = NA),
      axis.line.y = element_line(colour = NA),
      axis.ticks = element_blank(),
      panel.background = element_rect(fill = "white", colour = "white")
    )
}


jointplot_hist <- function(X, Y, result, nameDist1, nameDist2,
                           title = "Joint histogram") {
  cut_points1 <- result$cutpoints1
  cut_points2 <- result$cutpoints2

  # Custom density matrix for colour filling relative to 2D bin area in
  # the 2d histogram
  bin_count <- table(
    cut(X, cut_points1, include.lowest = T),
    cut(Y, cut_points2, include.lowest = T)
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

  hist2d <- ggplot(fill_density_flat) +
    geom_rect(
      aes(
        xmin = xstart,
        xmax = xend,
        ymin = ystart,
        ymax = yend,
        fill = log(density)
      ),
      na.rm = T,
      show.legend = F
    ) +
    scale_fill_gradient(
      low = "#f4f5fc",
      high = "#0013a3",
      position = "left",
      na.value = "white",
      limits = log(c(
        min_density,
        max(fill_density_flat$density,
          na.rm = T
        )
      ))
    ) +
    geom_vline(
      xintercept = cut_points1,
      linetype = "dashed",
      color = "grey"
    ) +
    geom_hline(
      yintercept = cut_points2,
      linetype = "dashed",
      color = "grey"
    ) +
    geom_point(
      data = data.frame(X, Y),
      aes(x = X, y = Y),
      shape = 21,
      alpha = .7,
      fill = "#ffef77",
      size = 2
    ) +
    xlab(nameDist1) + ylab(nameDist2) +
    theme_classic()

  g <- ggplot_build(hist2d)
  labels <- g$layout$panel_params[[1]]$y.labels


  side_hist_top <- ggplot(data.frame(X), aes(x = X)) +
    # Histogram with density instead of count on y-axis
    geom_histogram(
      aes(y = ..density..),
      breaks = cut_points1,
      colour = "black",
      fill = "white"
    ) +
    # Overlay with transparent density plot
    geom_density(
      adjust = 0.5,
      alpha = .5,
      fill = "#c1c6ee"
    ) +
    theme_side_hist() %+replace% theme(plot.margin = margin(
      5.5, 5.5, -25, 5.5,
      "pt"
    )) +
    # , expand = c(0, 0)) + #Pass hist2d's labels to align both X axes
    scale_y_continuous(
      labels = labels,
      breaks = seq(0, 0.1, length.out = length(labels))
    ) +
    ylab("X")
  # Histogram with density instead of count on y-axis
  side_hist_bot <- ggplot(data.frame(Y), aes(x = Y)) +
    geom_histogram(
      aes(y = ..density..),
      breaks = cut_points2,
      colour = "black",
      fill = "white"
    ) +
    # Overlay with transparent density plot
    geom_density(
      adjust = 0.5,
      alpha = .5,
      fill = "#c1c6ee"
    ) +
    theme_side_hist() %+replace% theme(plot.margin = margin(5.5, 5.5, 5.5, -30, "pt")) +
    scale_y_continuous(expand = c(0, 0)) +
    ylab("Y") +
    coord_flip()

  I2 <- result$info
  I2p <- result$infok

  empty <- ggplot() + geom_point(aes(1, 1), colour = "white") +
    geom_text(aes(
      x = 1,
      y = 0.5,
      label = paste0("I = ", round(I2, 3))
    )) +
    geom_text(aes(
      x = 1,
      y = 0,
      label = paste0("Ik = ", round(I2p, 3))
    )) +
    theme(
      axis.ticks = element_blank(),
      panel.background = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank()
    )


  return(
    gridExtra::grid.arrange(
      side_hist_top,
      empty,
      hist2d,
      side_hist_bot,
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

  barplot <- ggplot(plot_df, aes(x = X, y = Y)) +
    geom_boxplot(outlier.shape = NA, fill = "#c1c6ee") +
    geom_jitter(width = 0.2, height = 0) +
    geom_hline(
      yintercept = cut_points,
      linetype = "dashed",
      color = "grey",
      size = 1
    ) +
    theme_classic() +
    xlab(Xname) +
    ylab(Yname)

  return(barplot)
}
