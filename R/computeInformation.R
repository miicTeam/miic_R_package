#*******************************************************************************
# Filename   : computeInformation.R
#
# Description: Compute 2 and 3 point (conditional) mutual information
#*******************************************************************************

#===============================================================================
# FUNCTIONS (internal)
#===============================================================================
# compute_mi_batch
#-------------------------------------------------------------------------------
# Compute mutual information (MI) between a set of variables and variable(s)
# of interest.
#
# @param input_data [a data frame or a matrix, required]
#
# Expected layout is samples as rows and variables as columns. Column names
# must contain the names of the variables.
#
# @param var_of_interest_names [a string or vector of strings, optional,
# NULL by default]
#
# For the variable(s) of interest that are part of the \emph{input_data},
# you should supply their names here.
#
# @param var_of_interest_values [a data frame, optional, NULL by default]
#
# For the variables of interest that are not in \emph{input_data},
# a data frame can be supplied. The column names are the names
# of the variables of interest and rows are the samples
# ordered in the same way as the \emph{input_data}.
# Typically, such variables are metadata associated to samples but not stored
# in \emph{input_data}, e.g. a "Treatment" vs "Control" variable
# in an experiment and a count matrix with the expression of genes
# in \emph{input_data}.
#
# @param skip_cheks [a boolean, optional, FALSE by default]
#
# Before computing MI between the variable of interest and the features,
# \emph{input_data} is checked to filter out constant features and rows full
# of NAs. When the \emph{input_data} does not need such filtering,
# these checks can be skipped to speed up the process.
#
# @param precomputed_mis [a matrix, optional, NULL by default]
#
# if MIs has been previously computed between some variables
# in the \emph{input_data} and variable(s) of interest,
# supplying these precomputed MIs will speed up the process as the existing
# MIs (values different from NA) will not be recomputed.
# This matrix must have variables names from the \emph{input_data}
# as row names and variables of interest names as column names
# (the layout is the same as the \emph{mis} matrix returned).
#
# @param n_threads [a positive integer, optional, 1 by default]
#
# When set greater than 1, \emph{n_threads} parallel threads will be used for
# computation. Make sure your compiler is compatible with openmp
# if you wish to use multithreading.
#
# @param verbose [an integer, optional, 3 by default]
#
# Level of verbosity: 0=no display, 1=summary, 2=progress per variable of
# interest, 3=same as 2 with display of estimated time remaining.
#
# @return A matrix with the MIs (in bits) between the variables in
# \emph{input_data} (as rows) and the variable(s) of interest (as columns).
# Row and column names are sorted alphabetically.\cr
# If a precomputed MIs matrix was supplied, new values computed
# are added to the existing matrix.
#-------------------------------------------------------------------------------
compute_mi_batch <- function (input_data,
  var_of_interest_names=NULL, var_of_interest_values=NULL,
  skip_cheks=F, precomputed_mis=NULL, n_threads=1, verbose=3)
  {
  # LN_2 equivalent to constant for the function
  #
  LN_2 <- log(2)
  #
  # Check input_data
  #
  if (  ( ! is.data.frame (input_data) )
     && ( ! is.matrix(input_data) )
     && ( ! inherits(input_data, "Matrix") ) )
    miic:::miic_error  ("parameters",
      "the input data must be a data frame or a matrix.")
  if ( (ncol (input_data) <= 0) || (nrow (input_data) <= 0) )
    miic:::miic_error  ("parameters", "the input data is empty.")
  if ( is.data.frame (input_data) )
    # Ensure we have a true data frame, e.g. not a tibble
    # TODO evaluate run time impact on very large data frames
    # (let matrices unchanged to avoid warnings on large memory allocation)
    input_data <- as.data.frame (input_data)
  if ( is.null (colnames (input_data) ) )
    miic:::miic_error  ("parameters", "the input data must have column names.")
  if ( length(unique(colnames (input_data))) != ncol(input_data) )
    miic:::miic_error  ("parameters", "the input data have some column names duplicated.")
  #
  # Check variables of interest and group their values into a data frame
  #
  if ( is.null (var_of_interest_names) && is.null (var_of_interest_values) )
    miic:::miic_error  ("parameters", "the name of the variable(s) of interest",
      " or a data frame with the variable(s) of interest values must be supplied.")

  if (is.null (var_of_interest_names) )
    var_of_interest_names <- c()
  else
    {
    for (one_var_name in var_of_interest_names)
      if ( miic:::test_param_wrong_string (one_var_name, colnames(input_data) ) )
        miic:::miic_error ("parameters",  "Some of the variable of interest",
                           " names are incorrect or not in the input_data.")
    }

  extra_voi_names <- c()
  if ( ! is.null (var_of_interest_values) )
    {
    if ( ! is.data.frame (var_of_interest_values) )
      miic:::miic_error  ("parameters",
        "the var_of_interest_values must be a data frame.")
    # Ensure we have a true data frame, i.e. not a tibble
    var_of_interest_values <- as.data.frame (var_of_interest_values)
    if (ncol (var_of_interest_values) <= 0)
      {
      miic:::miic_warning  ("parameters",
        "the var_of_interest_values data frame has been supplied but is empty.")
      var_of_interest_values <- NULL
      }
    else if ( nrow (var_of_interest_values) != nrow (input_data) )
      miic:::miic_error  ("parameters",
        "the variable of interest values does not match the number of samples.")
    else
      {
      # Data frame OK, checks variables names not in data
      #
      extra_voi_names <- colnames (var_of_interest_values)
      #
      # Error or warning if voi requested as external in var_of_interest_values
      # are present in the input data
      #
      poss_wrong_idx = which ( extra_voi_names %in% colnames(input_data) )
      if (length (poss_wrong_idx) >= 1)
        {
        for (one_var_name in extra_voi_names[poss_wrong_idx])
          {
          one_var_voi_vals = var_of_interest_values[,one_var_name]
          one_var_input_vals = input_data[,one_var_name]

          if (any ( ( is.na (one_var_voi_vals) != is.na (one_var_input_vals) )
                  | (one_var_voi_vals[ !is.na(one_var_voi_vals) ] != one_var_input_vals[ !is.na(one_var_input_vals) ]) ) )
            miic:::miic_error  ("parameters",
              "the variable ", one_var_name, " is supplied both in input data",
              " and in variables of interest values.")
          #
          # Supplied in both in input data and in variables of interest values
          # and with identical values => just a warning, use input_data
          # and ignore variables of interest values
          #
          miic:::miic_warning  ("parameters",
            "the variable ", one_var_name, " is supplied both in input data",
            " and in variables of interest values.")
          var_of_interest_names = unique (c (var_of_interest_names, one_var_name) )
          var_of_interest_values[,one_var_name] = NULL
          }
        extra_voi_names <- colnames (var_of_interest_values)
        }
      }
    }
  all_voi_names <- unique ( c (var_of_interest_names, extra_voi_names) )
  #
  # Pre-computed MIs checks and preparation of returned value mat_mis
  #
  mat_mis = precomputed_mis
  if ( ! is.null (mat_mis) )
    {
    if (  ( ! is.matrix(mat_mis) )
       && ( ! inherits(mat_mis, "Matrix") ) )
      miic:::miic_error  ("parameters",
        "the precomputed MIs must be a matrix.")
    # Ensure we have a true data frame. i.e. not a tibble
    if ( is.null (colnames (mat_mis) ) )
      miic:::miic_error  ("parameters",
        "the precomputed MIs matrix has no column name.")
    if ( is.null (rownames (mat_mis) ) )
      miic:::miic_error  ("parameters",
        "the precomputed MIs matrix has no row name.")
    if ( (nrow (mat_mis) == 0) || (ncol (mat_mis) == 0) )
      mat_mis <- NULL
    }

  flag_mi_precomp = F
  if ( is.null (mat_mis) )
    mat_mis = matrix (NA_real_,
                     nrow = ncol (input_data),
                     ncol = length (all_voi_names),
                     dimnames = list ( sort (colnames (input_data)),
                                       sort (all_voi_names) ) )
  else
    {
    flag_mi_precomp = T
    #
    # Add missing row / columns
    #
    missing_row_names = colnames (input_data)[
      ! ( colnames (input_data) %in% rownames (mat_mis) ) ]
    if (length (missing_row_names) > 0)
      {
      n_rows_mis = nrow (mat_mis)
      range_to_add = (n_rows_mi+1):(n_rows_mis+length(missing_row_names))
      mat_mis[range_to_add, ] = NA_real_
      rownames (mat_mis)[range_to_add] = missing_row_names
      mat_mis = mat_mis[order( rownames (mat_mis) ), , drop=F]
      }
    missing_col_names = all_voi_names[
      ! ( all_voi_names %in% colnames (mat_mis) ) ]
    if (length (missing_col_names) > 0)
      {
      n_cols_mis = ncol (mat_mis)
      range_to_add = (n_cols_mis+1):(n_cols_mis+length(missing_col_names))
      mat_mis[, range_to_add] = NA_real_
      colnames (mat_mis)[range_to_add] = missing_col_names
      mat_mis = mat_mis[, order( colnames (mat_mis) ), drop=F]
      }
    }
  # print (mat_mis[1:5,1:4])
  #
  # Other parameters checks
  #
  skip_cheks <- miic:::check_param_logical (skip_cheks, "skip checks", F)
  n_threads <- miic:::check_param_int (n_threads, "number of threads", 1)
  verbose <- miic:::check_param_int (verbose, "verbose", 3, 0, 3)
  #
  # The bin size controls the number of features evaluated in one go
  #
  if ( (ncol (input_data) < 750) || (nrow(input_data) <= 2000) )
    bin_size = 100
  else if (nrow(input_data) <= 4000)
    bin_size = 75
  else if (nrow(input_data) <= 7000)
    bin_size = 50
  else
    bin_size = 20
  #
  # For each variable of interest (voi), compute MI
  #
  n_all_vois = length (all_voi_names)
  max_voi_name = max (unlist (lapply (all_voi_names, FUN=nchar) ) )
  n_vars = ncol (input_data)
  time_start = Sys.time()
  data_for_compute = input_data
  one_voi_idx = 3
  for (one_voi_idx in 1:n_all_vois)
    {
    one_voi_name = all_voi_names[[one_voi_idx]]
    str_progress_start = paste0 ("Computing MI for ", one_voi_name,
      paste (rep ( ' ', max_voi_name - nchar(one_voi_name) ), collapse="" ), " : ")
    if (verbose >= 2)
      cat (paste0 (str_progress_start, "0 %\r") )

    if (flag_mi_precomp)
      {
      var_to_recomp = rownames(mat_mis) [ is.na (mat_mis[, one_voi_name]) ]
      one_voi_name_in_recomp_idx = which(var_to_recomp == one_voi_name)
      if (length (one_voi_name_in_recomp_idx) > 0)
        var_to_recomp = var_to_recomp[ -one_voi_name_in_recomp_idx ]
      if (length (var_to_recomp) == 0)
        {
        if (verbose >= 2)
          cat (paste0 (str_progress_start, "already computed\n") )
        next
        }
      data_for_compute = input_data[, var_to_recomp, drop=F]
      n_vars = ncol (data_for_compute)
      }
    if (one_voi_name %in% var_of_interest_names)
      one_voi_values = data_for_compute[, one_voi_name]
    else
      one_voi_values = var_of_interest_values[, one_voi_name]
    #
    # Compute the mutual information by group of bin_size variables using miic
    #
    mat_mis[colnames(data_for_compute), one_voi_name] = 0
    start_idx <- 1
    while (start_idx < n_vars)
      {
      end_idx <- min (start_idx + bin_size - 1, n_vars)
      # print(paste0 ("From ", start_idx, " to ", end_idx, " (n_vars=", n_vars, ")") )
      time_str = ""
      if (verbose >= 3)
        {
        curr_time = Sys.time()
        elapsed_time = as.numeric (curr_time - time_start, units="secs")
        curr_progress = ( (one_voi_idx-1) + (start_idx - 1) / n_vars) / n_all_vois
        if (curr_progress > 0)
          {
          remain_time = (elapsed_time / curr_progress) - elapsed_time
          if (remain_time >= 3600)
            {
            time_str = paste0 (remain_time %/% 3600, "h " )
            remain_time = remain_time - (remain_time %/% 3600) * 3600
            }
          if (remain_time >= 60)
            {
            time_str = paste0 (time_str, remain_time %/% 60, "m " )
            remain_time = remain_time - (remain_time %/% 60) * 60
            }
         time_str = paste0 (", ", time_str, round (remain_time), "s to go" )
         }
        }
      if (verbose >= 2)
        cat (paste0 (str_progress_start,
          format (round ( ((start_idx-1) / n_vars) * 100, 2), nsmall=2), " %",
          time_str, "                \r") )

      data_loop <- data_for_compute [, start_idx:end_idx, drop=FALSE]
      if ( ! is.data.frame(data_loop) )
        data_loop = as.data.frame (data_loop)
      if (one_voi_name %in% colnames (data_loop))
        {
        # print ("voi in data loop, before:")
        # pos_col = which (colnames (data_loop) == one_voi_name)
        # print (miic:::list_to_str (colnames (data_loop)[(pos_col-1):(pos_col+1)] ) )
        data_loop[ , one_voi_name] <- NULL
        mat_mis[one_voi_name, one_voi_name] = NA_real_
        # print ("passe voi in data loop, after:")
        # print (miic:::list_to_str (colnames (data_loop)[(pos_col-1):(pos_col+1)] ) )
        }

      if (!skip_cheks)
        {
        # Remove rows full of NAs and constant variables
        # (would generate warnings if sent to miic function)
        #
        count_vals = unlist (apply (data_loop, MARGIN=2, FUN=function (x) {
          length (unique (x[!is.na(x)] ) ) }) )
        data_loop = data_loop[, count_vals >= 2, drop=F]

        count_nas = apply (data_loop, MARGIN=1, FUN=function(x) { sum (is.na(x) ) } )
        data_loop = data_loop[ count_nas < ncol(data_loop), , drop=F]
        }

      # print (paste0 ("nrow: ", nrow (data_loop),
      #               ", ncol: ", ncol (data_loop) ) )
      #
      if ( (nrow (data_loop) > 0) && (ncol (data_loop) > 0) )
        {
        so <- data.frame ("var_names" = c (colnames(data_loop), "var_interest"),
                          "is_consequence" = c (rep(1, ncol(data_loop)), 0),
                          stringsAsFactors = FALSE)
        data_loop$var_interest <- one_voi_values
        miic_res <- miic (data_loop, state_order=so,
          orientation=F, latent="no", n_threads=n_threads, verbose=0)
        miic_res <- miic_res$summary
        rownames (miic_res) <- NULL
        rownames (miic_res)[miic_res$x != "var_interest"] <- (
          miic_res[miic_res$x != "var_interest", "x"] )
        rownames (miic_res)[miic_res$y != "var_interest"] <- (
          miic_res[miic_res$y != "var_interest", "y"] )
        mat_mis[rownames(miic_res), one_voi_name] = round (
          (miic_res$info_shifted / miic_res$n_xy_ai) / LN_2, 6)
        }
      start_idx <- start_idx + bin_size
      }
    if (verbose >= 1)
      cat (paste0 (str_progress_start, "100 %                           \n") )
    }

  if (verbose >= 1)
    cat (paste0 (length (all_voi_names),
                 " variables of interest evaluated.\n") )
  return (mat_mis)
  }

#-------------------------------------------------------------------------------
# grid_plot
#-------------------------------------------------------------------------------
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

#===============================================================================
# FUNCTIONS (exported)
#===============================================================================
# computeMutualInfo
#-------------------------------------------------------------------------------
#' Compute (conditional) mutual information
#' @description For discrete or categorical variables, the (conditional)
#' mutual information is computed using the empirical frequencies minus a
#' complexity cost (computed as BIC or with the Normalized Maximum Likelihood).
#' When continuous variables are present, each continuous variable is
#' discretized for each mutual information estimate so as to maximize the
#' mutual information minus the complexity cost (see Cabeli 2020).
#'
#' @details For a pair of continuous variables \eqn{X} and \eqn{Y}, the mutual
#' information \eqn{I(X;Y)} will be computed iteratively. In each iteration, the
#' algorithm optimizes the partitioning of \eqn{X} and then of \eqn{Y},
#' in order to maximize
#' \deqn{Ik(X_{d};Y_{d}) = I(X_{d};Y_{d}) - cplx(X_{d};Y_{d})}
#' where \eqn{cplx(X_{d}; Y_{d})} is the complexity cost of the corresponding
#' partitioning (see Cabeli 2020).
#' Upon convergence, the information terms \eqn{I(X_{d};Y_{d})}
#' and \eqn{Ik(X_{d};Y_{d})}, as well as the partitioning of \eqn{X_{d}}
#' and \eqn{Y_{d}} in terms of cutpoints, are returned.
#'
#' For conditional mutual information with a conditioning set \eqn{U}, the
#' computation is done based on
#' \deqn{
#'   Ik(X;Y|U) = 0.5*(Ik(X_{d};Y_{d},U_{d}) - Ik(X_{d};U_{d})
#'                  + Ik(Y_{d};X_{d},U_{d}) - Ik(Y_{d};U_{d})),
#' }
#' where each of the four summands is estimated separately.
#'
#' @references
#' \itemize{
#' \item Cabeli \emph{et al.}, PLoS Comput. Biol. 2020, \href{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007866}{Learning clinical networks from medical records based on information estimates in mixed-type data}
#' \item Affeldt \emph{et al.}, UAI 2015, \href{https://auai.org/uai2015/proceedings/papers/293.pdf}{Robust Reconstruction of Causal Graphical Models based on Conditional 2-point and 3-point Information}
#' }
#'
#' @param x [a vector]
#' The \eqn{X} vector that contains the observational data of the first variable.
#' @param y [a vector]
#' The \eqn{Y} vector that contains the observational data of the second variable.
#' @param df_conditioning [a data frame]
#' The data frame of the observations of the conditioning variables.
#' @param maxbins [an integer]
#' When the data contain continuous variables, the maximum number of bins
#' allowed during the discretization. A smaller number makes the computation
#' faster, a larger number allows finer discretization.
#' @param cplx [a string]
#' The complexity model:
#' \itemize{
#' \item["bic"] Bayesian Information Criterion
#' \item["nml"] Normalized Maximum Likelihood, more accurate complexity cost
#' compared to BIC, especially on small sample size.
#' }
#' @param n_eff [an integer]
#' The effective number of samples. When there is significant autocorrelation
#' between successive samples, you may want to specify an effective number of
#' samples that is lower than the total number of samples.
#' @param sample_weights [a vector of floats]
#' Individual weights for each sample, used for the same reason as the effective
#' number of samples but with individual weights.
#' @param is_continuous [a vector of booleans]
#' Specify if each variable is to be treated as continuous (TRUE) or discrete
#' (FALSE), must be of length `ncol(df_conditioning) + 2`, in the order
#' \eqn{X, Y, U1, U2, ...}. If not specified, factors and character vectors are
#' considered as discrete, and numerical vectors as continuous.
#' @param plot [a boolean]
#' Specify whether the resulting XY optimum discretization is to be plotted
#' (requires `ggplot2` and `gridExtra`).
#'
#' @return A list that contains :
#' \itemize{
#' \item cutpoints1: Only when \eqn{X} is continuous, a vector containing
#'   the cutpoints for the partitioning of \eqn{X}.
#' \item cutpoints2: Only when \eqn{Y} is continuous, a vector containing
#'   the cutpoints for the partitioning of \eqn{Y}.
#' \item n_iterations: Only when at least one of the input variables is
#'   continuous, the number of iterations it takes to reach the convergence of
#'   the estimated information.
#' \item iteration1, iteration2, ... Only when at least one of the input
#'   variables is continuous, the list of vectors of cutpoints of each
#'   iteration.
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
#-------------------------------------------------------------------------------
computeMutualInfo <- function(x, y,
                              df_conditioning = NULL,
                              maxbins = NULL,
                              cplx = c("nml", "bic"),
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

  input_data = data.frame(x, y)
  if (!is.null(df_conditioning)) {
    input_data <- data.frame(input_data, df_conditioning)
  }

  if (!is.null(sample_weights) && length(sample_weights) != nrow(input_data)) {
    stop(paste(
      "Differing number of rows between `sample_weights` and input data:",
      length(sample_weights),
      length(x)
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
    result$n_iterations <- niterations
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
    nameDist1 <- deparse(substitute(x))
    nameDist2 <- deparse(substitute(y))
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

#-------------------------------------------------------------------------------
# computeThreePointInfo
#-------------------------------------------------------------------------------
#' Compute (conditional) three-point information
#' @description Three point information is defined and computed as the
#' difference of mutual information and conditional mutual information, e.g.
#' \deqn{I(X;Y;Z|U) = I(X;Y|U) - Ik(X;Y|U,Z)}
#' For discrete or categorical variables, the three-point information is
#' computed with the empirical frequencies minus a complexity cost
#' (computed as BIC or with the Normalized Maximum Likelihood).
#'
#' @details For variables \eqn{X}, \eqn{Y}, \eqn{Z} and a set of conditioning
#' variables \eqn{U}, the conditional three point information is defined as
#' \deqn{Ik(X;Y;Z|U) = Ik(X;Y|U) - Ik(X;Y|U,Z)}
#' where \eqn{Ik} is the shifted or regularized conditional mutual information.
#' See \code{\link{computeMutualInfo}} for the definition of \eqn{Ik}.
#'
#' @references
#' \itemize{
#' \item Cabeli \emph{et al.}, PLoS Comput. Biol. 2020, \href{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007866}{Learning clinical networks from medical records based on information estimates in mixed-type data}
#' \item Affeldt \emph{et al.}, UAI 2015, \href{https://auai.org/uai2015/proceedings/papers/293.pdf}{Robust Reconstruction of Causal Graphical Models based on Conditional 2-point and 3-point Information}
#' }
#'
#' @param x [a vector]
#' The \eqn{X} vector that contains the observational data of the first variable.
#' @param y [a vector]
#' The \eqn{Y} vector that contains the observational data of the second variable.
#' @param z [a vector]
#' The \eqn{Z} vector that contains the observational data of the third variable.
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
#' \item["bic"] Bayesian Information Criterion
#' \item["nml"] Normalized Maximum Likelihood, more accurate complexity cost
#' compared to BIC, especially on small sample size.
#' }
#' @param n_eff [an integer]
#' The effective number of samples. When there is significant autocorrelation
#' between successive samples, you may want to specify an effective number of
#' samples that is lower than the total number of samples.
#' @param sample_weights [a vector of floats]
#' Individual weights for each sample, used for the same reason as the effective
#' number of samples but with individual weights.
#' @param is_continuous [a vector of booleans]
#' Specify if each variable is to be treated as continuous (TRUE) or discrete
#' (FALSE), must be of length `ncol(df_conditioning) + 3`, in the order
#' \eqn{X, Y, Z, U1, U2, ...}. If not specified, factors and character vectors
#' are considered as discrete, and numerical vectors as continuous.
#'
#' @return A list that contains :
#' \itemize{
#' \item i3: The estimation of (conditional) three-point information without the
#' complexity cost.
#' \item i3k: The estimation of (conditional) three-point information with the
#' complexity cost (\emph{i3k = i3 - cplx}).
#' \item i2: For reference, the estimation of (conditional) mutual information
#' \eqn{I(X;Y|U)} used in the estimation of \emph{i3}.
#' \item i2k: For reference, the estimation of regularized (conditional) mutual
#' information \eqn{Ik(X;Y|U)} used in the estimation of \emph{i3k}.
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
#' message("I(X;Y;Z) = ", res$i3)
#' message("Ik(X;Y;Z) = ", res$i3k)
#'
#' \donttest{
#' # Independence, conditional dependence : X -> Z <- Y
#' X <- runif(N)
#' Y <- runif(N)
#' Z <- X + Y + rnorm(N, sd = 0.1)
#' res <- computeThreePointInfo(X, Y, Z)
#' message("I(X;Y;Z) = ", res$i3)
#' message("Ik(X;Y;Z) = ", res$i3k)
#' }
#-------------------------------------------------------------------------------
computeThreePointInfo <- function(x, y, z,
                              df_conditioning = NULL,
                              maxbins = NULL,
                              cplx = c("nml", "bic"),
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

  input_data = data.frame(x, y, z)
  if (!is.null(df_conditioning)) {
    input_data <- data.frame(input_data, df_conditioning)
  }

  if (!is.null(sample_weights) && length(sample_weights) != nrow(input_data)) {
    stop(paste(
      "Differing number of rows between `sample_weights` and input data:",
      length(sample_weights),
      length(x)
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

