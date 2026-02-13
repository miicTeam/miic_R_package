#*******************************************************************************
# Filename   : selectFeatures.R                Creation date: 17 October 2024
#
# Description: Features selection based on Mutual Information (around vois)
#
# Author     : Franck SIMON
#*******************************************************************************

#===============================================================================
# INTERNAL FUNCTIONS FOR selectFeatures (around vois => sfa_xx)
#===============================================================================
# Check pre-computed MIs
# Params:
# - precomputed_mis: a matrix containing the MIs between the variables of
#   interest (as columns) and variables evaluated (as rows)
# Return:
# - checked precomputed_mis
#-------------------------------------------------------------------------------
sfa_check_precomputed_mis <- function (precomputed_mis)
  {
  if ( is.null (precomputed_mis) )
    return (NULL)

  if ( ( ! is.matrix(precomputed_mis) )
    && ( ! inherits(precomputed_mis, "Matrix") ) )
    miic_error ("parameters", "the precomputed MIs must be a matrix.")
  #
  # Empty matrix is not an issue, we start with a new one
  #
  if ( (nrow (precomputed_mis) == 0) || (ncol (precomputed_mis) == 0) )
    return (NULL)
  #
  # If not empty, then rownames and colnames are expected
  #
  if ( is.null (colnames (precomputed_mis) ) )
    miic_error ("parameters",
      "the precomputed MIs matrix must have column names")
  if ( is.null (rownames (precomputed_mis) ) )
    miic_error ("parameters",
      "the precomputed MIs matrix must have row names")
  #
  # NB: colum names can be other variables of interest
  # and row names can be on variables not present for this round
  # but present in a full version of input_data (e.g. we pre-computed
  # MIs for all the genes in input_data and we are now looking only at
  # transcription factors)
  # => No cross check with input_data or variables of interest
  #
  return (precomputed_mis)
  }

#-------------------------------------------------------------------------------
# sfa_check_vois
#-------------------------------------------------------------------------------
# Check variables of interest for selectionFeatures
# Params:
# - input_data: a data frame or a matrix (including sparse matrix)
#   Expected layout is samples as rows and variables as columns.
#   Column names must contain the names of the variables.
# - var_of_interest_names: names of variables of interest (vois)
#   for vois present in data
# - var_of_interest_values: a data frame with vois not in input_data,
#   column names are the names of the vois and it must have
#   the same number of row than input_data
# Return: a named list with 4 items
# - "var_of_interest_names": a vector with the names of vois in input_data
# - "var_of_interest_values": a data frame with the vois not in input_data
# - "extra_voi_names": a vector with the names of the vois not in input_data
# - "all_voi_names": a vector with the name of all vois (in + not in input_data)
#-------------------------------------------------------------------------------
sfa_check_vois <- function (
  input_data, var_of_interest_names, var_of_interest_values)
  {
  if ( is.null (var_of_interest_names) && is.null (var_of_interest_values) )
    miic_error ("parameters", "the name of the variable(s) of interest or a",
      " data frame with the variable(s) of interest values must be supplied.")

  if ( is.null (var_of_interest_names) )
    var_of_interest_names <- c()
  else
    {
    for (one_var_name in var_of_interest_names)
      if ( test_param_wrong_string (one_var_name, colnames(input_data)) )
        miic_error ("parameters", "Some of the variable of interest",
                    " names are incorrect or not in the input_data.")
    }

  extra_voi_names <- rep ("", 0)
  if ( ! is.null (var_of_interest_values) )
    {
    if ( ! is.data.frame (var_of_interest_values) )
      miic_error ("parameters",
        "the var_of_interest_values must be a data frame.")
    # Ensure we have a true data frame, i.e. not a tibble
    var_of_interest_values <- as.data.frame (var_of_interest_values)
    if (ncol (var_of_interest_values) <= 0)
      {
      if ( is.null (var_of_interest_names) )
        miic_error ("parameters",
          "the var_of_interest_values data frame has been supplied but is empty.")
      else
        miic_warning ("parameters",
          "the var_of_interest_values data frame has been supplied but is empty.")
      var_of_interest_values <- NULL
      }
    else if ( nrow (var_of_interest_values) != nrow (input_data) )
      miic_error ("parameters",
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
      poss_wrong_idx <- which ( extra_voi_names %in% colnames(input_data) )
      if (length (poss_wrong_idx) >= 1)
        {
        for (one_var_name in extra_voi_names[poss_wrong_idx])
          {
          one_var_voi_vals <- var_of_interest_values[,one_var_name]
          one_var_input_vals <- input_data[,one_var_name]
          #
          # If not the same values => error
          #
          if (any ( ( is.na (one_var_voi_vals) != is.na (one_var_input_vals) )
                  | (one_var_voi_vals[ !is.na(one_var_voi_vals) ] != one_var_input_vals[ !is.na(one_var_input_vals) ]) ) )
            miic_error ("parameters",
              "the variable ", one_var_name, " is present in input data",
              " and has been supplied in variables of interest values.")
          #
          # Identical values => just a warning, use input_data
          # and ignore variables of interest values
          #
          miic_warning ("parameters",
            "the variable ", one_var_name, " is present in input data",
            " and has been supplied in variables of interest values.")
          var_of_interest_names <- unique (c (var_of_interest_names,
                                              one_var_name) )
          var_of_interest_values[,one_var_name] <- NULL
          }
        extra_voi_names <- colnames (var_of_interest_values)
        if (ncol (var_of_interest_values) <= 0)
          var_of_interest_values <- NULL
        }
      }
    }
  all_voi_names <- unique ( c (var_of_interest_names, extra_voi_names) )
  return ( list ("var_of_interest_names"=var_of_interest_names,
                 "var_of_interest_values"=var_of_interest_values,
                 "extra_voi_names"=extra_voi_names,
                 "all_voi_names"=all_voi_names) )
  }

#-------------------------------------------------------------------------------
# sfa_get_tops
#-------------------------------------------------------------------------------
# Select n_features from a list of variables evaluated against user supplied
# variables of interest (vois). The number of features selected per voi is
# intended to be balanced (~ n_features / n_vois) but this objective can
# not be met if some vois have to few features to select.
# Params:
# - n_features: total number of features to select
# - list_sorted: named list of named lists. The main list contain one item per
#   user variable of interest, each sub lists contains MIs of variables
#   evaluated, sorted in descending order, names of the sub lists are the
#   variables names.
# Return:
# - a vector with the features selected
#-------------------------------------------------------------------------------
sfa_get_tops <- function (n_features, list_sorted, verbose=3)
  {
  if (n_features <= 0)
    return ( c() )
  #
  # Select top features until we get enough or a bit too much
  #
  list_tops <- c()
  list_tops_prec <- c()
  n_vois <- length (list_sorted)
  n_tops_to_sel <- n_features %/% n_vois
  while (T)
    {
    list_tops <- unique (unlist (lapply (list_sorted, FUN=function(x) {
      if (length(x) > 0)
        names(x)[1:min(length(x), n_tops_to_sel)]
      } ) ) )
    n_missing <- n_features - length (list_tops)
    if (n_missing <= 0)
      break
    list_tops_prec <- list_tops
    n_tops_to_sel <- n_tops_to_sel + max (1, n_missing %/% n_vois)
    }
  #
  # Prepare the last positions used per voi for future check and warning
  #
  last_pos_per_voi <- rep (n_tops_to_sel, n_vois)
  #
  # If a bit too much, select in the last features added those with highest MI
  #
  if (n_missing < 0)
    {
    list_tops_to_test <- unlist (lapply (list_sorted,
      FUN=function(x) {
        if (length(x) > 0)
          names(x)[1:min(length(x), n_tops_to_sel)]
        } ) )
    list_vals_to_test <- unlist (lapply (list_sorted,
      FUN=function(x) {
        if (length(x) > 0)
          x[1:min(length(x), n_tops_to_sel)]
        } ) )
    list_tops_to_test <- list_tops_to_test[order (list_vals_to_test,
                                                  decreasing=T)]
    list_tops_to_test_uniq <- unique (list_tops_to_test)
    list_tops_to_test_uniq <- list_tops_to_test_uniq[
      ! (list_tops_to_test_uniq %in% list_tops_prec) ]

    n_missing <- n_features - length (list_tops_prec)
    tops_to_add <- list_tops_to_test_uniq[1:n_missing]
    list_tops <- c (list_tops_prec, tops_to_add)
    #
    # Decrease the last position for the variable not added
    #
    col_added <- unlist (lapply (list_sorted,
      FUN=function(x) {
        ifelse (length(x) < n_tops_to_sel, F,
                names(x)[n_tops_to_sel] %in% tops_to_add)
        } ) )
    last_pos_per_voi[!col_added] <- last_pos_per_voi[!col_added] - 1
    #
    # Warning if, for the last feature selected, a variable of another vois
    # has the same MI but was not selected
    #
    # list_tops_to_test was order by MI desc, order values in the same way
    list_vals_to_test <- sort (list_vals_to_test, decreasing=T)
    # remove variables/values selected before
    list_vals_to_test <- list_vals_to_test[
      ! (list_tops_to_test %in% list_tops_prec) ]
    list_tops_to_test <- list_tops_to_test[
      ! (list_tops_to_test %in% list_tops_prec) ]
    # min MI of variables selected
    min_mi <- min (list_vals_to_test[list_tops_to_test %in% tops_to_add])
    cnt_added <- length (which (
      list_vals_to_test[list_tops_to_test %in% tops_to_add] == min_mi) )
    # remove variables/values selected jsut above
    list_vals_to_test <- list_vals_to_test[! (list_tops_to_test %in% tops_to_add) ]
    list_tops_to_test <- list_tops_to_test[! (list_tops_to_test %in% tops_to_add) ]
    cnt_not_added <- length (which (list_vals_to_test == min_mi) )
    if (cnt_not_added > 0)
      miic_warning ("Features selection",
        "The last feature selected in the top ", n_tops_to_sel, "\n",
        "for the different variables of interest has a MI of ",
        round (min_mi,4), ".\n",
        cnt_added + cnt_not_added, " variables have the same MI",
        " for different variables of interest\n",
        "but to respect the number of features requested,\n",
        cnt_added, " have been selected while ",
        cnt_not_added, " have been discarded.\n",
        "You can check the returned MIs to include",
        " the discarded equivalent features.")
    }
  #
  # Issue warnings if, after the last position, we have other variables
  # with same MI
  #
  for (i in 1:n_vois)
    {
    vals_sorted_voi <- list_sorted[[i]]
    last_pos_voi <- last_pos_per_voi[[i]]
    if ( last_pos_voi >= length(vals_sorted_voi) )
      next
    last_val <- vals_sorted_voi[[last_pos_voi]]
    have_same_val_after <- (vals_sorted_voi [(last_pos_voi+1):length(vals_sorted_voi)] == last_val)
    if (sum (have_same_val_after) <= 0)
      next
    have_same_val_before <- (vals_sorted_voi [1:last_pos_voi] == last_val)
    miic_warning ("Features selection",
      paste0 ("The last feature selected for variable of interest ",
      names(list_sorted)[[i]], " has a MI of ", round(last_val,4), ".\n",
      sum (have_same_val_before) + sum(have_same_val_after),
      " variables have the same MI",
      " but to respect the number of features requested,\n",
      sum (have_same_val_before), " have been selected while ",
      sum (have_same_val_after), " have been discarded.\n",
      "You can check the returned MIs to include",
      " the discarded equivalent features.") )
    }
  if (verbose >= 1)
    miic_msg (length(list_tops), " features selected.")
  return (list_tops)
  }

#-------------------------------------------------------------------------------
# sfa_plot
#-------------------------------------------------------------------------------
# Barplot of the top features for one variable of interest (voi)
# Params:
# - mis: a matrix containing the MIs between the variables of
#   interest (as columns) and variables evaluated (as rows)
# - var_of_interest_name: a string, the voi to plot
# - n_plots: an integer, default 25, the number of features to plot
# - x_lab: a string, the x label, default "Top features for voi name"
# - y_lab: a string, the y label, default defined from unit and corrected
# - unit: a string, default "log_conf", possible values: "log_conf", "bits"
# - corrected: a boolean, default T, indicate if MIs include a correction
# - values: a boolean, default T, if T, display the MI of each feature plotted
# - annotate: a boolean, default T, if T, add an annotation about the number
#   of features
# - font_size: a boolean, default 11, the font size
# - box_fill: a string, default "#1F78B4", the box filling color
# - box_text: a string, default "white", the box text color
# Return: a ggplot2 barplot
#-------------------------------------------------------------------------------
sfa_plot <- function (mis, var_of_interest_name, n_plots=25,
  x_lab=NULL, y_lab=NULL, unit="log_conf", corrected=T, values=T, annotate=T,
  font_size=11, box_fill="#1F78B4", box_text="white")
  {
  # Check the parameters that can have been supplied by the user
  # (passed by the ... extra params of selectFeatures function)
  #
  n_plots <- check_param_int (
    n_plots, "number of features to plot", default=25, min=1)
  if ( is.null (x_lab) )
    x_lab <- paste0 ("Top features for ", var_of_interest_name)
  else
    x_lab <- as.character (x_lab)
  if ( is.null (y_lab) )
    {
    if (unit == "bits")
      y_lab <- ifelse (corrected, "MI' (bits)", "MI (bits)")
    else
      y_lab <- ifelse (corrected, "MI' (corrected log confidence)",
                                  "MI (uncorrected log confidence)")
    }
  else
    y_lab <- as.character (y_lab)
  values <- check_param_logical (values, "plotting of values", default=T)
  annotate <- check_param_logical (
    annotate, "plotting of annotation", default=T)
  font_size <- check_param_int (font_size, "font size", default=11, min=1)

  # TODO ? add a function for color checking
  if ( is.null (box_fill) )
    box_fill <- "#1F78B4"
  if ( is.null (box_text) )
    box_text <- "white"
  #
  # Extract and order desc the MIs for the var of interest to plot
  #
  v_mis <- mis[, var_of_interest_name, drop=T]
  v_mis <- v_mis[ (!is.na (v_mis)) & (v_mis > 0) ]
  v_mis <- sort (v_mis, decreasing=T)
  df <- data.frame ("Features"=names(v_mis), "MI"=v_mis)
  if (nrow (df) > n_plots)
    df <- df[1:n_plots, , drop=F]
  mi_min <- ifelse (nrow (df) > 0, min(df[,"MI"]), 0)
  mi_max <- ifelse (nrow (df) > 0, max(df[,"MI"]), 0)

  if (mi_max == 0)
    p <- ggplot2::ggplot () +
        ggplot2::theme_classic() +
        ggplot2::theme ( text=ggplot2::element_text(size=font_size),
                         axis.text.x=ggplot2::element_blank(),
                         axis.ticks.x=ggplot2::element_blank(),
                         axis.text.y=ggplot2::element_blank(),
                         axis.ticks.y=ggplot2::element_blank() ) +
        ggplot2::xlab (x_lab) +
        ggplot2::ylab (y_lab) +
        ggplot2::xlim (0, 1) +
        ggplot2::ylim (0, 1)
  else
    # with is a bad fix to avoid warnings from ggplot2 or notes from CRAN checks
    p <- with (df, ggplot2::ggplot (,
          mapping=ggplot2::aes (x=Features, y=MI) ) +
        ggplot2::geom_bar (stat="identity", fill=box_fill) +
        ggplot2::scale_x_discrete (limits=df$Features ) +
        ggplot2::theme_classic() +
        ggplot2::theme ( text=ggplot2::element_text(size=font_size),
          axis.text.x=ggplot2::element_text (angle=30, vjust=1, hjust=1) ) +
        ggplot2::xlab (x_lab) +
        ggplot2::ylab (y_lab) )

  if (values)
    {
    if (mi_min < 0.1)
      values_rounded <- round (df$MI, 3)
    else if (mi_min < 1)
      values_rounded <- round (df$MI, 2)
    else if (mi_min < 10)
      values_rounded <- round (df$MI, 1)
    else
      values_rounded <- round (df$MI, 0)
    p <- p + ggplot2::geom_text (ggplot2::aes (label=values_rounded),
      color=box_text, vjust=1.6, size=font_size*0.8/ggplot2::.pt )
    }

  if (annotate)
    {
    if (length(v_mis) > n_plots)
      label_text <- paste0 ("Top ", n_plots, " shown of ",
                            length(v_mis), " features with MI > 0")
    else
      label_text <- paste0 (length(v_mis), " features with MI > 0")
    if (mi_max == 0)
      annotation <- data.frame (x=1, y=1, label=label_text)
    else
      annotation <- data.frame (x=nrow(df), y=mi_max, label=label_text)

    # with is a bad fix to avoid warnings from ggplot2 or notes from CRAN checks
    p <- with (annotation,
      p + ggplot2::geom_text (annotation,
        mapping=ggplot2::aes(x=x, y=y, label=label),
        size=font_size*0.8/ggplot2::.pt, hjust=1, vjust=1) )
    }
  return (p)
  }

#===============================================================================
# FUNCTIONS (exported)
#===============================================================================
# selectFeatures
#-------------------------------------------------------------------------------
#' Features selection based on Mutual Information (MI)
#'
#' @description Select features by retrieving the set of variables having
#' the highest MIs with the variable(s) of interest.
#'
#' When multiple variables of interest are supplied, the function tries
#' to select an equal number of features per variable of interest.
#' Variables not sharing information at all (MI = 0) can not be selected
#' as features, so the number of features returned globally or for a variable
#' of interest can be lower than requested. This can also lead to
#' an unbalanced number of features per variable of interest.
#'
#' @param input_data [a data frame or a matrix, required]
#'
#' Expected layout is samples as rows and variables as columns.
#' Column names must contain the names of the variables.
#'
#' @param n_features [an integer, required]
#'
#' The total number of features to select, must be in the range
#' [0, number of variables in \emph{input_data}].
#' If several variables of interest are supplied,
#' the function tries to select an equal number of features per variable
#' and \emph{n_features} (if > 0) is expected to be greater or equal
#' than the number of variables of interest
#'
#' When \emph{n_features} is 0, no feature selection is performed
#' but the MIs between the variables in \emph{input_data} and the variable(s)
#' of interest are computed and returned.
#' This can be used to compute once for all the MIs and perform several
#' features selection at a later stage, supplying the MIs matrix in the
#' \emph{precomputed_mis} parameter.
#'
#' @param var_of_interest_names [a string or vector of strings, optional,
#' NULL by default]
#'
#' For the variable(s) of interest that are part of the \emph{input_data},
#' you should supply their names here.
#'
#' @param var_of_interest_values [a data frame, optional, NULL by default]
#'
#' For the variables of interest that are not in \emph{input_data},
#' a data frame can be supplied. The column names are the names
#' of the variables of interest and rows are the samples
#' ordered in the same way as \emph{input_data}.
#' Typically, such variables are metadata associated to samples but not
#' stored in \emph{input_data}, e.g. a "Treatment" vs "Control" variable
#' in an experiment and a count matrix with the expression of genes
#' in \emph{input_data}.
#'
#' @param unit [a string, optional, "log_conf" by default]
#'
#' Indicates the "unit" of MIs returned and used for plotting if requested.
#' Possible values are "log_conf" or "bits".
#'
#' @param corrected [a boolean, optional, TRUE by default]
#'
#' When set to TRUE, the mutual information values are corrected by subtracting
#' a complexity term (computed with the Normalized Maximum Likelihood).
#' For datasets having very few samples, the complexity term can have
#' a disproportionate impact. Setting \emph{corrected} to FALSE switches
#' to the use of non corrected mutual information.
#'
#' @param precomputed_mis [a matrix, optional, NULL by default]
#'
#' if MIs have been previously computed between some variables
#' in the \emph{input_data} and variable(s) of interest,
#' supplying these precomputed MIs speeds up the process as the existing MIs
#' (values present and different from NA) are not recomputed.
#' This matrix must have variables names from the \emph{input_data}
#' as row names and variables of interest names as column names
#' (the layout is the same as the \emph{mis} matrix returned).
#' To be valid, the pre-computed MIs must have been computed using
#' the same \emph{unit} and \emph{corrected} parameters.
#'
#' @param skip_cheks [a boolean, optional, FALSE by default]
#'
#' Before computing MIs between the variable(s) of interest and the features,
#' \emph{input_data} is checked to filter out constant features and rows full
#' of NAs. When the \emph{input_data} does not need such filtering,
#' these checks can be skipped to speed up the process.
#'
#' @param n_threads [a positive integer, optional, 1 by default]
#'
#' When set greater than 1, \emph{n_threads} parallel threads are used for
#' computation. Make sure your compiler is compatible with openmp
#' if you wish to use multithreading.
#'
#' @param verbose [an integer, optional, 3 by default]
#'
#' Level of verbosity: 0=no display, 1=summary, 2=progress per variable of
#' interest, 3=same as 2 with display of estimated time remaining.
#'
#' @param plot [a boolean, optional, FALSE by default]
#'
#' If set to TRUE, a plot with the top features for each variable of
#' interest is generated (requires `ggplot2`).
#'
#' @param ...
#'
#' If plotting is requested, extra parameters can be used to customize the plot
#' rendering:
#'
#' \itemize{
#' \item \emph{n_plots}: an integer, optional, 25 by default, the number of
#'   features to plot.
#' \item \emph{x_lab}: a string, optional, "Top features for
#'   variable_of_interest_name" by default, the X axis label.
#' \item \emph{y_lab}: a string, optional, determined from the \emph{unit}
#'   and \emph{corrected} parameters by default, the Y axis label.
#' \item \emph{values}: a boolean, optional, TRUE by default,
#'   displays the MI of each feature when activated.
#' \item \emph{annotate}: a boolean, optional, TRUE by default,
#'   displays the shown/total number of features with MI > 0 when activated.
#' \item \emph{font_size}: an integer, optional, 11 by default, the font size.
#' \item \emph{box_fill}: a string, optional, "#1F78B4" by default
#'   (darker blue of the Paired brewer palette), the bar plot filling color.
#' \item \emph{box_text}: a string, optional, "white" by default,
#'   the text color of MI values in the bar plot.
#' }
#'
#' @return A named list with three items:
#'
#' \itemize{
#' \item \emph{features}: a vector with the features selected.
#' \item \emph{mis}: a matrix with the MIs between the variables in
#'   \emph{input_data} (as rows) and the variable(s) of interest (as columns).
#'   Row and column names are sorted alphabetically.\cr
#'   The MIs can be expressed, depending on the \emph{unit} parameter, as
#'   log confidence or in bits and, depending on the \emph{corrected} parameter,
#'   include a correction or not.\cr
#'   If a precomputed MIs matrix was supplied, new values computed
#'   are added to the existing matrix.
#' \item \emph{plots}: when the \emph{plot} parameter is turned to TRUE,
#'   a list of plots with the top features for each variable of interest.
#'   The list is empty if the \emph{plot} parameter is set to FALSE.
#' }
#'
#' @examples
#' library(miic)
#'
#' # Get the 10 top features related to "tp53" (gene mutation)
#' # and "TP53" (gene expression) in the dataset
#' ret <- selectFeatures (cosmicCancer,
#'                        var_of_interest_names=c("TP53", "tp53"),
#'                        n_features=10)
#' message ("Tops 10 features: ", paste (ret$features, collapse=", ") )
#'
#' # Features selection with reuse of above computed MIs
#' ret <- selectFeatures (cosmicCancer,
#'                        var_of_interest_names=c("TP53", "tp53"),
#'                        n_features=20,
#'                        precomputed_mis=ret$mis)
#' message ("Tops 20 features: ", paste (ret$features, collapse=", ") )
#'
#' # Features selection using an external metadata "Ploidy"
#' # (simulation by extracting "Ploidy" out of the dataset)
#' df_external_meta <- data.frame ("Ploidy"=cosmicCancer$Ploidy)
#' df_data <- cosmicCancer[ , ! (colnames(cosmicCancer) == "Ploidy") ]
#' ret <- selectFeatures (df_data,
#'                        var_of_interest_values=df_external_meta,
#'                        n_features=10)
#' message ("Tops 10 features: ", paste (ret$features, collapse=", ") )
#'
#' # Same features selection with plotting of the result
#' ret <- selectFeatures (df_data,
#'                        var_of_interest_values=df_external_meta,
#'                        n_features=10,
#'                        plot=TRUE)
#' print (ret$plots[[1]])
#'
#' # Features selection using variables of interest coming both from
#' # the dataset ("tp53", "TP53") and an external metadata ("Ploidy")
#' # (simulation by extracting "Ploidy" out of the dataset)
#' df_external_meta <- data.frame ("Ploidy"=cosmicCancer$Ploidy)
#' df_data <- cosmicCancer[ , ! (colnames(cosmicCancer) == "Ploidy") ]
#' ret <- selectFeatures (df_data,
#'                        var_of_interest_names=c("TP53", "tp53"),
#'                        var_of_interest_values=df_external_meta,
#'                        n_features=20)
#' message ("Tops 20 features: ", paste (ret$features, collapse=", ") )
#'
#' # Same features selection with plotting using a customized rendering
#' if ( require(gridExtra) ) {
#'   ret <- selectFeatures (df_data,
#'                          var_of_interest_names=c("TP53", "tp53"),
#'                          var_of_interest_values=df_external_meta,
#'                          n_features=20,
#'                          plot=TRUE,
#'                          n_plots=10,
#'                          box_fill="lightgreen",
#'                          box_text="black")
#'   do.call ("grid.arrange", c(ret$plots, ncol=2,
#'     top="Top features on cosmicCancer dataset" ) )
#' }
#'
#' @export
#-------------------------------------------------------------------------------
selectFeatures <- function (input_data, n_features, var_of_interest_names=NULL,
  var_of_interest_values=NULL, unit="log_conf", corrected=T,
  precomputed_mis=NULL, skip_cheks=F, n_threads=1, verbose=3, plot=F, ...)
  {
  # Check parameters
  #
  input_data <- sf_check_input_data (input_data)
  n_features <- check_param_int (
    n_features, "number of features", default=0, min=0)
  vois <- sfa_check_vois (input_data,
    var_of_interest_names, var_of_interest_values)
  unit <- check_param_string (unit, "unit", c("log_conf", "bits") )
  corrected <- check_param_logical (corrected, "corrected", T)
  precompured_mis <- sfa_check_precomputed_mis (precomputed_mis)
  skip_cheks <- check_param_logical (skip_cheks, "skip checks", F)
  n_threads <- check_param_int (n_threads, "number of threads", 1, min=1)
  verbose <- check_param_int (verbose, "verbose", 3, min=0, max=3)
  plot <- check_param_logical (plot, "plot", F)
  #
  # Cross checks
  #
  if ( n_features > ncol(input_data) )
    {
    n_features <- ncol(input_data)
    miic_warning ("parameter",  "the number of features can not be greater",
      " than the number of variables in input data.",
      " It has been reduced to ", n_features, ".")
    }
  if ( (n_features > 0) && ( n_features < length (vois$all_voi_names) ) )
    {
    n_features <- length (vois$all_voi_names)
    miic_warning ("parameter",  "the number of features, if not 0,",
      " must be >= number of variables of interest.",
      " It has been increased to ", n_features, ".")
    }
  #
  # Compute MIs
  #
  mat_mis <- compute_mi_batch (input_data=input_data,
    var_of_interest_names=vois$var_of_interest_names,
    var_of_interest_values=vois$var_of_interest_values,
    unit=unit, corrected=corrected, precomputed_mis=precomputed_mis,
    skip_cheks=skip_cheks, n_threads=n_threads, verbose=verbose)
  #
  # The mat_mis can contain more rows than the features to select
  # e.g. we pre-computed the MI for all genes and now we want select only
  # the TFs. In this case, in input_data, the variables are only the TFs
  # while mat_mis would contain all genes. Same for the columns as
  # we can have precomputed more variables of interest than the ones we use now
  #
  mat_mis_filt <- mat_mis[rownames(mat_mis) %in% colnames(input_data),
                          colnames(mat_mis) %in% vois$all_voi_names,
                          drop=F]
  #
  # Plot if requested
  #
  list_plots <- list ()
  # one_col <- colnames (mat_mis_filt)[[1]]
  if (plot)
    {
    if ( base::requireNamespace("ggplot2", quietly=TRUE) )
      for ( one_col in colnames (mat_mis_filt) )
        list_plots[[one_col]] <- sfa_plot (
          mis=mat_mis_filt, var_of_interest_name=one_col,
          unit=unit, corrected=corrected, ...)
    else
      miic_warning ("Features selection", "Plotting requires ggplot2.")
    }
  #
  # If no feature selection, can end here
  #
  if (n_features <= 0)
    return (list ("features"=c(), "mis"=mat_mis, "plots"=list_plots) )
  #
  # For each voi requested, order mis decreasing and filter on MI' > 0
  #
  list_mis_sorted <- list()
  uniq_poss_feats <- c()
  for (one_voi in vois$all_voi_names)
    {
    list_mis_sup_0 <- mat_mis_filt[, one_voi] [mat_mis_filt[, one_voi] > 0]
    list_mis_sorted[[one_voi]] <- sort (list_mis_sup_0, decreasing=T)
    uniq_poss_feats <- unique (c (uniq_poss_feats,
                                  names(list_mis_sup_0) ) )
    }
  #
  # Check if too much features to select
  #
  if (n_features >= length(uniq_poss_feats) )
    {
    if (n_features > length(uniq_poss_feats) )
      miic_warning ("features selection",
        "too few variables share information with the variable(s)",
        " of interest to select the top ", n_features,
        ", returning all the ", length(uniq_poss_feats),
        " variables sharing information.")
    return (list ("features"=uniq_poss_feats, "mis"=mat_mis, "plots"=list_plots) )
    }
  #
  # Select top features until we get enough or a bit too much
  #
  list_tops <- sfa_get_tops (
    n_features=n_features, list_sorted=list_mis_sorted, verbose=verbose)

  return (list ("features"=list_tops, "mis"=mat_mis, "plots"=list_plots) )
  }
