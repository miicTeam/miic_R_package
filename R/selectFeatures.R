#*******************************************************************************
# Filename   : selectFeatures.R                 Creation date: 17 October 2024
#
# Description: Features selection based on Mutual Information
#
# Author     : Franck SIMON
#
# TODO ? warning discrete number of levels (as miic) ?")
#*******************************************************************************

################################################################################
# FUNCTIONS (internal)
################################################################################
# COMMON FUNCTIONS TO selectFeatures AND selectFeaturesPath (sf_xx)
#===============================================================================
# sf_check_input_data
#-------------------------------------------------------------------------------
# Check input data
# Params:
# - input_data: a data frame or a matrix (including sparse matrix)
#   Expected layout is samples as rows and variables as columns.
#   Column names must contain the names of the variables.
#-------------------------------------------------------------------------------
sf_check_input_data <- function (input_data)
  {
  if (  ( ! is.data.frame (input_data) )
     && ( ! is.matrix(input_data) )
     && ( ! inherits(input_data, "Matrix") ) )
    miic_error ("parameters",
      "the input data must be a data frame or a matrix.")
  if ( (ncol (input_data) <= 0) || (nrow (input_data) <= 0) )
    miic_error ("parameters", "the input data is empty.")
  if ( is.data.frame (input_data) )
    # Ensure we have a true data frame, e.g. not a tibble
    # (but let matrices unchanged to avoid warnings on large memory allocation)
    # TODO evaluate run time impact on very large data frames
    input_data <- as.data.frame (input_data)
  if ( is.null (colnames (input_data) ) )
    miic_error ("parameters", "the input data must have column names.")
  return (input_data)
  }

#-------------------------------------------------------------------------------
# sf_check_precomputed_mis
#-------------------------------------------------------------------------------
# Check pre-computed MIs
# Params:
# - precomputed_mis: a matrix containing the MIs between the variables of
#   interest (as columns) and variables evaluated (as rows)
# Return:
# - checked precomputed_mis
#-------------------------------------------------------------------------------
sf_check_precomputed_mis <- function (precomputed_mis)
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

#===============================================================================
# FUNCTIONS FOR selectFeatures only (sfo_xx)
#===============================================================================
# sfo_check_vois
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
sfo_check_vois <- function (
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
# sfo_get_tops
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
sfo_get_tops <- function (n_features, list_sorted, verbose=3)
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
    min_mi = min (list_vals_to_test[list_tops_to_test %in% tops_to_add])
    cnt_added = length (which (
      list_vals_to_test[list_tops_to_test %in% tops_to_add] == min_mi) )
    # remove variables/values selected jsut above
    list_vals_to_test <- list_vals_to_test[! (list_tops_to_test %in% tops_to_add) ]
    list_tops_to_test <- list_tops_to_test[! (list_tops_to_test %in% tops_to_add) ]
    cnt_not_added = length (which (list_vals_to_test == min_mi) )
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
# sfo_plot
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
sfo_plot <- function (mis, var_of_interest_name, n_plots=25,
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
# FUNCTIONS FOR selectFeaturesPath (sfp_xx)
#===============================================================================
# sfp_check_vois
#-------------------------------------------------------------------------------
# Check variables of interest for selectFeaturesPath
# Params:
# - input_data: a data frame or a matrix (including sparse matrix)
#   Expected layout is samples as rows and variables as columns.
#   Column names must contain the names of the variables.
# - var_of_interest_names_side1: names of variables of interest (vois)
#   for vois present in data
# - var_of_interest_values_side1: a data frame with vois not in input_data,
#   column names are the names of the vois and it must have
#   the same number of row than input_data
# - var_of_interest_names_side2: same as var_of_interest_names_side2
#   for the other side of the path
# - var_of_interest_values_side2: same as var_of_interest_values_side1
#   for the other side of the path
# Return: list of 3 items, "side1", "side2" and "all".
# "side1" and "side2" are nested lists with:
# - "var_of_interest_names": a vector with the names of vois in input_data
# - "var_of_interest_values": a data frame with the vois not in input_data
# - "extra_voi_names": a vector with the names of the vois not in input_data
# - "all_voi_names": a vector with the name of all vois (in + not in input_data)
# "all" is a nested list with:
# - "all_voi_names: a vector with all the voi names (in + not in input_data)
#   from all sides
# - "var_of_interest_names": a vector with the names of vois in input_data
#   from all sides
# - "extra_voi_names": a vector with the names of the vois not in input_data
#   from all sides
#-------------------------------------------------------------------------------
sfp_check_vois <- function (input_data,
    var_of_interest_names_side1, var_of_interest_values_side1,
    var_of_interest_names_side2, var_of_interest_values_side2)
  {
  # Check VOIs of each side
  #
  vois <- list ()
  for (i in 1:2)
    {
    if (i == 1)
      {
      var_of_interest_names <- var_of_interest_names_side1
      var_of_interest_values <- var_of_interest_values_side1
      }
    else
      {
      var_of_interest_names <- var_of_interest_names_side2
      var_of_interest_values <- var_of_interest_values_side2
      }
    #
    # Same kind of tests as compute_mi_batch
    #
    if ( is.null (var_of_interest_names) && is.null (var_of_interest_values) )
      miic_error ("parameters", "the name of the variable(s) of",
        " interest or a data frame with the variable(s) of interest values",
        " must be supplied for side ", i, ".")

    if (is.null (var_of_interest_names) )
      var_of_interest_names <- c()
    else
      {
      for (one_var_name in var_of_interest_names)
        if ( test_param_wrong_string (one_var_name, colnames(input_data) ) )
          miic_error ("parameters",  "Some of the variable of interest",
            " names for side ", i, " are incorrect or not in the input_data.")
      }

    extra_voi_names <- rep ("", 0)
    if ( ! is.null (var_of_interest_values) )
      {
      if ( ! is.data.frame (var_of_interest_values) )
        miic_error ("parameters",
          "the var_of_interest_values for side ", i, " must be a data frame.")
      # Ensure we have a true data frame, i.e. not a tibble
      var_of_interest_values <- as.data.frame (var_of_interest_values)
      if (ncol (var_of_interest_values) <= 0)
        {
        if ( is.null (var_of_interest_names) )
          miic_error ("parameters",
            "the var_of_interest_values data frame for side ", i,
            " has been supplied but is empty.")
        else
          miic_warning ("parameters",
            "the var_of_interest_values data frame for side ", i,
            " has been supplied but is empty.")
        var_of_interest_values <- NULL
        }
      else if ( nrow (var_of_interest_values) != nrow (input_data) )
        miic_error ("parameters",
          "the variable of interest values for side ", i,
          " does not match the number of samples.")
      else
        {
        # Data frame OK, checks variables names not in data
        #
        extra_voi_names <- colnames (var_of_interest_values)
        #
        # Error or warning if voi requested as external in
        # var_of_interest_values are present in the input data
        #
        poss_wrong_idx <- which ( extra_voi_names %in% colnames(input_data) )
        if (length (poss_wrong_idx) >= 1)
          {
          for (one_var_name in extra_voi_names[poss_wrong_idx])
            {
            one_var_voi_vals <- var_of_interest_values[,one_var_name]
            one_var_input_vals <- input_data[,one_var_name]

            if (any ( ( is.na (one_var_voi_vals) != is.na (one_var_input_vals) )
                    | (one_var_voi_vals[ !is.na(one_var_voi_vals) ] != one_var_input_vals[ !is.na(one_var_input_vals) ]) ) )
              miic_error ("parameters",
                "the variable ", one_var_name, " is present both in input data",
                " and in variables of interest values of side ", i, ".")
            #
            # Supplied in both in input data and in variables of interest values
            # and with identical values => just a warning, use input_data
            # and ignore variables of interest values
            #
            miic_warning ("parameters",
              "the variable ", one_var_name, " is present both in input data",
              " and in variables of interest values of side ", i, ".")
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

    vois[[paste0 ("side", i)]] = list ("var_of_interest_names" = var_of_interest_names,
      "var_of_interest_values" = var_of_interest_values,
      "extra_voi_names" = extra_voi_names,
      "all_voi_names" = c(var_of_interest_names, extra_voi_names) )
    }
  #
  # Test on each side done, now test one side against the other
  #
  all_voi_names <- c(vois[[1]]$var_of_interest_names, vois[[1]]$extra_voi_names,
                     vois[[2]]$var_of_interest_names, vois[[2]]$extra_voi_names)
  are_duplicated <- duplicated (all_voi_names)
  if (any (are_duplicated))
    miic_error ("parameters",
      "Some variable(s) have been supplied in both sides: ",
      list_to_str (all_voi_names[are_duplicated], n_max=10), ".")
  #
  # Create variables that store the vois of both sides
  #
  var_of_interest_names_all_sides <- unique (c(vois[[1]]$var_of_interest_names,
                                               vois[[2]]$var_of_interest_names))
  if (  is.null (vois[[1]]$var_of_interest_values)
     && is.null (vois[[2]]$var_of_interest_values) )
    var_of_interest_values_all_sides <- NULL
  else if (  is.null (vois[[1]]$var_of_interest_values) )
    var_of_interest_values_all_sides <- vois[[2]]$var_of_interest_values
  else if (  is.null (vois[[2]]$var_of_interest_values) )
    var_of_interest_values_all_sides <- vois[[1]]$var_of_interest_values
  else
    var_of_interest_values_all_sides <- cbind (vois[[1]]$var_of_interest_values,
                                               vois[[2]]$var_of_interest_values)
  extra_voi_names_all_sides <- unique (c (vois[[1]]$extra_voi_names,
                                          vois[[2]]$extra_voi_names) )

  vois[["all"]] = list (
    "var_of_interest_names" = var_of_interest_names_all_sides,
    "var_of_interest_values" = var_of_interest_values_all_sides,
    "extra_voi_names" = extra_voi_names_all_sides,
    "all_voi_names" = all_voi_names)
  return (vois)
  }

#-------------------------------------------------------------------------------
# sfp_prepare_couples
#-------------------------------------------------------------------------------
# Prepare a data frame with the couples between the variables of interest (vois)
# from each side of the path. For each couple, the function will pick the MI
# from the mis matrix if present, and if not (case between 2 vois not in
# input_data), call computeMutualInfo
# Params:
# - vois: the list returned by sfp_check_vois
# - mis: a matrix containing the MIs between the variables of
#   interest (as columns) and variables evaluated (as rows)
# - unit: a string, possible values "log_conf", "bits"
# - corrected: a boolean, indicated if the MI includes a correction
# Return:
# - a data frame with all the couples as rows. Information available are
#   "x" (first side of the couple), "y" (first side) and "mi"
#-------------------------------------------------------------------------------
sfp_prepare_couples <- function (vois, mis, unit, corrected)
  {
  LN_2 <- log (2)
  couples <- expand.grid (vois$side1$all_voi_names, vois$side2$all_voi_names,
                          stringsAsFactors=F)
  colnames (couples) <- c ("x", "y")
  couples$mi <- unlist (apply (couples, MARGIN=1, function (x) {
    if ( ! (x[[2]] %in% vois$side2$extra_voi_names) )
      return (mis[ x[[2]], x[[1]] ])
    if ( ! (x[[1]] %in% vois$side1$extra_voi_names) )
      return (mis[ x[[1]], x[[2]] ])
    #
    # 2 variables given as metadata => the MI has not been computed
    #
    list_vois_vals <- list (
      "voi1"=vois$side1$var_of_interest_values[ , x[[1]] ],
      "voi2"=vois$side2$var_of_interest_values[ , x[[2]] ])
    are_continuous <- unlist (lapply (list_vois_vals, FUN=function(y) {
      return (  is.numeric (y)
        && (length (unique (y[!is.na(y)]) ) >= MIIC_CONTINUOUS_TRESHOLD) )
      } ) )
    completes_samples <- ( (!is.na (list_vois_vals[[1]]))
                         & (!is.na (list_vois_vals[[2]])) )
    ret <- computeMutualInfo (list_vois_vals[[1]][completes_samples],
                              list_vois_vals[[2]][completes_samples],
                              is_continuous=are_continuous, plot=F)
    if (unit == "bits")
      {
      nb_completes_samples <- sum ( (!is.na (list_vois_vals[[1]]))
                                  & (!is.na (list_vois_vals[[2]])) )
      mi_val <- ifelse (corrected,
                        (ret$infok / nb_completes_samples) / LN_2,
                        (ret$info  / nb_completes_samples) / LN_2)
      }
    else
      mi_val <- ifelse (corrected, ret$infok, ret$info)
    return (mi_val)
    } ) )

  couples_mi_0_test <- (couples$mi <= 0)
  if ( any(couples_mi_0_test) )
    {
    couples_mi_0 <- couples[couples_mi_0_test, , drop=F]
    miic_warning ("path feature selection", "MI = 0 for ",
      paste ( apply (couples_mi_0, MARGIN=1, FUN=function(x) {
                paste0 (x[[1]], "-", x[[2]]) } ), collapse=", "),
      ", no feature selection possible on these couple(s)." )
    }
  return (couples)
  }

#-------------------------------------------------------------------------------
# sfp_recurs
#-------------------------------------------------------------------------------
# Recursion for selectFeaturesPath
#
# Principe:
# At depth 1, the variables of interest (vois) are the ones from the users.
# For each couple of vois, sfp_recurs look for contributors
# (selection first on dpi, then on ni3 if method is "score")
# From these contributors, the n_selected most probable are kept as features,
# then these n_selected features become vois at depth + 1.
# e.g. at depth 1, we have 2 calls per couple of vois:
#  - one vois from user for side 1 - features selected => depth 2
#  - features selected - one vois from user for side 2 => depth 2
# and so on at depth 2, 3, ..., depth_max or no feature can be selected
#
# Specific params (see selectFeaturesPath for common parameters)
# - depth: an integer >= 1, current depth of the recursion
# - progress: a real >= 0, progress achieved before the call
# - progress_inc: a real <= 100, the progression achievable by this call:
#   set by the caller at 100 / 2 ^ depth_caller / number of couples
# - plot: a data frame to store features selected with their plot position
# - plot_start: a real between 0 and 100, the beginning of the plot area
# - plot_end: a real between 0 and 100, the end of the plot area
#   features selected will be plotted at the middle of [plot_start, plot_end]
# - all_couples: the list of couples of vois evaluated
# - all_scores: the list of scores used to select the features
#
# Return: a list with 5 items. The list are completed as the recursion move
# forward.
# - "features": a vectot with the features selected
# - "mis": the MI matrix
# - "couples": the couples evaluated (x, y, mi, features)
# - "scores": the variables evaluated (x, y, z, mis)
# - "plot": if plotting is requested, a data frame with the plot information
#-------------------------------------------------------------------------------
sfp_recurs <- function (input_data,
  var_of_interest_names_side1, var_of_interest_values_side1,
  var_of_interest_names_side2, var_of_interest_values_side2,
  method, n_selected, corrected, precomputed_mis, skip_cheks,
  n_threads, verbose, depth_max, depth, progress, progress_inc,
  plot, plot_start, plot_end, all_couples, all_scores)
  {
  if (depth > depth_max)
    return (list ("mis"=precomputed_mis, "couples"=all_couples,
                  "scores"=all_scores, "plot"=plot) )
  vois <- sfp_check_vois (input_data=input_data,
    var_of_interest_names_side1=var_of_interest_names_side1,
    var_of_interest_values_side1=var_of_interest_values_side1,
    var_of_interest_names_side2=var_of_interest_names_side2,
    var_of_interest_values_side2=var_of_interest_values_side2)

  if (verbose >= 2)
    {
    str_disp1 <- paste (vois$side1$all_voi_names, collapse=",")
    if (nchar (str_disp1) > 30)
      str_disp1 <- paste0 (substr(str_disp1, 1, 27), "...")
    str_disp2 <- paste (vois$side2$all_voi_names, collapse=",")
    if (nchar (str_disp2) > 30)
      str_disp2 <- paste0 (substr(str_disp2, 1, 27), "...")
    str_display <- paste0 ("Depth ", depth, ", ", str_disp1, "-", str_disp2)
    if (depth == 1)
      miic_msg (str_display, ", computing MIs...")
    else if ( (progress == 0) && (verbose >= 3) )
      cat (paste0 (str_display, ", progress ", round (progress, 2), " %...",
                   # TODO other way to get 40 spaces ?
                   paste (rep(" ", 40), collapse=""), "\r") )
    }
  #
  # Compute MIs
  #
  mat_mis <- compute_mi_batch (input_data=input_data,
    var_of_interest_names=vois[["all"]]$var_of_interest_names,
    var_of_interest_values=vois[["all"]]$var_of_interest_values,
    unit="log_conf", corrected=corrected, precomputed_mis=precomputed_mis,
    skip_cheks=skip_cheks, n_threads=n_threads,
    verbose=ifelse (depth==1, min(verbose, 2), 0) )
  #
  # The mat_mis can contain more rows than the features to select
  # e.g. we pre-computed the MI for all genes and now we want select only
  # the TFs. In this case, in input_data, the variables are only the TFs
  # while mat_mis would contain all genes. Same for the columns as
  # we can have precomputed more variables of interest than the ones we use now
  #
  mat_mis_filt <- mat_mis[rownames(mat_mis) %in% colnames(input_data),
                          colnames(mat_mis) %in% vois[["all"]]$all_voi_names,
                          drop=F]
  #
  # Prepare the couples of vois (one from each side),
  # For each, pick or compute the MI between the vois
  #
  if ( (verbose >= 2) && (depth == 1) )
    miic_msg (str_display, ", evaluate couples...")
  couples <- sfp_prepare_couples (
    vois=vois, mis=mat_mis_filt, unit="log_conf", corrected=corrected)
  couples$depth <- depth
  couples$features <- NA_character_
  couples <- couples[, colnames(all_couples), drop=F]
  #
  # Filter out couple already done, update all_couples
  #
  couples_done <- apply (couples, MARGIN=1, FUN=function(one_row) {
    any ( (all_couples$x == one_row[["x"]])
        & (all_couples$y == one_row[["y"]]) )
    } )
  couples <- couples[ ! couples_done, , drop=F]
  all_couples <- rbind (all_couples, couples)
  #
  # For the next steps, keep only couples with MI > 0
  #
  couples <- couples[couples$mi > 0, , drop=F]
  if (nrow (couples) <= 0)
    return (list ("mis"=mat_mis, "couples"=all_couples,
                  "scores"=all_scores, "plot"=plot) )
  #
  # For each voi requested, filter variables on MI > 0
  #
  list_mis_sup_0 <- list()
  for (one_voi in vois[["all"]]$all_voi_names)
    {
    mis_tmp <- mat_mis_filt[, one_voi]
    mis_tmp <- mis_tmp[ ( ! is.na (mis_tmp) ) & (mis_tmp > 0) ]
    list_mis_sup_0[[one_voi]] <- mis_tmp
    }
  #
  # For each couple, filter variables passsing the basic DPI test:
  # for a pair xy, all z so that MI xz >= MI xy and MI zy >= MI xy
  # Then compute a DPI value and, if method is "score", the NI3 and score
  #
  list_scores_sorted <- list()
  for ( i in 1:nrow (couples) )
    {
    x_name <- couples[i, "x"]
    y_name <- couples[i, "y"]
    if (  (verbose >= 2) && (depth == 1) )
      miic_msg ("Depth ", depth, ", ", x_name, "-", y_name,
                       ", computing DPIs...")
    #
    # Basic DPI check : keep only features z when Ixz > Ixy and Iyz > Ixy
    #
    feat_voi1 <- list_mis_sup_0[[ x_name ]]
    feat_voi2 <- list_mis_sup_0[[ y_name ]]
    mi_threshold <- couples[i, "mi"]
    feat_voi1 <- feat_voi1[feat_voi1 >= mi_threshold]
    feat_voi2 <- feat_voi2[feat_voi2 >= mi_threshold]

    feat_kept <- names(feat_voi1) [names(feat_voi1) %in% names(feat_voi2)]
    feat_kept <- feat_kept[ ! (feat_kept %in% c(x_name, y_name) ) ]

    df_scores <- all_scores[FALSE, , drop=F]
    if (length (feat_kept) <= 0)
      {
      list_scores_sorted[[i]] <- df_scores
      next
      }
    df_scores[feat_kept, "z"] <- feat_kept
    rownames(df_scores) <- feat_kept
    df_scores[, c("depth", "x", "y", "mi_xy")] <- list (
      depth, x_name, y_name, couples[i, "mi"])
    #
    # Store values explaining DPI test
    #
    df_scores [feat_kept, "mi_xz"] <- mat_mis_filt[feat_kept, x_name]
    df_scores [feat_kept, "mi_zy"] <- mat_mis_filt[feat_kept, y_name]
    xz <- df_scores$mi_xz - df_scores$mi_xy
    zy <- df_scores$mi_zy - df_scores$mi_xy
    df_scores$dpi <- ifelse (xz < zy, xz, zy)
    # TODO: check why this term ? To be comparable with NI3 ?
    df_scores$dpi <- df_scores$dpi - log1p ( exp ( -abs (xz - zy) ) )
    #
    # Refine with ni3 if method is "score"
    #
    if (method == "score")
      {
      if (  (verbose >= 2) && (depth == 1) )
        miic_msg ("Depth ", depth, ", ", x_name, "-", y_name,
                  ", computing scores...")
      if (x_name %in% vois$all$var_of_interest_names)
        x <- input_data[ , x_name]
      else
        x <- vois$all$var_of_interest_values[ , x_name]
      x_continuous <- ( is.numeric (x)
        && (length (unique (x[!is.na(x)]) ) >= MIIC_CONTINUOUS_TRESHOLD) )

      if (y_name %in% vois$all$var_of_interest_names)
        y <- input_data[ , y_name]
      else
        y <- vois$all$var_of_interest_values[ , y_name]
      y_continuous <- ( is.numeric (y)
        && (length (unique (y[!is.na(y)]) ) >= MIIC_CONTINUOUS_TRESHOLD) )
      #
      # Compute NI3 for all remaining features
      #
      df_scores [, "i3"] <- apply (df_scores, MARGIN=1, FUN=function(one_row) {
        z_name <- one_row[["z"]]
        z <- input_data[ , z_name, drop=T]
        #
        # Filter rows with 1 NA as computeThreePointInfo needs complete samples
        #
        df_tmp <- data.frame ("x"=x, "y"=y, "z"=z, stringsAsFactors=F)
        has_na <- apply (df_tmp, MARGIN=1, FUN=function(x) { anyNA(x) } )
        df_tmp <- df_tmp[ !has_na, , drop=F]
        if (nrow (df_tmp) <= 0)
          return (NA_real_)

        z_continuous <- ( is.numeric (df_tmp$z)
          && (length (unique (df_tmp$z[!is.na(df_tmp$z)]) ) >= MIIC_CONTINUOUS_TRESHOLD) )
        are_continuous <- c (x_continuous, y_continuous, z_continuous)
        ret_ni3 <- computeThreePointInfo (x=df_tmp$x, y=df_tmp$y, z=df_tmp$z,
                                          is_continuous=are_continuous)
        ifelse (corrected, ret_ni3$i3k, ret_ni3$i3)
        } )
      df_scores$score <- ifelse (df_scores$dpi < df_scores$i3,
                                 df_scores$dpi, df_scores$i3)
      }
    #
    # Order and memorize scores/dpi for next round
    #
    df_scores_sorted <- df_scores[ ( ! is.na (df_scores[, method]) )
                                 & (df_scores[, method] > 0), , drop=F]
    df_scores_sorted <- df_scores_sorted[order (df_scores_sorted[, method],
                                                decreasing=T), , drop=F]
    list_scores_sorted[[i]] <- df_scores_sorted
    #
    # Align scores of this round to all_scores to be returned back
    #
    rownames(df_scores) <- NULL
    all_scores <- rbind (all_scores, df_scores)
    }
  #
  # DPIs or/and scores computed, select the variables that will be kept:
  # the most probable contributors (ais)
  #
  list_ais_couples <- list()
  for ( i in 1:nrow (couples) )
    {
    if (nrow (list_scores_sorted[[i]]) <= 0)
      {
      list_ais_couples[[i]] <- NA_character_
      next
      }
    x_name <- couples[i, "x"]
    y_name <- couples[i, "y"]
    list_ais <- rownames (list_scores_sorted[[i]])[1:min (
      n_selected, nrow(list_scores_sorted[[i]]) )]
    all_couples [ (all_couples$depth == depth)
                & (all_couples$x == x_name)
                & (all_couples$y == y_name),
                "features" ] <- paste (list_ais, collapse=",")
    list_ais_couples[[i]] <- list_ais
    }
  #
  # Recursive calls with variables found to be in the path as new vois
  #
  # If we have only 1 couple of vois to evaluate at each level of recursion,
  # the progress increment would be divided by 2 at each level
  # as we have 2 sides of the dichotomy to investigate
  # at depth 1, 100 / 2 => 50
  # - recursion on left done = 50%, on right = 100%
  # at depth 2,  50 / 2 => 25
  # - if progress was  0, recursion done on left = 25% and on right = 50%
  # - if progress was 50, recursion done on left = 75% and on right = 100%
  # As we can have several couples of variables, we divide also
  # by the number of couples.
  #
  progress_inc <- (progress_inc / 2) / nrow (couples)
  plot_middle <- (plot_start + plot_end) / 2
  if ( (depth == 1) && (verbose == 2) )
    miic_msg (str_display, ", recursing...")
  features <- c()
  for ( i in 1:nrow (couples) )
    {
    list_ais <- list_ais_couples[[i]]
    list_ais <- list_ais[ !is.na (list_ais) ]
    if (length (list_ais) <= 0)
      {
      progress <- progress + progress_inc * 2
      next
      }
    if ( ! is.null(plot) )
      {
      # Add the most probable contributors as the features selected,
      # Plot position is the middle of the recursion step:
      # at depth 1: 50, at depth 2: 25 or 75, ...
      #
      idxs <- ( (nrow(plot)+1):(nrow(plot)+length(list_ais)) )
      plot[idxs, ] <- list (plot_middle, depth, list_ais)
      }
    features <- unique (c (features, list_ais) )
    #
    # Recurs with 1st voi and selected features (= the most probable contributors)
    #
    x_name <- couples[i, "x"]
    var_of_interest_names_side1 <- NULL
    var_of_interest_values_side1 <- NULL
    if (x_name %in% vois[["side1"]]$var_of_interest_names)
      var_of_interest_names_side1 <- x_name
    else
      var_of_interest_values_side1 <- vois[["side1"]]$var_of_interest_values[
        , x_name, drop=F]

    input_data_rec <- input_data[,
      unique ( c (var_of_interest_names_side1,
                  rownames(list_scores_sorted[[i]])) ),
      drop=F]

    ret <- sfp_recurs (input_data=input_data_rec,
      var_of_interest_names_side1=var_of_interest_names_side1,
      var_of_interest_values_side1=var_of_interest_values_side1,
      var_of_interest_names_side2=list_ais,
      var_of_interest_values_side2=NULL,
      method=method, n_selected=n_selected, corrected=corrected,
      precomputed_mis=mat_mis, skip_cheks=skip_cheks,
      n_threads=n_threads, verbose=verbose,
      depth_max=depth_max, depth=depth+1,
      progress=progress, progress_inc=progress_inc,
      plot=plot, plot_start=plot_start, plot_end=plot_middle,
      all_couples=all_couples, all_scores=all_scores)

    features <- unique ( c (features, ret$features) )
    mat_mis <- ret$mis
    all_couples <- ret$couples
    all_scores <- ret$scores
    plot <- ret$plot

    progress <- progress + progress_inc
    if (verbose >= 3)
      cat (paste0 ("Depth ", depth, ", ", x_name, "-", y_name,
        ", progress ", round (progress, 2), " %...",
        paste (rep(" ", 50), collapse=""), "\r") )
    #
    # Recurs with selected features (= the most probable ais)
    # and voi from the other side
    #
    y_name <- couples[i, "y"]
    var_of_interest_names_side2 <- NULL
    var_of_interest_values_side2 <- NULL
    if (y_name %in% vois[["side2"]]$var_of_interest_names)
      var_of_interest_names_side2 <- y_name
    else
      var_of_interest_values_side2 <- vois[["side2"]]$var_of_interest_values[,
        y_name, drop=F]

    input_data_rec <- input_data[,
      unique ( c (var_of_interest_names_side2,
                  rownames(list_scores_sorted[[i]])) ),
      drop=F]
    ret <- sfp_recurs (input_data=input_data_rec,
      var_of_interest_names_side1=list_ais,
      var_of_interest_values_side1=NULL,
      var_of_interest_names_side2=var_of_interest_names_side2,
      var_of_interest_values_side2=var_of_interest_values_side2,
      method=method, n_selected=n_selected, corrected=corrected,
      precomputed_mis=mat_mis, skip_cheks=skip_cheks,
      n_threads=n_threads, verbose=verbose,
      depth_max=depth_max, depth=depth+1,
      progress=progress, progress_inc=progress_inc,
      plot=plot, plot_start=plot_middle, plot_end=plot_end,
      all_couples=all_couples, all_scores=all_scores)

    features <- unique ( c( features, ret$features) )
    mat_mis <- ret$mis
    all_couples <- ret$couples
    all_scores <- ret$scores
    plot <- ret$plot

    progress <- progress + progress_inc
    if (verbose >= 3)
      cat (paste0 ("Depth ", depth, ", ", x_name, "-", y_name,
        ", progress ", round (progress, 2), " %...",
        paste (rep(" ", 50), collapse=""), "\r") )
    }
  if ( (depth == 1) && (verbose >= 3) )
    miic_msg (str_display, ", progress 100 %",
              paste (rep(" ", 50), collapse="") )

  return (list ("features"=features, "mis"=mat_mis,
                "couples"=all_couples, "scores"=all_scores,
                "plot"=plot) )
  }

#-------------------------------------------------------------------------------
# sfp_plot
#-------------------------------------------------------------------------------
# Params (internal):
# - vois: a list, the vois returned by sfp_check_vois
# - df_plots: a data frame with the features and their position,
#   returned by sfp_recurs
# Params (can be provided by the user via the ... in selectFeaturesPath):
# - depth_plot: an integer in the range [1:10], default 4,
#   the maximum depth plotted
# - annotate: a boolean, default T, indicates the depth of the recursion
# - x_lab: the x label, default "Features for vois_side1-vois_side2"
# - font_size: an integer >= 1, default 11, the font size
#-------------------------------------------------------------------------------
sfp_plot <- function (vois, df_plots, depth_plot=4, annotate=T,
                      x_lab=NULL, font_size=11)
  {
  depth_plot <- check_param_int (depth_plot, "plot depth", default=4, min=1, max=10)
  annotate <- check_param_logical (annotate, "plotting of annotation", default=T)
  if ( is.null (x_lab) )
    x_lab <- paste0 ("Features for ",
                     paste (vois$side1$all_voi_names, collapse=", "), " - ",
                     paste (vois$side2$all_voi_names, collapse=", ") )
  else
    x_lab <- as.character (x_lab)
  font_size <- check_param_int (font_size, "font size", default=11, min=1)

  depth_max <- ifelse ( nrow(df_plots) <= 0, 0, max(df_plots$depth) )
  all_feats <- unique (df_plots$features)
  df_plots <- df_plots[df_plots$depth <= depth_plot, , drop=F]
  df_plots$depth <- NULL
  df_plots <- unique (df_plots)

  y_feat <- 0.5

  side1_names <- paste (sort (vois$side1$all_voi_names), collapse="\n")
  side2_names <- paste (sort (vois$side2$all_voi_names), collapse="\n")
  g <- ggplot2::ggplot() +
    ggplot2::theme_classic() +
    ggplot2::theme ( text=ggplot2::element_text (size=font_size),
      axis.text.x=ggplot2::element_text (angle=30, vjust=1, hjust=1),
      axis.line.y=ggplot2::element_blank(),
      axis.text.y=ggplot2::element_blank(),
      axis.ticks.y=ggplot2::element_blank() ) +
    ggplot2::labs (x=x_lab, y="") +
    ggplot2::ylim (0, 1) +
    ggplot2::geom_text (ggplot2::aes (
        x=c(0,1), y=y_feat, label=c(side1_names, side2_names) ),
      size=font_size*0.8/ggplot2::.pt, hjust=0.5, vjust=0.5, color="black") +
    ggplot2::geom_vline (xintercept=c(0,1), linewidth=0.5,
                         linetype="dashed", color="darkgrey")
  xticks_text <- c ("0", "1")
  xticks_pos <- c (0, 1)
  #
  # Pick the features position at the highest depth
  #
  plot_positions <- c()
  plot_texts <- c()
  feats_displayed <- c()
  for (i in 1:depth_plot)
    {
    pos_inc <- 100 / (2 ^ i)
    for ( one_pos in seq (from=pos_inc, to=100-pos_inc, by=pos_inc) )
      {
      if (! (one_pos %in% df_plots$pos) )
        next
      if (one_pos %in% plot_positions)
        next
      feats_to_add <- unique (df_plots$features[df_plots$pos == one_pos])
      feats_to_add <- feats_to_add[ ! (feats_to_add %in% feats_displayed)]
      if (length(feats_to_add) <= 0)
        next
      plot_positions <- c (plot_positions, one_pos)
      plot_texts <- c (plot_texts, paste (sort (feats_to_add), collapse="\n") )
      feats_displayed <- c(feats_displayed, feats_to_add)
      }
    }
  plot_positions <- round (plot_positions / 100, 3)
  g <- g +
    ggplot2::geom_text (
      ggplot2::aes (x=plot_positions, y=y_feat, label=plot_texts),
      size=font_size*0.8/ggplot2::.pt, hjust=0.5, vjust=0.5, color="black") +
    ggplot2::geom_vline ( xintercept=plot_positions, linewidth=0.6,
                          linetype="dotted", color="darkgrey")
  xticks_text <- c (xticks_text, as.character(plot_positions) )
  xticks_pos <- c (xticks_pos, plot_positions)
  #
  # Set ticks
  #
  g <- g + ggplot2::scale_x_continuous (breaks=xticks_pos, label=xticks_text)
  #
  # Add annotation if plot depth > depth search
  #
  if (annotate)
    {
    if (depth_max == 0)
      label_txt <- paste0 ("No feature found")
    else
      {
      if (depth_max > depth_plot)
        label_txt <- paste0 (depth_plot,  " levels shown of ", depth_max)
      else
        label_txt <- paste0 ("All ", depth_max,  " levels shown")

      shown_feats <- unique (df_plots$features)
      n_feat_not_shown <- sum ( ! (all_feats %in% shown_feats) )
      if (n_feat_not_shown <= 0)
        label_txt <- paste0 (label_txt, ",\nall ", length(shown_feats),
                             " features shown")
      else
        label_txt <- paste0 (label_txt, ",\n", n_feat_not_shown,
                             " feature(s) not shown")
      }
    g <- g + ggplot2::geom_text (ggplot2::aes (x=0.875, y=1, label=label_txt),
      size=font_size*0.8/ggplot2::.pt, hjust=0.5, vjust=1)
    }
  return (g)
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
  vois <- sfo_check_vois (input_data,
    var_of_interest_names, var_of_interest_values)
  unit <- check_param_string (unit, "unit", c("log_conf", "bits") )
  corrected <- check_param_logical (corrected, "corrected", T)
  precompured_mis <- sf_check_precomputed_mis (precomputed_mis)
  skip_cheks <- check_param_logical (skip_cheks, "skip checks", F)
  n_threads <- check_param_int (n_threads, "number of threads", 1, min=1)
  verbose <- check_param_int (verbose, "verbose", 3, min=0, max=3)
  plot <- check_param_logical (plot, "plot", F)
  #
  # Cross checks
  #
  if ( n_features > ncol(input_data) )
    {
    n_features = ncol(input_data)
    miic_warning ("parameter",  "the number of features can not be greater",
      " than the number of variables in input data.",
      " It has been reduced to ", n_features, ".")
    }
  if ( (n_features > 0) && ( n_features < length (vois$all_voi_names) ) )
    {
    n_features = length (vois$all_voi_names)
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
        list_plots[[one_col]] <- sfo_plot (
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
  list_tops <- sfo_get_tops (
    n_features=n_features, list_sorted=list_mis_sorted, verbose=verbose)

  return (list ("features"=list_tops, "mis"=mat_mis, "plots"=list_plots) )
  }

#-------------------------------------------------------------------------------
# selectFeaturesPath
#-------------------------------------------------------------------------------
#' Features selection on the path between variables of interest
#'
#' @description Select the variables that are the most likely to be in the path
#' between two set of variable(s) of interest.
#' Variables selection can be performed using a score combining Data Processing
#' inequality (DPI) and 3 points information or only using the DPI.\cr
#' The search is recursive by dichotomy: in the first round, a set of features
#' on the path between the main variables of interest is selected.
#' From then, the path is split in two parts:
#' \itemize{
#'   \item variables of interest of side 1 + set of features selected
#'   \item set of features selected + variables of interest of side 2
#'   }
#' Each part of the path is investigated with the same principle:
#' finding a set of features in the path subsection, splitting the path in two,
#' and so on.
#' The recursion ends when no feature can be found
#' or when the maximum level of recursion is reached.
#'
#' @references
#' \itemize{
#' \item Affeldt \emph{et al.}, UAI 2015, \href{https://auai.org/uai2015/proceedings/papers/293.pdf}{Robust Reconstruction of Causal Graphical Models based on Conditional 2-point and 3-point Information}
#' }
#'
#' @param input_data [a data frame or a matrix, required]
#'
#' Expected layout is samples as rows and variables as columns. Column names
#' must contain the names of the variables.
#'
#' @param var_of_interest_names_side1 [a string or vector of strings, optional,
#' NULL by default]
#'
#' For the variable(s) of interest that are part of the \emph{input_data}
#' on the first side of the path, you should supply their names here.
#'
#' @param var_of_interest_values_side1 [a data frame, optional, NULL by default]
#'
#' For the variables of interest that are not in \emph{input_data}
#' and on the first side of the path, a data frame can be supplied.
#' The column names are the names of the variables of interest
#' and rows are the samples ordered in the same way as the \emph{input_data}.
#' Typically, such variables are metadata associated to samples but not
#' stored in \emph{input_data}, e.g. a "Treatment" vs "Control" variable
#' in an experiment and a count matrix with the expression of genes
#' in \emph{input_data}.
#'
#' @param var_of_interest_names_side2 [a string or vector of strings, optional,
#' NULL by default]
#'
#' Same as \emph{var_of_interest_names_side1} for the second side of the path.
#'
#' @param var_of_interest_values_side2 [a data frame, optional, NULL by default]
#'
#' Same as \emph{var_of_interest_values_side1} for the second side of the path.
#'
#' @param method [a string, optional, "score" by default,
#' possible values: "score", "dpi"]
#'
#' When set to "score", the variables selection is performed using a score
#' combining DPI and 3 points information (see Affeldt 2015).
#' By turning it to "dpi", the selection is based only on the DPI,
#' which speeds up the process but is less discriminating.
#'
#' @param n_selected [a positive integer, optional, 10 by default]
#'
#' The number of features selected at each step of the recursion.
#' Decreasing this value speeds up the process while reducing the number
#' of features selected. Increasing it has opposite effect, more features
#' at the cost of an increased processing time.
#'
#' @param corrected [a boolean, optional, TRUE by default]
#'
#' When set to TRUE, the mutual information and 3 points information
#' are corrected by subtracting a complexity term
#' (computed with the Normalized Maximum Likelihood).
#' For dataset having very few samples, the complexity term can have
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
#' To be valid, the pre-computed MI values must have been computed using
#' the same \emph{corrected} parameter.
#'
#' @param skip_cheks [a boolean, optional, FALSE by default]
#'
#' Before computing MI between the variable of interest and the features,
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
#' @param depth_max [a positive integer, optional, 10 by default]
#'
#' The maximum depth of the recursion.
#'
#' @param verbose [an integer, optional, 3 by default]
#'
#' Level of verbosity: 0=no display, 1=summary, 2=progress per couple of
#' variable of interest, 3=include more detail on the progress during
#' the recursion.
#'
#' @param plot [a boolean, optional, FALSE by default]
#'
#' If set to TRUE, a plot with the top features is generated
#' (requires `ggplot2`).
#'
#' Please note on the plot rendering that positions are indicative with
#' a tendency to be displayed to the left side: a feature can appear multiple
#' times during the recursion, but couples already investigated are skipped
#' to speed up the process. For the plotting, the rule applied to select the
#' position of each feature is (from the positions returned by the recursion),
#' to pick the one corresponding to the minimal depth.
#'
#' @param ...
#'
#' If plotting is requested, extra parameters can be used to customize the plot
#' rendering:
#'
#' \itemize{
#' \item \emph{x_lab}: a string, optional, "Features for
#'   variables_of_interest_side1-variables_of_interest_side2" by default,
#'   the X axis label.
#' \item \emph{depth_plot}: an integer between 1 and 10, optional, 4 by default,
#'   maximal depth used for the plot.
#' \item \emph{annotate}: a boolean, optional, TRUE by default,
#'   When activated, displays the maximal depth of the recursion
#'   and the number of features not displayed on the plot.
#' \item \emph{font_size}: an integer, optional, 11 by default, the font size.
#' }
#'
#' @return A named list with five items:
#'
#' \itemize{
#' \item \emph{features}: a vector with the features selected.
#' \item \emph{mis}: a matrix with the MIs between the variables in
#'   \emph{input_data} (as rows) and the variable(s) of interest (as columns).
#'   Row and column names are sorted alphabetically.\cr
#'   The MIs can be corrected or not, depending on the \emph{corrected}
#'   parameters.
#'   If a pre-computed MIs matrix was supplied, new values computed
#'   are added to the existing matrix.
#' \item \emph{couples}: a data frames with, at each depth, the couples
#'   of variables used with their MIs and features selected.
#'   Please note that couples already encountered are not included
#'   as no recursion is performed on duplicates.
#' \item \emph{scores}: a data frame with the information about all pairs of
#'   variables passing the dpi test, containing the depth, the MI values,
#'   the DPI test and, if \emph{method} is "score", the 3 points information
#'   and score.
#' \item \emph{plot}: when the \emph{plot} parameter is turned to TRUE,
#'   a plot with the features of the top levels of recursion, NULL otherwise.
#' }
#'
#' @examples
#' library(miic)
#'
#' \donttest{
#' # Features selection on the path between an external metadata "Ploidy"
#' # (simulation by extracting "Ploidy" out of the dataset)
#' # and the gene expression of TP53
#' df_external_meta <- data.frame ("Ploidy"=cosmicCancer$Ploidy)
#' df_data <- cosmicCancer[ , ! (colnames(cosmicCancer) == "Ploidy") ]
#' ret <- selectFeaturesPath (df_data,
#'                            var_of_interest_values_side1=df_external_meta,
#'                            var_of_interest_names_side2="TP53")
#' message ("Features selected: ", paste (ret$features, collapse=", ") )
#'
#' # Same features selection with reuse of the MIs computed above and plot
#' ret <- selectFeaturesPath (df_data,
#'                            var_of_interest_values_side1=df_external_meta,
#'                            var_of_interest_names_side2="TP53",
#'                            precomputed_mis=ret$mis,
#'                            plot=TRUE)
#' print (ret$plot)
#'
#' # Features selection using multiple variables on each side
#' # (reusing the MIs computed above and selecting only one feature per round
#' # to speed up the example)
#' ret <- selectFeaturesPath (cosmicCancer,
#'                            var_of_interest_names_side1=c("TP53", "Ploidy"),
#'                            var_of_interest_names_side2=c("FOXM1", "AURKA"),
#'                            n_selected=1, precomputed_mis=ret$mis)
#' message ("Features selected: ", paste (ret$features, collapse=", ") )
#'
#' # Similar features selection with plotting using a customized rendering
#' # (reusing the MIs computed above and limiting the recursion depth
#' # to speed up the example)
#' ret <- selectFeaturesPath (cosmicCancer,
#'                            var_of_interest_names_side1=c("TP53", "Ploidy"),
#'                            var_of_interest_names_side2=c("FOXM1", "AURKA"),
#'                            depth_max=2, plot=TRUE, font_size=12,
#'                            x_lab="My features selection on CosmicCancer")
#' print (ret$plot)
#' }
#'
#' @export
#-------------------------------------------------------------------------------
selectFeaturesPath <- function (input_data,
  var_of_interest_names_side1=NULL, var_of_interest_values_side1=NULL,
  var_of_interest_names_side2=NULL, var_of_interest_values_side2=NULL,
  method="score", n_selected=10, corrected=T, precomputed_mis=NULL,
  skip_cheks=F, n_threads=1, depth_max=10, verbose=3, plot=F, ...)
  {
  # Check parameters
  #
  input_data <- sf_check_input_data (input_data=input_data)
  vois <- sfp_check_vois (input_data=input_data,
    var_of_interest_names_side1=var_of_interest_names_side1,
    var_of_interest_values_side1=var_of_interest_values_side1,
    var_of_interest_names_side2=var_of_interest_names_side2,
    var_of_interest_values_side2=var_of_interest_values_side2)
  method <- check_param_string ( method, "method", c("score", "dpi") )
  n_selected <- check_param_int (n_selected,
    "number of selected feature per round", default=10, min=1)
  corrected <- check_param_logical (corrected, "corrected", T)
  precomputed_mis <- sf_check_precomputed_mis (precomputed_mis)
  skip_cheks <- check_param_logical (skip_cheks, "skip checks", F)
  n_threads <- check_param_int (
    n_threads, "number of threads", default=1, min=1)
  depth_max <- check_param_int (
    depth_max, "maximum depth", default=10, min=1)
  verbose <- check_param_int (verbose, "verbose", 3, 0, 3)
  plot <- check_param_logical (plot, "plot", F)

  all_couples <- data.frame (
    "depth"=integer(), "x"=character(), "y"=character(),
    "mi"=numeric(), "features"=character(), stringsAsFactors=F)
  all_scores <- data.frame (
    "depth"=integer(), "x"=character(), "y"=character(), "z"=character(),
    "mi_xy"=numeric(), "mi_xz"=numeric(), "mi_zy"=numeric(), "i3"=numeric(),
    "dpi"=numeric(),  "score"=numeric(), stringsAsFactors=F)
  if (plot)
    df_plots <- data.frame ("pos"=numeric(), "depth"=integer(), "features"=character(),
                            stringsAsFactors=F)
  else
    df_plots <- NULL

  if (verbose >= 1)
    miic_msg ("Selecting features on path between ",
      paste0 (vois$side1$all_voi_names, collapse=","),
      " and ", paste0 (vois$side2$all_voi_names, collapse=","), "...")

  ret_recurs <- sfp_recurs (input_data=input_data,
    var_of_interest_names_side1 =vois[["side1"]]$var_of_interest_names,
    var_of_interest_values_side1=vois[["side1"]]$var_of_interest_values,
    var_of_interest_names_side2 =vois[["side2"]]$var_of_interest_names,
    var_of_interest_values_side2=vois[["side2"]]$var_of_interest_values,
    method=method, n_selected=n_selected, corrected=corrected,
    precomputed_mis=precomputed_mis, skip_cheks=skip_cheks,
    n_threads=n_threads, verbose=verbose,
    depth_max=depth_max, depth=1, progress=0, progress_inc=100,
    plot=df_plots, plot_start=0, plot_end=100,
    all_couples=all_couples, all_scores=all_scores)

  if (plot)
    {
    if ( base::requireNamespace("ggplot2", quietly=TRUE) )
      ret_recurs$plot <- sfp_plot (vois, unique (ret_recurs$plot), ...)
    else
      miic_warning ("Path features selection", "Plotting requires ggplot2.")
    }

  if (verbose >= 1)
    miic_msg (length (ret_recurs$features), " features selected.")
  return (ret_recurs)
  }
