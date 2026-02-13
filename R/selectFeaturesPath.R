#*******************************************************************************
# Filename   : selectFeaturesPath.R             Creation date: 10 February 2026
#
# Description: Features selection on the path between two sets of variables of
#              interest
#
# Author     : Franck SIMON
#
# TODO ? warning discrete number of levels (as miic) ?")
#*******************************************************************************

#===============================================================================
# INTERNAL FUNCTIONS FOR selectFeaturesPath (sfp_xx)
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
#   the same number of rows than input_data
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
            " names for side ", i, " (e.g. ", one_var_name,
            ") are incorrect or not in the input_data.")
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

    vois[[paste0 ("side", i)]] <- list (
      "var_of_interest_names" = var_of_interest_names,
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

  vois[["all"]] <- list (
    "var_of_interest_names" = var_of_interest_names_all_sides,
    "var_of_interest_values" = var_of_interest_values_all_sides,
    "extra_voi_names" = extra_voi_names_all_sides,
    "all_voi_names" = all_voi_names)
  return (vois)
  }

#-------------------------------------------------------------------------------
# sfp_check_precomputed_values
#-------------------------------------------------------------------------------
# Check values coming from a previous run
# Params:
# - precomputed_values: a list containing the precomputed MIs and scores
# Return:
# - checked precomputed_values list
#-------------------------------------------------------------------------------
sfp_check_precomputed_values <- function (precomputed_values)
  {
  # Basic checks of precomputed_values
  #
  if ( is.null (precomputed_values) )
    return (NULL)

  if ( ! is.list (precomputed_values) )
    {
    miic_warning ("parameters",
      "Pre-computed values must be provided as a list. They will be ignored.")
    return (NULL)
    }
  item_names <- names (precomputed_values)
  #
  # Check mis
  #
  if ( ! ("mis" %in% item_names) )
    {
    miic_warning ("parameters", "Pre-computed values must contain a 'mis' item.",
      " They will be ignored.")
    return (NULL)
    }
  precomputed_mis <- precomputed_values[["mis"]]
  if ( ! is.list (precomputed_mis) )
    {
    miic_warning ("parameters",
      "Pre-computed values 'mis' item must be a list of numeric matrices.",
      " The pre-computed values will be ignored.")
    return (NULL)
    }
  #
  # Utility function for mis matrices check
  #
  check_one_mis_mat <- function (mat)
    {
    is.matrix (mat) && ncol (mat) >= 2 && is.numeric (mat)
    }
  if ( ! all (vapply (precomputed_mis, check_one_mis_mat, logical(1) ) ) )
    {
    miic_warning ("parameters",
      "Pre-computed values 'mis' item contains invalid matrix/matrices.",
      " The pre-computed values will be ignored.")
    return (NULL)
    }
  #
  # Check scores
  #
  if ( ! ("scores" %in% item_names) )
    {
    miic_warning ("parameters",
      "Pre-computed values must contain a 'scores' item.",
      " They will be ignored.")
    return (NULL)
    }
  precomputed_scores <- precomputed_values[["scores"]]
  if ( ! is.list (precomputed_scores) )
    {
    miic_warning ("parameters",
      "Pre-computed values 'scores' item must be a list of data frames.",
      " The pre-computed values will be ignored.")
    return (NULL)
    }
  #
  # Utility function for scores data frames check
  #
  check_df_structure <- function (df, expected_df)
    {
    (   is.data.frame(df)
    && length (colnames(df)) == length (colnames(expected_df))
    && all (colnames(df) == colnames (expected_df))
    && identical ( vapply(df, typeof, ""), vapply(expected_df, typeof, "") ) )
    }
  expected_scores <- data.frame (
    "x"=character(), "y"=character(), "z"=character(), "ais"=character(),
    "mi_xy"=numeric(), "mi_xz"=numeric(), "mi_zy"=numeric(),
    "dpi"=numeric(),  "i3"=numeric(), "score"=numeric(), stringsAsFactors=F)

  if ( ! all (vapply (precomputed_scores, check_df_structure, logical(1), expected_scores)))
    {
    miic_warning ("parameters",
      "Pre-computed values 'scores' item contains invalid data frame(s).",
      " The pre-computed values will be ignored.")
    return (NULL)
    }

  return (precomputed_values)
  }

#-------------------------------------------------------------------------------
# sfp_init_couples_deeper
#-------------------------------------------------------------------------------
# Prepare a data frame for the first depth exploration containing
# all the couples between the variables of interest (vois),
# Each row of the data frame is associated with a list of remaining variables
# (all variables of the dataset except vois at the init)
#
# Params:
# - vois: the list returned by sfp_check_vois
# - vars_to_test: a list of variables
# Return:
# - "couples": a data frame with all the couples as rows
# - "vars_remain": a list with the variables (candidate features) for each
#   couple
#-------------------------------------------------------------------------------
sfp_init_couples_deeper <- function (vois1, vois2, vars_to_test)
  {
  couples_deeper <- expand.grid ("voi1"= vois1, "voi2"=vois2, stringsAsFactors=F)
  vars_to_test <- vars_to_test[ ! vars_to_test %in% c(vois1, vois2) ]
  #
  # Always order vois alphabetically (to find couples easily:
  # not need to test x-y and y-x, it is always x-y)
  #
  swap <- couples_deeper$voi1 > couples_deeper$voi2
  couples_deeper[swap, c("voi1", "voi2")] <- couples_deeper[swap, c("voi2", "voi1")]
  #
  # Associate each couple to its variables remaining
  #
  vars_remain <- replicate (nrow (couples_deeper), vars_to_test, simplify=FALSE)

  return ( list ("couples" = couples_deeper, "vars_remain" = vars_remain) )
  }

#-------------------------------------------------------------------------------
# sfp_prepare_couples_deeper
#-------------------------------------------------------------------------------
# Prepare a data frame for the next depth exploration, containing all the
# couples between the variables of interest (vois) and the variables selected
# as features.
# For each couple, the remaining list of variables is filtered
# to keep only variables sharing information with voi1 or voi2:
#     if couple is voi1-feat_x, then keep only the z | I'(voi1;z) > 0
# and if couple is voi2-feat_x, then keep only the z | I'(voi2;z) > 0)
#
# Params:
# - couple: the couple of vois to evaluate
# - features: the list of features selected for this couple of vois
# - mis: a matrix containing the MIs between the variables of
#   interest (as columns) and variables evaluated (as rows)
#   It should contain only pertinent variables for the couple voi1-voi2
# Return: a list containing:
# - "couples": a data frame with all the couples as rows
# - "vars_remain": a list with the remaining variables (candidate features)
#   for each couple
#-------------------------------------------------------------------------------
sfp_prepare_couples_deeper <- function (couple, features, mis)
  {
  # If no feature selected, no couple to explore at next depth
  #
  n_feats <- length (features)
  if (n_feats <= 0)
    return ( list ("couples" = data.frame ("voi1"=character(),
                     "voi2"=character(), stringsAsFactors=F),
                   "vars_remain" = list() ) )
  #
  # Features selected, prepare the list of couples
  #
  couples_deeper <- expand.grid (
    "voi1"= couple$voi1, "voi2"=features, stringsAsFactors=F)
  couples_deeper <- rbind (couples_deeper, expand.grid (
    "voi1"=features, "voi2"=couple$voi2, stringsAsFactors=F) )
  #
  # Always order vois alphabetically (to find couples easily:
  # not need to test x-y and y-x, it is always x-y)
  #
  swap <- couples_deeper$voi1 > couples_deeper$voi2
  couples_deeper[swap, c("voi1", "voi2")] <- couples_deeper[swap, c("voi2", "voi1")]
  #
  # Associate each couple to its variables remaining (having MI > 0 for the voi)
  #
  vars_remain_voi1 <- rownames(mis)[ !is.na (mis[, couple$voi1])
                                   & mis[, couple$voi1] > 0]
  vars_remain_voi2 <- rownames(mis)[ !is.na (mis[, couple$voi2])
                                   & mis[, couple$voi2] > 0]
  vars_remain <- c ( replicate (n_feats, vars_remain_voi1, simplify=FALSE),
                     replicate (n_feats, vars_remain_voi2, simplify=FALSE) )
  #
  # For each couple, remove the vois from the vars remaining.
  # It is true for the vois of the couple:
  # e.g. if the voi2 of the couple is "Col3" and "Col3" is part of input_data
  # => "Col3a1" is in the mis matrix rows
  # => "Col3a1" must be removed from the remaining variables for next depth
  # It is true also for the new vois, the variables selected as features:
  # e.g. "Tcf4" was a remaining variable and has been selected as feature
  # => "Tcf4" is in the mis matrix rows and has MI with the two vois > 0
  # => "Tcf4" must be removed from the remaining variables for next depth
  #
  for ( i in seq_len (nrow (couples_deeper)) )
    vars_remain[[i]] <- vars_remain[[i]] [
      ! (vars_remain[[i]] %in% c( couple [c("voi1", "voi2")],
                                  couples_deeper[i, c("voi1", "voi2")]) ) ]
  if ( nrow (unique (couples_deeper)) != nrow (couples_deeper) )
    {
    print (couples_deeper)
    stop ("couples_deeper not unique !")
    }
  return ( list ("couples" = couples_deeper, "vars_remain" = vars_remain) )
  }

#-------------------------------------------------------------------------------
# sfp_filter_data
#-------------------------------------------------------------------------------
# Utility function to filter data for one couple evaluation
# Params:
# - voi1: the first voi name
# - voi2: the second voi name
# - input_data: the (unfiltered) input data
# - input_names: variable to keep in input_data
# - extra_data: the (unfiltered) extra_data
# - extra_names: variable to keep in extra_data
# - ais_vect: the list of variables used in the condionning set
# Return: a named list with
# - "voi1_vals": the voi1 values
# - "voi2_vals": the voi2 values
# - "input_data": the input_data filtered on input_names
# - "extra_data": the extra_data filtered on extra_names
# - "cond_data": the cond_data filtered on ais_vect
#-------------------------------------------------------------------------------
sfp_filter_data <- function (
  voi1, voi2, input_data, input_names, extra_data, extra_names, ais_vect)
  {
  if (voi1 %in% extra_names)
    voi1_vals <- extra_data[, voi1]
  else
    voi1_vals <- input_data[, voi1]
  if (voi2 %in% extra_names)
    voi2_vals <- extra_data[, voi2]
  else
    voi2_vals <- input_data[, voi2]
  input_data_filt <- input_data[ , input_names, drop=F]
  extra_data_filt <- NULL
  if (length (extra_names) > 0)
    extra_data_filt <- extra_data[, extra_names, drop=F]
  cond_data_filt <- NULL
  if (length (ais_vect) > 0)
    cond_data_filt <- input_data[, ais_vect, drop=F]

  return ( list ("voi1_vals"=voi1_vals, "voi2_vals"=voi2_vals,
                 "input_data"=input_data_filt, "extra_data"=extra_data_filt,
                 "cond_data"=cond_data_filt) )
  }

#-------------------------------------------------------------------------------
# sfp_get_mi
#-------------------------------------------------------------------------------
# Return the MI between two variables
# Params:
# - mis: a matrix with pre-computed MIs
# - x_name: the first variable name
# - x_vals: the first variable values
# - y_name: the second variable name
# - y_vals: the second variable values
# - cond_data: data of conditioning set applied for alternate path exploration
# - corrected: a boolean indicating if we use corrected/non corrected MI
# Return:
# - the MI between the two variables
#-------------------------------------------------------------------------------
sfp_get_mi <- function (mis, x_name, x_vals, y_name, y_vals, cond_data, corrected=T)
  {
  # Get MI from mis matrix if already computed
  #
  if ( (x_name %in% rownames (mis)) && (y_name %in% colnames (mis)) )
    return (mis[x_name, y_name])
  if ( (y_name %in% rownames (mis)) && (x_name %in% colnames (mis)) )
    return (mis[y_name, x_name])
  #
  # No pre-computed MI => Call to computeMutualInfo
  #
  res <- computeMutualInfo (
    x=x_vals, y=y_vals, df_conditioning=cond_data, maxbins=50)

  mi_xy <- ifelse (corrected, res$infok, res$info)
  return (mi_xy)
  }

#-------------------------------------------------------------------------------
# sfp_compute_dpis
#-------------------------------------------------------------------------------
# Compute the DPI on a set of variables and a couple of variables of interest
# Params:
# - couple: a couple of two vois
# - mis: the MI computed for the set of variables and each of the vois
# Return:
# - a data frame, ranked by DPI desc with all the variables passing the DPI test
#-------------------------------------------------------------------------------
sfp_compute_dpis <- function (couple, mis)
  {
  x_name <- couple$voi1
  y_name <- couple$voi2
  mi_xy  <- couple$mi
  #
  # If no mi between the vois => return empty scores (should not happen)
  #
  if ( is.na(mi_xy) || mi_xy == 0 )
    {
    return (data.frame (x     = character(),
                        y     = character(),
                        z     = character(),
                        ais   = character(),
                        mi_xy = numeric(),
                        mi_xz = numeric(),
                        mi_zy = numeric(),
                        dpi   = numeric(),
                        i3    = numeric(),
                        score = numeric(),
                        stringsAsFactors = FALSE) )
    }
  #
  # Extract columns from mis
  #
  mis_x <- mis[, x_name]
  mis_y <- mis[, y_name]
  #
  # Filter rows having mi_xz >= mi_xy and mi_yz >= mi_xy
  #
  keep <- ( !is.na(mis_x) & mis_x >= mi_xy
          & !is.na(mis_y) & mis_y >= mi_xy )
  var_names <- rownames (mis)
  keep[var_names == x_name | var_names == y_name] <- FALSE
  #
  # If no var passes the DPI test => return empty scores
  #
  if ( ! any(keep) )
    {
    return (data.frame (x     = character(),
                        y     = character(),
                        z     = character(),
                        ais   = character(),
                        mi_xy = numeric(),
                        mi_xz = numeric(),
                        mi_zy = numeric(),
                        dpi   = numeric(),
                        i3    = numeric(),
                        score = numeric(),
                        stringsAsFactors = FALSE) )
    }
  #
  # Prepare vectors for the returned data frame
  #
  mi_xz <- mis_x[keep]
  mi_zy <- mis_y[keep]
  z     <- var_names[keep]
  #
  # DPI
  #
  xz  <- mi_xz - mi_xy
  zy  <- mi_zy - mi_xy
  # TODO: check why log1p ... term ? To be comparable with NI3 ?
  dpi <- pmin (xz, zy) - log1p (exp ( -abs(xz - zy) ) )
  # TODO: For now, we keep all results
  # => dpi can be < 0, due to log1p(exp(-abs(xz - zy)))
  # Do we filter out in the future ?
  #
  # Build data frame ordered by dpi desc and return
  #
  ord <- order (dpi, decreasing = TRUE)
  ret <- data.frame (x     = x_name,
                     y     = y_name,
                     z     = z[ord],
                     ais   = couple$ais,
                     mi_xy = mi_xy,
                     mi_xz = mi_xz[ord],
                     mi_zy = mi_zy[ord],
                     dpi   = dpi[ord],
                     i3    = NA_real_,
                     score = NA_real_,
                     stringsAsFactors = FALSE)
  rownames(ret) <- NULL
  return (ret)
  }

#-------------------------------------------------------------------------------
# sfp_update_one_score_with_i3
#-------------------------------------------------------------------------------
sfp_update_one_score_with_i3 <- function (score, x_vals, y_vals,
  input_data, cond_data, corrected, precomputed_scores, has_na)
  {
  z_name <- score$z
  i3 <- NA_real_
  #
  # If previous scores are available, extract the i3 for the z
  #
  if ( ! is.null (precomputed_scores) )
    {
    idx_z <- which (precomputed_scores$z == z_name)
    if (length (idx_z) == 1)
      i3 <- precomputed_scores$i3 [[idx_z]]
    #
    #
    # NB can be NA because of data or change in parameter:
    # - e.g. NA x-z + NA z-y + NA df_cond => no sample => i3 is NA
    # e.g. n_selected = 25 before, new call 50 => only 25 i3 in precomputed
    # As we can know if the NA is a true NA or not, we recompute all NA i3s
    #
    if ( ! is.na (i3) )
      {
      score$i3 <- i3
      score$score <- min (score$dpi, i3, na.rm=T)
      return (score)
      }
    }
  #
  # If no pre-computed i3, compute it
  #
  if ( is.na (i3) )
    {
    z_vals <- input_data[ , z_name, drop=T]
    #
    # TODO not necessary as computeThreePointInfo Check for complete samples ?
    # Ensure we get at least one sample as computeThreePointInfo crashes
    # when called on empty data
    #
    has_na_fct <- has_na | is.na (z_vals)
    if ( ! all (has_na_fct) )
      {
      ret_ni3 <- computeThreePointInfo (
        x=x_vals[!has_na_fct], y=y_vals[!has_na_fct], z=z_vals[!has_na_fct],
        df_conditioning=cond_data[!has_na_fct, , drop=F], maxbins=50)
      i3 <- ifelse (corrected, ret_ni3$i3k, ret_ni3$i3)
      }
    }

  score$i3 <- i3
  # NB: i3 NA could happen e.g. NA x-z + NA z-y + NA df_cond => no sample
  score$score <- ifelse ( is.na (i3), NA_real_, min (score$dpi, i3, na.rm=T) )
  return (score)
  }

#-------------------------------------------------------------------------------
# sfp_compute_top_scores
#-------------------------------------------------------------------------------
# Refine the scoring of top DPIs with I3 scores (like miic for contributors)
# until we get n_selected variables with the highest scores possible
# Params:
# - scores: the scores (only with dpis for now)
# - n_selected, an integer, the number of features to be selected
# - x_vals: the first variable values
# - y_vals: the second variable values
# - input_data: input_data (filtered on a subset of variables)
# - cond_data: data of conditioning set applied for alternate path exploration
# - corrected: a boolean indicating if we used corrected/non corrected I3
# - precomputed_scores: TODO
# Return:
# - the score data frame, ranked by top scores desc and dpi desc
#-------------------------------------------------------------------------------
sfp_compute_top_scores <- function (scores, n_selected,
  x_vals, y_vals, input_data, cond_data, corrected, precomputed_scores)
  {
  if (nrow (scores) <= 0)
    return (scores)

  EPSILON <- 1e-6
  #
  # Pre-identify rows with NA as compute3PointsInfo does not support NA
  #
  has_na <- is.na (x_vals) | is.na (y_vals)
  if ( ! is.null (cond_data) )
    {
    has_na_cond <- apply (cond_data, 1, anyNA)
    has_na <- has_na | has_na_cond
    }
  #
  # Look at the top dpis and compute i3 till we have n_selected
  #
  for ( i in 1:nrow (scores) )
    {
    if ( is.na (scores[i, "score"]) )
      scores[i, ] <- sfp_update_one_score_with_i3 (scores[i, ], x_vals, y_vals,
        input_data, cond_data, corrected, precomputed_scores, has_na)
    #
    # As scores can be lower than dpis, the order of candidate Z can change
    # we can not stop after n_selected as the true tops n_selected scores
    # can be further
    #
    if (i < n_selected)
      next
    tmp_scores <- scores$score[1:i]
    tmp_scores <- tmp_scores[ (!is.na(tmp_scores)) & (tmp_scores > 0) ]
    if (length (tmp_scores) < n_selected)
      next
    #
    # If we are able to find n_selected z about to be selected that explain
    # all the mi between vois, we won't find better z, we can stop
    #
    zs_explain_all <- abs (tmp_scores - scores$mi_xy[[1]])
    are_ok <- (zs_explain_all < EPSILON)
    if (sum (are_ok) >= n_selected)
      break
    #
    # We have n_selected possible z with i3 > 0, but we need to make sure
    # that these n_selected possible z are the best ones,
    # we can stop if the min of scores > next dpi: no more better score to come
    # e.g. n_selected = 2, scores= 10, 7, next dpi 5, next score <= 5, stop
    #      n_selected = 2, scores= 10, 7, next dpi 8, next score <= 8, continue
    #
    # At first, take the n_selected top scores only to find the correct min
    # Avoid case like: (n_selected = 1)
    #  dpi i3 score
    #  10  6  6      => for i = 1, score 6 < next dpi 9 => continue
    #   9  8  8      => for i = 2, top scores [6, 8] => min 6 < next dpi 7 => would continue
    #   7  x  x         but with n_selected = 1, correct test is
    #                   top n_selected scores [8] > next dpi 7 => stop
    #
    tmp_scores <- sort (tmp_scores, decreasing=T)
    tmp_scores <- tmp_scores[1:n_selected]
    min_score <- min (tmp_scores)
    if (  (nrow (scores) <= i)              # no more row for next dpi
       || (scores$dpi[i+1] < min_score) )   # no better score possible
      break
    }
  #
  # Now, we can reorder the top i lines on the score desc
  #
  scores[1:i, ] <- scores[ order (scores$score[1:i], decreasing=T), ]
  return (scores)
  }

#-------------------------------------------------------------------------------
# sfp_update_couple_for_alt
#-------------------------------------------------------------------------------
sfp_update_couple_for_alt <- function (couple, scores)
  {
  # If no z passes the dpi test, no candidate top_ai to condition
  # and evaluate the need to explore an alternate path.
  # As top_ai_for_alt and mis_remain_for_alt are initialized with NA
  # => nothing to do => return the couple unchanged
  #
  if (nrow (scores) <= 0)
    return (couple)

  EPSILON <- 1e-6
  #
  # The scores are already ordered down by top contributors => pick the first
  #
  row_top_ai <- 1
  #
  # If i3 is na, weird case and if i3 < 0 => common child, can not be use to
  # As top_ai_for_alt and mis_remain_for_alt are initialized with NA
  # => nothing to do => return the couple unchanged
  #
  if ( is.na (scores[row_top_ai, "i3"]) || (scores[row_top_ai, "i3"] < 0) )
    return (couple)
  #
  # Check remaining MI when subtracting i3
  #
  if (abs (scores[row_top_ai, "mi_xy"] - scores[row_top_ai, "i3"]) <= EPSILON)
    # Top ai found, MI fully explained => no alternate path exploration
    couple[ c("top_ai_for_alt", "mis_remain_for_alt") ] <- list (
      scores[row_top_ai, "z"], 0)
  else # Top ai found, MI not fully explained => alternate path exploration
    couple[ c("top_ai_for_alt", "mis_remain_for_alt") ] <- list (
      scores[row_top_ai, "z"], scores[row_top_ai, "mi_xy"] - scores[row_top_ai, "i3"])

  return (couple)
  }

#-------------------------------------------------------------------------------
# sfp_eval_one_couple
#-------------------------------------------------------------------------------
sfp_eval_one_couple <- function(
  depth, couple, input_data, extra_data, vars_remain,
  precomputed_mis, precomputed_scores, n_selected,
  corrected, skip_cheks, n_threads, verbose, verbose_str)
  {
  voi1 <- couple$voi1
  voi2 <- couple$voi2
  mi_vois <- couple$mi
  ais <- couple$ais
  ais_vect <- c()
  if ( ! is.na (ais) )
    ais_vect <- strsplit (ais, ",")[[1]]
  #
  # Prepare data for MI estimation
  #
  voi_names <- c()
  extra_names <- c()
  if ( voi1 %in% colnames(extra_data) )
    extra_names <- voi1
  else
    voi_names <- voi1
  if ( voi2 %in% colnames(extra_data) )
    extra_names <- c(extra_names, voi2)
  else
    voi_names<- c(voi_names, voi2)
  vars_to_keep <- c (voi_names, vars_remain)
  if ( length(vars_to_keep) != length (unique (vars_to_keep)) )
    stop ("TODO length(vars_to_keep) != unique (length(vars_to_keep) )")
  data_filt <- sfp_filter_data (voi1=voi1, voi2=voi2,
    input_data=input_data, input_names=vars_to_keep,
    extra_data=extra_data, extra_names=extra_names, ais_vect=ais_vect)
  #
  # Compute MIs with the 2 VOIs
  #
  if (verbose >= 3)
    cat_for_rewrite (paste0 (verbose_str, ", computing MIs ...") )
  mis <- compute_mi_batch (input_data=data_filt[["input_data"]],
    var_of_interest_names=voi_names, var_of_interest_values=data_filt[["extra_data"]],
    df_conditioning=data_filt[["cond_data"]], unit="log_conf", corrected=corrected,
    precomputed_mis=precomputed_mis, skip_cheks=skip_cheks,
    n_threads=n_threads, verbose=ifelse (verbose >= 2, 2, 0),
    verbose_start=verbose_str, verbose_end="")
  #
  # The mis can contain more rows/columns than we are interested now
  #
  mis_filt <- mis [rownames(mis) %in% vars_to_keep,
                   colnames(mis) %in% c(voi1, voi2),
                   drop=F]
  #
  # Pick (or compute) reference MI for voi1-voi2
  #
  couple$mi <- sfp_get_mi (mis=mis_filt,
    x_name=voi1, x_vals=data_filt[["voi1_vals"]],
    y_name=voi2, y_vals=data_filt[["voi2_vals"]],
    cond_data=data_filt[["cond_data"]], corrected=corrected)
  #
  # Ends if MI estimation returns an error or no MI between vois
  #
  if ( is.na (couple$mi) || (couple$mi <= 0) )
    {
    scores_empty <- data.frame (x     = character(),
                                y     = character(),
                                z     = character(),
                                ais   = character(),
                                mi_xy = numeric(),
                                mi_xz = numeric(),
                                mi_zy = numeric(),
                                dpi   = numeric(),
                                i3    = numeric(),
                                score = numeric(),
                                stringsAsFactors = FALSE)
    couples_deeper_empty <- list ("couples" = data.frame ("voi1" = character(),
                                    "voi2" = character(), stringsAsFactors=F),
                                  "vars_remain" = list() )
    ret <- list ("couple" = couple,
                 "mis" = mis,
                 "scores" = scores_empty,
                 "features" = c(),
                 "couples_deeper" = couples_deeper_empty )
    return (ret)
    }
  #
  # Compute DPTs
  #
  if (verbose >= 3)
    cat_for_rewrite (paste0 (verbose_str, ", computing DPIs ...") )
  scores <- sfp_compute_dpis (couple=couple, mis=mis_filt)
  #
  # Compute I3 and scores till we find the top n_selected vars
  #
  if (verbose >= 3)
    cat_for_rewrite (paste0 (verbose_str, ", computing top scores ...") )
  scores <- sfp_compute_top_scores (
    scores=scores, n_selected=n_selected,
    x_vals=data_filt[["voi1_vals"]], y_vals=data_filt[["voi2_vals"]],
    input_data=data_filt[["input_data"]], cond_data=data_filt[["cond_data"]],
    corrected=corrected, precomputed_scores=precomputed_scores)
  #
  # Evaluate if alternate path exploration would be needed
  #
  if (verbose >= 3)
    cat_for_rewrite (paste0 (verbose_str, ", test if alternate needed ...") )
  couple <- sfp_update_couple_for_alt (couple=couple, scores=scores)
  #
  # Using scores, take top n_selected vars
  #
  if (verbose >= 3)
    cat_for_rewrite (paste0 (verbose_str, ", selecting tops ...") )
  new_features <- c()
  couples_deeper <- list ("couples" = data.frame ("voi1" = character(),
                            "voi2" = character(), stringsAsFactors=F),
                          "vars_remain" = list() )
  #
  # Scores are ordered desc by score then dpi but no guarantee that the
  # n_selected tops are ok: especially if we have few variables remaining,
  # dpi can be a bit < 0 with the - log1p (exp (..?) ) term
  # and scores could be < 0 for common child
  #
  tmp_scores <- scores[scores$dpi > 0 & (is.na (scores$score) | scores$score > 0),
                       , drop=F]
  #
  # If nrow (tmp_scores == 0), no var passing "corrected" dpi and/or i3 test
  # => No feature to select
  # As couple$features is initialized with NA, nothing to do => return
  # NB: we return scores and not tmp_scores to see if regularization affect dpi
  #
  if ( nrow (tmp_scores) <= 0 )
    {
    couples_deeper_empty <- list ("couples" = data.frame ("voi1" = character(),
                                    "voi2" = character(), stringsAsFactors=F),
                                  "vars_remain" = list() )
    ret <- list ("couple" = couple,
                 "mis" = mis,
                 "scores" = scores,
                 "features" = c(),
                 "couples_deeper" = couples_deeper_empty )
    return (ret)
    }
  #
  # We can select new features => update the couple and prepare couples deeper
  #
  new_features <- tmp_scores$z [1:min (n_selected, nrow(tmp_scores) )]
  couple$features <- paste (new_features, collapse=",")
  if (verbose >= 3)
    cat_for_rewrite (paste0 (verbose_str, ", preparing couples for next depth ...") )
  couples_deeper <- sfp_prepare_couples_deeper (couple, new_features, mis_filt)
  #
  # Prepare structure to be returned
  #
  return ( list ("couple" = couple,
                 "mis" = mis,
                 "scores" = scores,
                 "features" = new_features,
                 "couples_deeper" = couples_deeper) )
  }

#-------------------------------------------------------------------------------
# sfp_plot
#-------------------------------------------------------------------------------
# Params (internal):
# - vois: a list, the vois returned by sfp_check_vois
# - features: the list of selected features
# - mis: the MI matrix computed with the user vois and without condtionning
# Params (can be provided by the user via the ... in selectFeaturesPath):
# - n_bins: an integer >= 1, default 20, the maximum number of bins plotted,
#   used to place the features along the x axis
# - max_feats_per_bin: an integer >= 1, default 25, the maximum number of
#   features displayed per bin.
# - annotate: a boolean, default T, indicates the total features found
#   and, if needed, the number not displayed on the plot
# - x_lab: the x label, default "Features for vois_side1-vois_side2"
# - font_size: an integer >= 1, default 11, the font size
#-------------------------------------------------------------------------------
sfp_plot <- function (vois, features, mis, n_bins=20, max_feats_per_bin=25,
                      annotate=T, x_lab=NULL, font_size=11)
  {
  # Check params
  #
  n_bins <- check_param_int (n_bins,
    "maximum plot number of bins", default=10, min=1)
  max_feats_per_bin <- check_param_int (max_feats_per_bin,
    "maximum features plotted per bin", default=25, min=1)
  annotate <- check_param_logical (annotate, "plotting of annotation", default=T)
  if ( is.null (x_lab) )
    x_lab <- paste0 ("Features for ",
                     paste (vois$side1$all_voi_names, collapse=", "), " - ",
                     paste (vois$side2$all_voi_names, collapse=", ") )
  else
    x_lab <- as.character (x_lab)
  font_size <- check_param_int (font_size, "font size", default=11, min=1)

  if (length(features) > 0)
    {
    # Extract mis of features with the vois
    #
    mis_main_path <- mis[rownames(mis) %in% features,
      colnames(mis) %in% vois$all$all_voi_names, drop=F]
    mis_main_side1 <- mis_main_path[,
      colnames(mis_main_path) %in% vois$side1$all_voi_names, drop=F]
    mis_main_side2 <- mis_main_path[,
      colnames(mis_main_path) %in% vois$side2$all_voi_names, drop=F]
    #
    # Compute base position for each feature
    #
    feat_pos <- unlist (lapply (features, FUN=function (one_feat) {
      mean_side1 <- mean (mis_main_side1[ one_feat, ], na.rm = T)
      mean_side2 <- mean (mis_main_side2[ one_feat, ], na.rm = T)
      if ( is.na (mean_side1) || (mean_side1 == 0) )
        {
        if ( is.na (mean_side2) || (mean_side2 == 0) )
          return (0.5) # Unknown => to the middle
        return (0.99)  # On side2 => 0.99
        }
      if ( is.na (mean_side2) || (mean_side2 == 0) )
        return (0.01)  # On side1 => 0.01
      pos <- round (mean_side2 / (mean_side1 + mean_side2), 3)
      pos <- max (0.01, min (0.99, pos) )
      return (pos)
      } ) )
    #
    # Use the number of bins to position the feaures along the x axis
    #
    bin_pos <- seq (from=0, to=1, length.out=n_bins+2)

    bin_ranges <- bin_pos[1:(length(bin_pos)-1)] + bin_pos[[2]] / 2
    bin_ranges[ c (1, length(bin_ranges) ) ] <- c (0,1)

    feat_ind <- unlist (lapply (feat_pos, FUN=function(pos) {
      which (bin_ranges > pos)[[1]] }) )

    feat_counts <- table (feat_ind)
    #
    # The plot locations are known for each feature
    # Look if too much to display at each location
    #
    df_plot <- data.frame ("feature"=features,
                           "index"=feat_ind,
                           "position"=bin_pos[feat_ind],
                           stringsAsFactors=F)
    df_plot <- df_plot[order(df_plot$feature), , drop=F]
    one_grp <- names(feat_counts)[[1]]
    for (one_grp in names(feat_counts))
      {
      rows_grp <- which(df_plot$index == one_grp)
      if (length(rows_grp) > max_feats_per_bin)
        {
        df_plot[rows_grp[max_feats_per_bin], "feature"] <- "..."
        df_plot[rows_grp[(max_feats_per_bin+1):length(rows_grp)],
                "feature"] <- NA_character_
        }
      }
    df_plot <- df_plot[ ! is.na(df_plot$feature), , drop=F]
    }
  #
  # Base plot with the 2 set of vois
  #
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
  # Plot the features at their respective position
  #
  if (length(features) > 0)
    {
    plot_positions <- c()
    plot_texts <- c()
    one_pos <- 0.5
    for (one_pos in unique (df_plot$position) )
      {
      feats_to_plot <- df_plot$feature[df_plot$position == one_pos]
      plot_positions <- c (plot_positions, one_pos)
      plot_texts <- c (plot_texts, paste (feats_to_plot, collapse="\n") )
      }
    g <- g +
      ggplot2::geom_text (
        ggplot2::aes (x=plot_positions, y=y_feat, label=plot_texts),
        size=font_size*0.8/ggplot2::.pt, hjust=0.5, vjust=0.5, color="black") +
      ggplot2::geom_vline ( xintercept=plot_positions, linewidth=0.6,
                            linetype="dotted", color="darkgrey")
    xticks_text <- c (xticks_text, as.character(round (plot_positions, 4) ) )
    xticks_pos <- c (xticks_pos, plot_positions)
    }
  #
  # Set ticks
  #
  g <- g + ggplot2::scale_x_continuous (breaks=xticks_pos, label=xticks_text)
  #
  # Annotation
  #
  if (annotate)
    {
    if (length(features) <= 0)
      label_txt <- paste0 ("No feature found")
    else
      {
      shown_feats <- df_plot$feature[ df_plot$feature != "..." ]
      n_feat_not_shown <- length (features) - length (shown_feats)
      if (n_feat_not_shown <= 0)
        label_txt <- paste0 (length(features), " features")
      else
        label_txt <- paste0 (length(features), " features",
                             ",\n(", n_feat_not_shown, " not shown)")
      }
    g <- g + ggplot2::geom_text (ggplot2::aes (x=0.875, y=1, label=label_txt),
      size=font_size*0.8/ggplot2::.pt, hjust=0.5, vjust=1)
    }
  return (g)
  }

#===============================================================================
# FUNCTIONS (exported)
#===============================================================================
# selectFeaturesPath
#-------------------------------------------------------------------------------
#' Features selection on the path between variables of interest
#'
#' @description Select the variables that are the most likely to be in the path
#' between two set of variable(s) of interest.
#' Variables selection is performed using a score combining Data Processing
#' inequality (DPI) and 3 points information (see Affeldt 2015).\cr
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
#' @param n_selected [a positive integer, optional, 1 by default]
#'
#' The number of features selected at each step, for each path explored
#' and at each iteration of the recursion. Increasing this value
#' will identify more variables as possible features at each step
#' with the consequence to reach quicker a target \emph{max_features}
#' while increasing the false positives.
#'
#' @param n_alternate [a positive integer, optional, 2 by default]
#'
#' The number of alternate paths to explore at each step of the recursion.
#' Decreasing this value speeds up the process while reducing the number
#' of possible paths discovered. Increasing it has opposite effect,
#' more features from alternate possible paths at the cost of an increased
#' processing time.
#'
#' @param max_depth [a positive integer, optional, 4 by default]
#'
#' The maximum depth of the recursion.
#' Decreasing this value speeds up the process while reducing the granularity
#' of the path discovery. Increasing it has opposite effect, more features
#' can be collected along the path at the cost of an increased processing time.
#'
#' @param max_features [a positive integer, optional, 250 by default]
#'
#' The maximum total number of features to select. Set by default to a large
#' number of features as the goal is to not reach this maximal number while
#' finding all variables on the path between the variables of interest.
#' Consider to tune the \emph{max_depth} and \emph{n_alternate} parameters
#' to identify more or less features.
#'
#' @param precomputed_values [a list, optional, NULL by default]
#'
#' If path features selection have been previously performed on the same
#' conditions (same variables of interest, same \emph{input_data},
#' same \emph{corrected} parameter),
#' Supplying the result of a previous run speeds up the process
#' as the existing computed values will be reused.
#' The expected format of the list is the same as the function returned value,
#' with at least an item "mis" containing a sub-list of \emph{mis} matrices
#' and an "scores" item containing a sub-list of \emph{scores} data frames.
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
#' @param verbose [an integer, optional, 2 by default]
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
#' Please note on the plot rendering that positions are indicative
#' as they are only based on the ratio of MI shared between each feature
#' and the variables of interest of each side.
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
#' \item \emph{mis}: a list of matrices with the MIs between the variables in
#'   \emph{input_data} (as rows) and the variable(s) of interest (as columns).
#'   Row and column names are sorted alphabetically.\cr
#'   The MIs can be corrected or not, depending on lthe \emph{corrected}
#'   parameters. If pre-computed values were supplied, new values computed
#'   are added to the existing matrices.
#' \item \emph{couples}: a data frames with, at each depth, the couples
#'   of variables used with their MIs and features selected.
#'   Please note that couples already encountered are not included
#'   as no recursion is performed on duplicates.
#' \item \emph{scores}: a list of data frames with, for all pairs of variables
#'   passing the dpi test, the information about the contionning set used,
#'   the MI values, the DPI test, the 3 points information and the scores.
#' \item \emph{plot}: when the \emph{plot} parameter is turned to TRUE,
#'   a plot with the features selected, NULL otherwise.
#'   Please note on the plot rendering that positions are indicative
#'   as they are only based on the ratio of MI shared between each feature
#'   and the variables of interest of each side.
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
#' ret <- selectFeaturesPath (input_data=df_data,
#'                            var_of_interest_values_side1=df_external_meta,
#'                            var_of_interest_names_side2="TP53")
#' message ("Features selected: ", paste (ret$features, collapse=", ") )
#'
#' # Same features selection with reuse of the previous run above and plot
#' ret <- selectFeaturesPath (input_data=df_data,
#'                            var_of_interest_values_side1=df_external_meta,
#'                            var_of_interest_names_side2="TP53",
#'                            precomputed_values=ret, plot=TRUE)
#' print (ret$plot)
#'
#' # Features selection using multiple variables on each side
#' # (limiting the recursion depth to speed up the example)
#' # and selecting two features per step
#' ret <- selectFeaturesPath (input_data=cosmicCancer,
#'                            var_of_interest_names_side1=c("TP53", "Ploidy"),
#'                            var_of_interest_names_side2=c("FOXM1", "AURKA"),
#'                            max_depth=2, n_selected=2)
#' message ("Features selected: ", paste (ret$features, collapse=", ") )
#'
#' # Similar features selection with plotting using a customized rendering
#' # (reusing the previous run above and limiting the recursion depth
#' # to speed up the example)
#' ret <- selectFeaturesPath (cosmicCancer,
#'                            var_of_interest_names_side1=c("TP53", "Ploidy"),
#'                            var_of_interest_names_side2=c("FOXM1", "AURKA"),
#'                            max_depth=2, precomputed_values=ret,
#'                            plot=TRUE, font_size=12,
#'                            x_lab="My features selection on CosmicCancer")
#' print (ret$plot)
#' }
#'
#' @export
#-------------------------------------------------------------------------------
selectFeaturesPath <- function (input_data,
  var_of_interest_names_side1=NULL, var_of_interest_values_side1=NULL,
  var_of_interest_names_side2=NULL, var_of_interest_values_side2=NULL,
  n_selected=1, n_alternate=2, max_features=250, max_depth=4,
  precomputed_values=NULL, corrected=T, skip_cheks=F, n_threads=1, verbose=2,
  plot=F, ...)
  {
  # Check parameters
  #
  input_data <- sf_check_input_data (input_data=input_data)
  vois <- sfp_check_vois (input_data=input_data,
    var_of_interest_names_side1=var_of_interest_names_side1,
    var_of_interest_values_side1=var_of_interest_values_side1,
    var_of_interest_names_side2=var_of_interest_names_side2,
    var_of_interest_values_side2=var_of_interest_values_side2)
  n_selected <- check_param_int (n_selected,
    "number of selected features at each step", default=1, min=1)
  n_alternate <- check_param_int (n_alternate,
    "number of alternate paths", default=2, min=1)
  max_depth <- check_param_int (
    max_depth, "maximum depth", default=4, min=1)
  max_features <- check_param_int (max_features,
    "maximum number of features", default=250, min=1)
  precomputed_values <- sfp_check_precomputed_values (precomputed_values=precomputed_values)
  corrected <- check_param_logical (corrected, "corrected", T)
  skip_cheks <- check_param_logical (skip_cheks, "skip checks", F)
  n_threads <- check_param_int (
    n_threads, "number of threads", default=1, min=1)
  verbose <- check_param_int (verbose, "verbose", 2, 0, 4)
  plot <- check_param_logical (plot, "plot", F)
  # TODO ADD CHECK ON ... params names
  #
  # Init of structures to be returned
  #
  feat_kept <- c()
  all_couples <- NULL
  #
  # Initialize the couples of vois (one from each side) to be evaluated
  # The structures are reused as we go deeper
  #
  list_deeper <- sfp_init_couples_deeper (vois1=vois$side1$all_voi_names,
    vois2=vois$side2$all_voi_names, vars_to_test=colnames(input_data) )
  #
  # Pre-computed values
  #
  all_mis <- list()
  all_scores <- list()
  if ( "mis" %in% names (precomputed_values) )
    all_mis <- precomputed_values[["mis"]]
  if ( "scores" %in% names (precomputed_values) )
    all_scores <- precomputed_values[["scores"]]
  #
  # Recursive loop
  #
  too_much_feats <- F
  for (depth in 1:max_depth)
    {
    # Init variables for the new depth
    #
    couples <- list_deeper$couples
    couples_vars_remain <- list_deeper$vars_remain
    if (nrow (couples) <= 0)
      break
    #
    # Reset variables for the next depth
    #
    list_deeper <- list ("couples" = data.frame ("voi1" = character(),
                           "voi2" = character(), stringsAsFactors=F),
                         "vars_remain" = list() )
    #
    # Update couples for this depth, ready for main path exploration
    #
    couples <- transform (couples,
                          depth = depth,
                          alt = NA_integer_,   # Just here for the column order:
                          mi = NA_real_,       # depth -> alt -> mi -> ais
                          ais = NA_character_) # ! Init ais to NA for main path !
    #
    # Search for main path (idx=0), 1st alternate, 2nd alternate, ...
    #
    for (alt_idx in 0:n_alternate)
      {
      # Update/reset couples for this alternate with values to "not evaluated"
      #
      couples <- transform (couples,
                            alt = alt_idx,
                            mi = NA_real_,
                            features = NA_character_,
                            top_ai_for_alt = NA_character_,
                            mis_remain_for_alt = NA_real_)
      #
      # Init the list of couples evaluation (oce = One Couple Evaluation)
      #
      list_oce <- vector ( "list", nrow(couples) )
      couples_done <- rep ( F, nrow(couples) )
      #
      # feat_kept are the features identified and selected (< max_features)
      # feat_new will be used to test if > max_features, then return feat_kept
      #
      feat_new <- feat_kept
      for ( cpl_idx in 1:nrow(couples) )
        {
        verbose_str <- ""
        if (verbose >= 2)
          {
          ai_str <- ""
          if ( ! is.na (couples[cpl_idx, "ais"]) )
            ai_str <- paste0 ("|", couples[cpl_idx, "ais"])
          verbose_str <- paste0 ("Depth ", depth,
            ", ", couples[cpl_idx, "voi1"], "-", couples[cpl_idx, "voi2"], ai_str,
            " (", cpl_idx, "/", nrow(couples), ")")
          if (verbose == 3)
            verbose_str <- paste0 (verbose_str,
              ", ", length(feat_kept),
              " + ", length(feat_new) - length(feat_kept), " features")
          else if (verbose >= 4)
            verbose_str <- paste0 (verbose_str,
              ", ", length(couples_vars_remain[[cpl_idx]]), " vars remaining",
              ", ", length(feat_kept),
              " + ", length(feat_new) - length(feat_kept), " features")
          cat_for_rewrite (verbose_str)
          }
        #
        # If the couple with the same ais has been done, no need to redo it
        # (as voi1 is always < voi2 alphabetically, no need to test voi2-voi1)
        #
        if (any (  all_couples$voi1 == couples[cpl_idx, "voi1"]
                 & all_couples$voi2 == couples[cpl_idx, "voi2"]
                 & ( (  is.na (all_couples$ais) &  is.na (couples[cpl_idx, "ais"]) )
                   | ( !is.na (all_couples$ais) & !is.na (couples[cpl_idx, "ais"])
                     & all_couples$ais == couples[cpl_idx, "ais"] ) ) ) )
          {
          if (!is.na (couples[cpl_idx, "ais"]))
            {
            print ("")
            print ("===================================================================")
            print ("Couple with ai already evaluated:")
            print (couples[cpl_idx, ])
            aze = which ( ( (all_couples$voi1 == couples[cpl_idx, "voi1"])
                  & (all_couples$voi2 == couples[cpl_idx, "voi2"]) )
                  & ( (   is.na (all_couples$ais)  &   is.na (couples[cpl_idx, "ais"])   )
                    | ( (!is.na (all_couples$ais)) & (!is.na (couples[cpl_idx, "ais"]) )
                      & (all_couples$ais == couples[cpl_idx, "ais"]) ) ) )[[1]]
            print ("Couple found:")
            print (all_couples[aze, , drop=F])
            print ("===================================================================")
            }

          couples_done[[cpl_idx]] <- T
          list_oce[[cpl_idx]] <- NA
          next
          }
        #
        # Extract mis and scores from previous run if supplied
        #
        if (is.na (couples[cpl_idx, "ais"]) )
          one_key_mis <- "|[none]"
        else
          one_key_mis <- paste0 ("|'", couples[cpl_idx, "ais"], "'")
        one_key_scores <- paste0 ("'" , couples[cpl_idx, "voi1"],
                                 "'-'", couples[cpl_idx, "voi2"],
                                  "'" , one_key_mis)
        mis_loop <- NULL
        if ( one_key_mis %in% names (all_mis) )
          mis_loop <- all_mis[[one_key_mis]]
        scores_loop <- NULL
        if ( one_key_scores %in% names (all_scores) )
          scores_loop <- all_scores[[one_key_scores]]
        #
        # Evaluate the couple
        #
        eoc_ret <- sfp_eval_one_couple (depth=depth,
          couple=couples[cpl_idx, , drop=T], input_data=input_data,
          extra_data=vois$all$var_of_interest_values,
          vars_remain=couples_vars_remain[[cpl_idx]],
          precomputed_mis=mis_loop, precomputed_scores=scores_loop,
          n_selected=n_selected, corrected=corrected,
          skip_cheks=skip_cheks, n_threads=n_threads,
          verbose=verbose, verbose_str=verbose_str)
        #
        # Update pre-computed structures and couple
        #
        all_mis[[one_key_mis]] <- eoc_ret$mis
        all_scores[[one_key_scores]] <- eoc_ret$scores
        couples[cpl_idx, ] <- eoc_ret$couple
        #
        # Report on the couple
        #
        feat_new <- unique ( c (feat_new, eoc_ret$features) )
        if (verbose >= 2)
          {
          verbose_str <- paste0 ("Depth ", depth,
            ", ", couples[cpl_idx, "voi1"], "-", couples[cpl_idx, "voi2"], ai_str,
            " (", cpl_idx, "/", nrow(couples), ")")
          if (verbose == 3)
            verbose_str <- paste0 (verbose_str,
              ", ", length(feat_kept),
              " + ", length(feat_new) - length(feat_kept), " features")
          else if (verbose >= 4)
            verbose_str <- paste0 (verbose_str,
              ", ", length(couples_vars_remain[[cpl_idx]]), " vars remaining",
              ", ", length(feat_kept),
              " + ", length(feat_new) - length(feat_kept), " features")
          verbose_str <- paste0 (verbose_str, " => Done\n")
          cat_for_rewrite (verbose_str)
          }
        #
        # Test if the max_features param is exceeded to stop the process
        #
        if (length (feat_new) > max_features)
          {
          too_much_feats <- T
          break
          }
        #
        # Ok to continue, memorize the eoc, they will be processed after all
        # couples are done (once we know that, by adding all the features found
        # for this depth and alternate path, we do not exceed the max_features)
        #
        list_oce[[cpl_idx]] <- eoc_ret
        }
      #
      # All couples for this depth + main or alternate path done
      #
      # Filter out couples that were done previously as we simply skip them
      # when iterating over the couples
      #
      couples_vars_remain <- couples_vars_remain [!couples_done]
      couples <- couples[!couples_done, , drop=F]
      list_oce <- list_oce[!couples_done]
      #
      # If no couple remains => we skipped all couples for this depth / path
      # => Nothing to do more (no new feature found, no new couple evaluated)
      # => Go to the next depth
      #
      if (nrow (couples) <= 0)
        break
      #
      # This round has evaluated some couples => Check if we exceeded the
      # max features => if so, we don't add the last features found
      #
      if (too_much_feats)
        break
      #
      # Some couples evaluated + max features not exceeded
      # => Update the list of features selected and list of all couples
      #
      if ( is.null(all_couples) )
        all_couples <- couples
      else
        all_couples <- rbind (all_couples, couples)
      feat_kept <- feat_new
      #
      # Memorize the structures to go deeper before looking for alternate path
      #
      for ( cpl_idx in 1:nrow(couples) )
        {
        # For next depth
        #
        if (nrow (list_oce[[cpl_idx]]$couples_deeper$couples) > 0)
          {
          list_start <- nrow (list_deeper$couples) + 1
          list_end <- (list_start - 1) + nrow (list_oce[[cpl_idx]]$couples_deeper$couples)
          list_deeper$couples <- rbind (list_deeper$couples,
                                        list_oce[[cpl_idx]]$couples_deeper$couples)
          list_start <- length (list_deeper$vars_remain)
          for (couples_deeper_idx in 1:nrow(list_oce[[cpl_idx]]$couples_deeper$couples))
            list_deeper$vars_remain[[list_start+couples_deeper_idx]] <-
              list_oce[[cpl_idx]]$couples_deeper$vars_remain[[couples_deeper_idx]]
          }
        }
      #
      # Prepare the structures for next alternate path
      # Filter couples on those still having MI not explained by the top ai
      #
      couples_vars_remain <- couples_vars_remain[
          ( ! is.na (couples$mis_remain_for_alt) )
        & (couples$mis_remain_for_alt > 0) ]
      couples <- couples[ ( ! is.na (couples$mis_remain_for_alt) )
                        & (couples$mis_remain_for_alt > 0), , drop=F]
      #
      # If no more couple with a remaining MI when conditioning in the top ai
      # No alternate path exploration needed for any couple => go to next depth
      #
      if ( nrow (couples) <= 0)
        break
      #
      # Update the ais to use for the alternate path exploration
      #
      for ( i in 1:nrow(couples))
        couples[i, "ais"] <- ifelse (is.na (couples[i, "ais"]),
          couples[i, "top_ai_for_alt"],
          paste0 (couples[i, "ais"], ",", couples[i, "top_ai_for_alt"]) )
      }
    if (too_much_feats)
      break
    #
    # Report when depth fully done
    #
    if ( (verbose >= 1) )
      cat_for_rewrite (paste0 ("Depth ", depth,  " fully explored, ",
        length (feat_kept), " features selected\n") )
    }
  #
  # Final displays
  #
  if (verbose >= 1)
    {
    # NB for alt_idx 0 => means no feature was added since the previous depth
    # => we already printed full exploration of previous depth
    #
    if ( (too_much_feats) )
      {
      if (verbose >= 2)
        {
        if (alt_idx == 0)
          cat_for_rewrite (paste0 ("Depth ", depth, ", main path",
            ", maximum number of features exceeded\n") )
        else
          cat_for_rewrite (paste0 ("Depth ", depth, ", alternate path ", alt_idx,
            ", maximum number of features exceeded\n") )
        }
      }
    cat_for_rewrite (paste0 (length (feat_kept), " features selected\n") )
    }
  #
  # Prepare returned values
  #
  ret <- list ("features"=feat_kept, "couples"=all_couples,
               "mis"=all_mis, "scores"=all_scores)
  #
  # Optional plot
  #
  if (plot)
    {
    if ( base::requireNamespace("ggplot2", quietly=TRUE) )
      ret$plot <- sfp_plot (
        vois=vois, features=feat_kept, mis=all_mis[["|[none]"]], ...)
    else
      miic_warning ("Path features selection", "Plotting requires ggplot2.")
    }
  return (ret)
  }

