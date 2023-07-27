#===============================================================================
# CONSTANTS
#===============================================================================
MIIC_VALID_LATENT <- c ("orientation", "yes", "no")
MIIC_VALID_CONSISTENT <- c ("no", "orientation", "skeleton")

MIIC_CONTINUOUS_TRESHOLD <- 5

STATE_ORDER_VALID_COLUMS <- c ("var_names", "var_type",
    "levels_increasing_order", "is_contextual", "is_consequence")

#===============================================================================
# FUNCTIONS
#===============================================================================
# list_to_str
#-------------------------------------------------------------------------------
# Utility function to transform the first n_max items of a list into a
# displayable (comma + space separated) string: "item1, item2, item3, ..."
# - list: a list
# - n_max: int, optional, NULL by default, maximum of items used. If NULL,
#   all items are used. If the list has more than n_max items, ", ..." is added
# return: string, the displayable string
#-------------------------------------------------------------------------------
list_to_str <- function (list, n_max=NULL) {
  if (is.null (list) )
    return ("NULL")
  if (length (list) == 0)
    return ("")
  if ( (! is.null (n_max)) && (length (list) > n_max) )
    list <- c (list[1:n_max], "...")
  ret <- paste (unlist(list), collapse=", ")
  return (ret)
}

#-------------------------------------------------------------------------------
# miic_warning
#-------------------------------------------------------------------------------
# Utility function to raise a warning
#-------------------------------------------------------------------------------
miic_warning <- function (context, ...) {
  warning (paste0 ("Warning in ", context, ": ", ...), call.=FALSE)
}

#-------------------------------------------------------------------------------
# Check input data
#-------------------------------------------------------------------------------
# input_data: a dataframe with variables as columns and rows as samples
#-------------------------------------------------------------------------------
check_input_data <- function (input_data) {
  if ( is.null(input_data) )
    stop ("The input data is required.", call.=FALSE)
  if ( ! is.data.frame (input_data) )
    stop ("The input data must be a dataframe.", call.=FALSE)
  if (nrow (input_data) == 0)
    stop ("The input data is empty.", call.=FALSE)
  if (ncol (input_data) == 0)
    stop ("The input data has no variable.", call.=FALSE)
  rows_only_na <- rowSums (is.na (input_data)) == ncol (input_data)
  input_data <- input_data [!rows_only_na, ]
  if (nrow (input_data) == 0)
    stop ("The input data contains only NAs", call.=FALSE)
  if ( any (rows_only_na) )
    miic_warning ("input data", "the input data contains ", sum(rows_only_na),
             " row(s) with only NAs. These row(s) will be removed.")

  n_unique_vals <- unlist (lapply (input_data, function (x) {
    length (unique (x[!is.na(x)] ) ) } ) )
  vars_constant = (n_unique_vals <= 1)
  if ( any (vars_constant) ) {
    msg_str <- list_to_str (colnames (input_data)[vars_constant], n_max=10)
    if (sum (vars_constant) == 1)
      miic_warning ("input data", "the variable ", msg_str, " is constant.",
        " Such variable can not be connected and should be removed.")
    else
      miic_warning ("input data", sum(vars_constant), " variables (", msg_str,
        ") are constant. Such variables can not be connected",
        " and should be removed.")
  }
  return (input_data)
}

#-------------------------------------------------------------------------------
# check_state_order
#-------------------------------------------------------------------------------
# Check the state_order: verify that the var_names column is supplied, that all
# variables are present in the var_names. For each extra columns, check the
# values supplied and correct them if necessary. Wrong or not supplied values
# are assigned to default values.
# Stop if all variables are contextual or consequences, other issues will
# raise a warning and be "fixed" by using default values.
# Parameters:
# - input_data: a dataframe with variables as columns and rows as samples
# - state_order: a dataframe, can be NULL.
#   possible/expected columns are:
#   * var_names: the list of all variables in the input data
#   * var_type: 0=discrete, 1=continuous (default deduced from the input data,
#     note that numerical variables with less than x unique values will be
#     discrete)
#   * levels_increasing_order: NA for continuous variables. For discrete,
#     can be NA or the full ordered list of the unique values. (default NA)
#   * is_contextual: 0=not contextual, 1=contextual (default 0)
#   * is_consequence: 0=not consequence, 1=consequence (default 0)
# Return: the checked and eventually generated or completed state order dataframe
#-------------------------------------------------------------------------------
check_state_order <- function (input_data, state_order) {
  data_var_names <- colnames (input_data)
  n_vars <- length (data_var_names)
  #
  # Basic checks
  #
  if ( is.null (state_order) )
    state_order <- data.frame ("var_names"=data_var_names, stringsAsFactors=F)
  if ( ! is.data.frame(state_order) ) {
    miic_warning ("state order",
      "the supplied state_order is not a dataframe and will be ignored.")
    state_order <- data.frame ("var_names"=data_var_names, stringsAsFactors=F)
  }

  if ( ! ("var_names" %in% colnames (state_order)) ) {
    miic_warning ("state order", "the column var_names is missing,",
                  " the supplied state_order will be ignored.")
    state_order <- data.frame ("var_names"=data_var_names)
  }
  #
  # Check if the state_order columns are valid
  #
  valid_cols <- STATE_ORDER_VALID_COLUMS
  mismatch <- is.na (match (colnames (state_order), valid_cols))
  if ( any (mismatch) ) {
    msg_str <- list_to_str (colnames (state_order)[mismatch], n_max=10)
    if (sum (mismatch) == 1)
      miic_warning ("state order", "the column ", msg_str,
        " is not valid and will be ignored.")
    else
      miic_warning ("state order", sum (mismatch), " columns (", msg_str,
        ") are not valid and will be ignored.")
    state_order <- state_order[, !mismatch, drop=FALSE]
  }
  #
  # We ensure that the var_names column is the first
  #
  idx_var_names <- which (colnames(state_order) == "var_names")
  idx_others <- 1:ncol (state_order)
  idx_others <- idx_others[idx_others != idx_var_names]
  state_order <- state_order[, c(idx_var_names, idx_others), drop=FALSE]
  #
  # Check variables in state_order not in data
  #
  mismatch <- is.na (match (state_order$var_names, data_var_names))
  if ( any (mismatch) ) {
    msg_str <- list_to_str (state_order$var_names[mismatch], n_max=10)
    if (sum (mismatch) == 1)
      miic_warning ("state order", "the variable ", msg_str,
        " does not match any name in input data and will be ignored.")
    else
      miic_warning ("state order", sum (mismatch), " variables (", msg_str,
        ") do not match any name in input data and will be ignored.")
    state_order <- state_order[!mismatch, ]
  }
  #
  # Before checking variables in data not in the state_order
  # if var_type, is_contextual or is_consequence are present, we flag NA
  # in these columns as "NA" to be able to display correct warnings later.
  # ( !! this change the column type to character, even if no NA is detected !! )
  #
  if ("var_type" %in% colnames (state_order) )
    state_order$var_type[ is.na (state_order$var_type) ] <- "NA"
  if ("is_contextual" %in% colnames (state_order) )
    state_order$is_contextual[ is.na (state_order$is_contextual) ] <- "NA"
  if ("is_consequence" %in% colnames (state_order) )
    state_order$is_consequence[ is.na (state_order$is_consequence) ] <- "NA"
  #
  # Check variable in data not in the state_order
  #
  not_found <- is.na (match (data_var_names, state_order$var_names))
  if ( any (not_found) ) {
    msg_str <- list_to_str (data_var_names[not_found], n_max=10)
    if ( sum (not_found) == 1)
      miic_warning ("state order", "the variables ", msg_str,
        " in input data can not be found in the state order. Default values",
        " will be used for this variable.")
    else
      miic_warning ("state order", sum (not_found), " variables (", msg_str,
        ") in input data can not be found in the state order. Default values",
        " will be used for these variables.")
    #
    # Add missing variable names with NA in the other columns
    #
    na_vals <- rep (NA, ncol(state_order) - 1)
    for (i in which (not_found))
      state_order[nrow(state_order)+1,] <- c(data_var_names[i], na_vals)
  }
  #
  # The state_order rows are ordered as the variables in the data
  #
  state_order <- state_order [order (match(state_order$var_names, data_var_names)),,
                              drop=FALSE]
  rownames (state_order) <- NULL
  #
  # var_type (0=discrete / 1=continuous)
  #
  data_is_num <- unlist (lapply (input_data, is.numeric) )
  n_unique_vals <- unlist (lapply (input_data, function (x) {
    length (unique (x[!is.na(x)] ) ) } ) )
  var_type_specified <- rep (F, n_vars)
  if ( ! ("var_type" %in% colnames (state_order) ) ) {
    state_order$var_type <- as.integer (data_is_num)
    #
    # Continuous Variables with less than MIIC_CONTINUOUS_TRESHOLD are
    # considered as discrete
    #
    state_order$var_type [ n_unique_vals < MIIC_CONTINUOUS_TRESHOLD ] = 0
  } else {
    var_type_specified <- rep (T, n_vars)
    #
    # Exclude NAs from the warning (NA = row added because var name missing)
    #
    non_valid <- ( ( ! (is.na (state_order$var_type) ) )
                 & ( ! (state_order$var_type %in% c(0,1)) ) )
    if ( any (non_valid) ) {
      msg_str <- list_to_str (state_order$var_names[non_valid], n_max=10)
      if ( sum (non_valid) == 1)
        miic_warning ("state order", "the variable ", msg_str,
          " does not have a valid value in the var_type column,",
          " the invalid value be ignored and type determined from data.")
      else
        miic_warning ("state order", sum(non_valid), " variables (", msg_str,
          ") do not have a valid value in the var_type column,",
          " the invalid values will be ignored and types determined from data.")
    }
    #
    # All non 0 or 1 need to be fixed
    #
    non_valid <- ! (state_order$var_type %in% c(0,1))
    if ( any (non_valid) ) {
      state_order$var_type[non_valid] <- as.integer(data_is_num)[non_valid]
      var_type_specified[non_valid] <- F
    }
    #
    # Ensure the type of var_type is numerical
    # (because when looking for NAs present before, the column type has been
    # shifted to character. Now, we are sure that we have only O and 1 => as.int
    #
    state_order$var_type = as.integer(state_order$var_type)
    #
    # Check var_type against data
    #
    pb_continuous <- (state_order$var_type == 1) & (!data_is_num)
    if ( any (pb_continuous) ) {
      msg_str <- list_to_str (state_order$var_names[pb_continuous], n_max=10)
      if ( sum (pb_continuous) == 1)
        miic_warning ("state order", "the variable ", msg_str,
          " is declared continuous in the var_type column but is not numeric.",
          " This variable will be considered as discrete.")
      else
        miic_warning ("state order", sum (pb_continuous), " variables (", msg_str,
          ") are declared continuous in the var_type column but these variables",
          " are not numeric. These variables will be considered as discrete.")
      state_order$var_type[pb_continuous] <- 0
    }
  }
  #
  # Check the number of unique values versus var_type
  #
  for (i in 1:n_vars) {
    if (state_order[i, "var_type"] == 1) { # Continuous
      #
      # Less than 3 unique values does not make sense for a continuous variable
      #
      if (n_unique_vals[[i]] <= 2) {
        if (var_type_specified[[i]])
          miic_warning ("state order", "variable ", data_var_names[[i]],
              " specified as continuous has only ", n_unique_vals[[i]],
              " non-NA unique values. It will be processed as discrete.")
        state_order$var_type[[i]] <- 0
      }
      #
      # Less than MIIC_CONTINUOUS_TRESHOLD unique variables can be discretized
      # but may not be truly continuous
      #
      else if (n_unique_vals[[i]] < MIIC_CONTINUOUS_TRESHOLD)
        miic_warning ("state order", "numerical variable ", data_var_names[[i]],
                 " is treated as continuous but has only ", n_unique_vals[[i]],
                 " non-NA unique values.")
    } else { # discrete var
      if ( data_is_num[[i]] && (n_unique_vals[[i]] >= MIIC_CONTINUOUS_TRESHOLD * 2) )
        miic_warning ("state order", "numerical variable ", data_var_names[[i]],
          " is treated as discrete but has ", n_unique_vals[[i]], " levels.")
    }
  }
  #
  # is_contextual
  #
  if ( ! ("is_contextual" %in% colnames (state_order) ) )
    state_order$is_contextual <- rep (0, n_vars)
  else {
    #
    # Exclude NAs from the warning (NA = row added because var name missing)
    #
    non_valid <- ( ( ! (is.na (state_order$is_contextual) ) )
                 & ( ! (state_order$is_contextual %in% c(0,1)) ) )
    if (any (non_valid)) {
      msg_str <- list_to_str (state_order$var_names[non_valid], n_max=10)
      if (sum (non_valid) == 1)
        miic_warning ("state order", "the variable ", msg_str,
          " does not have a valid value in the is_contextual column,",
          " this variable will be considered as not contextual.")
      else
        miic_warning ("state order", sum (non_valid), " variables (", msg_str,
          ") do not have a valid value in the is_contextual column,",
          " these variables will be considered as not contextual.")
    }
    #
    # All non 0 or 1 are not valid => set to not contextual
    #
    non_valid <- ! (state_order$is_contextual %in% c(0,1))
    if (any (non_valid))
      state_order$is_contextual[non_valid] <- 0
    #
    # Ensure the type of is_contextual is numerical
    # (because when looking for NAs present before, the column type has been
    # shifted to character. Now, we are sure that we have only O and 1 => as.int
    #
    state_order$is_contextual = as.integer(state_order$is_contextual)
    #
    # Stop if all variables are contextual
    #
    if (all (state_order$is_contextual == 1))
      stop (paste0 ("All variables have been defined as contextual in the",
        " state_order. No network can be infered with these settings."), call.=FALSE)
  }
  #
  # is_consequence
  #
  if ( ! ("is_consequence" %in% colnames (state_order) ) )
    state_order$is_consequence <- rep (0, n_vars)
  else {
    #
    # Exclude NAs from warnings (NA = row added because var name missing)
    #
    non_valid <- ( ( ! (is.na (state_order$is_consequence) ) )
                 & ( ! (state_order$is_consequence %in% c(0,1)) ) )
    if (any (non_valid)) {
      msg_str <- list_to_str (state_order$var_names[non_valid], n_max=10)
      if (sum (non_valid) == 1)
        miic_warning ("state order", "the variable ", msg_str,
          " does not have a valid value in the is_consequence column,",
          " this variable will be considered as not consequence")
      else
        miic_warning ("state order", sum (non_valid), " variables (", msg_str,
          ") do not have a valid value in the is_consequence column,",
          " these variables will be considered as not consequence")
    }
    #
    # All non 0 or 1 are not valid => set to not consequence
    #
    non_valid <- ! (state_order$is_consequence %in% c(0,1))
    if (any (non_valid))
      state_order$is_consequence[non_valid] <- 0
    #
    # Ensure the type of is_consequence is numerical
    # (because when looking for NAs present before, the column type has been
    # shifted to character. Now, we are sure that we have only O and 1 => as.int
    #
    state_order$is_consequence = as.integer(state_order$is_consequence)
    #
    # Stop if all variables are consequences
    #
    if (all (state_order$is_consequence == 1))
      stop (paste0 ("All variables have been defined as consequences in the",
        " state_order. No network can be infered with these settings."), call.=FALSE)
  }
  #
  # levels_increasing_order
  #
  if ( ! ("levels_increasing_order" %in% colnames (state_order) ) )
    state_order$levels_increasing_order <- NA
  else {
    for (i in 1:n_vars) {
      order_str <- state_order[i, "levels_increasing_order"]
      if ( is.na (order_str) )
        next
      if (order_str == "") {
        state_order[i, "levels_increasing_order"] <- NA
        next
      }
      if (state_order[i, "var_type"] == 1) {
        miic_warning ("state order", "variable ", state_order[i, "var_names"],
                      " is considered as a continuous variable,",
                      " the provided levels order will be ignored.")
        state_order[i, "levels_increasing_order"] <- NA
        next
      }
      #
      # Discrete var, check the match of unique values in data and values
      # in levels_increasing_order
      #
      orders <- trimws (unlist (strsplit (as.character (order_str), ",") ) )
      values <- unique (input_data[!is.na(input_data[,i]),i])
      #
      # Convert values in state order using the same type as data.
      # It will avoid issues when comparing TRUE/FALSE with T/F or 1.0 with 1
      # If the values comming from the state_order can not be converted,
      # leave the value unchanged to display a meaningful warning laterly
      #
      if (is.logical (values)) {
        suppressWarnings ( { orders_log <- as.logical(orders) } )
        orders[!is.na (orders_log)] <- orders_log[!is.na (orders_log)]
      }
      else if (is.integer (values)) {
        suppressWarnings ( { orders_int <- as.integer(orders) } )
        orders[!is.na (orders_int)] <- orders_int[!is.na (orders_int)]
      }
      else if (is.numeric (values)) {
        suppressWarnings ( { orders_num <- as.numeric(orders) } )
        orders[!is.na (orders_num)] <- orders_num[!is.na (orders_num)]
      }
      orders <- as.character(orders)
      values <- as.character (values)
      #
      # If the column in input_data does not contain a "NA" string
      # and if the levels_increasing_order contains "NA" string,
      # issue a specific warning about NA
      # NB : we test here only "NA" string as only "NA" is converted as NA in R
      # by default when using read.table or read.csv, other #NA, N/A, ...
      # will however be discarded later with a less specific warning
      #
      if ( (! ("NA" %in% values)) && ("NA" %in% orders) ) {
        miic_warning ("state order", "variable ", state_order[i, "var_names"],
          " has a NA value in the provided levels order. NA can not be used to",
          " order levels and should not be included in the provided levels order.")
        orders <- orders[ orders != "NA"]
        if ( length (orders) == 0 ) {
          state_order[i, "levels_increasing_order"] <- NA
          next
        }
      }
      #
      # Check if some provided levels are not in the data
      #
      not_in_data <- is.na (match (orders, values) )
      if ( any (not_in_data) ) {
        msg_str <- list_to_str (orders[not_in_data], n_max=10)
        if (sum (not_in_data) == 1)
          miic_warning ("state order", "variable ", state_order[i, "var_names"],
            " has value ", msg_str, " in the provided levels order not present",
            " in the data. This value will be ignored.")
        else
          miic_warning ("state order", "variable ", state_order[i, "var_names"],
            " has values ", msg_str, " in the provided levels order not present",
            " in the data. These values will be ignored.")
        orders <- orders[!not_in_data]
        if ( length (orders) == 0 ) {
          state_order[i, "levels_increasing_order"] <- NA
          next
        }
      }
      #
      # Check if missing levels compared to data
      #
      absent <- is.na (match (values, orders) )
      if ( any (absent) ) {
        msg_str <- list_to_str (values[absent], n_max=10)
        if (sum (absent) == 1)
          miic_warning ("state order", "variable ", state_order[i, "var_names"],
            " has value ", msg_str, " in the data that can not be found",
            " in the provided levels order.",
            " The provided levels order for this variable will be ignored.")
        else
          miic_warning ("state order", "variable ", state_order[i, "var_names"],
            " has values ", msg_str, " in the data that can not be found",
            " in the provided levels order.",
            " The provided levels order for this variable will be ignored.")
        state_order[i, "levels_increasing_order"] <- NA
        next
      }
      #
      # If the levels_increasing_order was not turned into NA,
      # update the levels_increasing_order to have a clean string without
      # leading or trailing blanks and same type format between data and state
      # order (i.e. if data column is logical (TRUE/FALSE), values (T/F) in the
      # state_order will be converted as TRUE/FALSE )
      #
      state_order[i, "levels_increasing_order"] <- paste0 (orders, collapse=",")
    }
  }
  #
  # Cross checks : check that no var is both contextual and consequence
  #
  ctx_and_csq = state_order$is_contextual + state_order$is_consequence
  ctx_and_csq = (ctx_and_csq >= 2)
  if (any (ctx_and_csq)) {
    msg_str <- list_to_str (state_order$var_names[ctx_and_csq], n_max=10)
    if (sum (ctx_and_csq) == 1)
      miic_warning ("state order", "the variable ", msg_str,
        " can not be defined as both contextual and consequence. This variable",
        " will be considered as neither contextual nor consequence.")
    else
      miic_warning ("state order", sum (ctx_and_csq), " variables (", msg_str,
        ") can not be defined as both contextual and consequence. These",
        " variables will be considered as neither contextual nor consequence.")
    state_order$is_contextual[ctx_and_csq] = 0
    state_order$is_consequence[ctx_and_csq] = 0
  }
  return (state_order)
}

#-------------------------------------------------------------------------------
# check_other_df
#-------------------------------------------------------------------------------
# input_data: a dataframe with variables as columns and rows as samples
# - df: the datafame to check, expected to be a 2 columns dataframe. All values
#   in the dataframe are expected to be variables names.
#   An invalid dataframe will be ignored, Invalid rows will be discarded
# - df_name: the datafame name (i.e. :"black box", "true edges")
#   This value is used only to display messages
# return: the dataframe checked
#-------------------------------------------------------------------------------
check_other_df <- function (input_data, df, df_name) {
  if ( is.null(df) )
    return (NULL)
  #
  # Basic checks
  #
  if ( ! is.data.frame(df) ) {
    miic_warning (df_name, "The ", df_name, " parameter, if provided,",
      " must be a dataframe. The ", df_name, " will be ignored.")
    return (NULL)
  }

  n_cols <- 2
  if (ncol(df) != n_cols) {
    miic_warning (df_name, "The expected dataframe must have ", n_cols,
      " columns but the provided one has ", ncol(df), " and will be ignored.")
    return (NULL)
  }

  if (nrow(df) == 0) {
    miic_warning (df_name, "The provided dataframe is empty.")
    return (df)
  }

  data_var_names <- colnames (input_data)
  rows_with_warning <- c()
  for ( row_idx in 1:nrow(df) ) {
    for (col_idx in 1:2) {
      one_var_name <- df[row_idx, col_idx]
      if (! (one_var_name %in% data_var_names) ) {
        miic_warning (df_name, "The variable ", one_var_name,
          " is not present in the input data. The row ", row_idx, " will be ignored.")
        rows_with_warning[[length(rows_with_warning)+1]] <- row_idx
      }
    }
    if (  (length(rows_with_warning) > 0)
       && (rows_with_warning[[length(rows_with_warning)]] == row_idx) )
      next
    if (df[row_idx, 1] == df[row_idx, 2]) {
      miic_warning (df_name, "the variables must be different for each row (found ",
        df[row_idx, 1], " two times at row ", row_idx, "). This row will be ignored.")
      rows_with_warning[[length(rows_with_warning)+1]] <- row_idx
    }
  }
  rows_ok <- unlist (lapply (1:nrow(df), FUN=function (x) { ! (x %in% rows_with_warning) } ) )
  df <- df [rows_ok, ]
  #
  # Remove duplicate row
  n_rows_sav = nrow(df)
  # Equal rows
  df = unique (df)
  rownames(df) = NULL
  # Equal rows but with variable names swapped
  rows_kept = rep (T, nrow(df))
  for (i in 1:nrow(df)) {
    if ( ! rows_kept[[i]] )
      next
    dup_inverse = ( (df[,1] == df[i,2])
                  & (df[,2] == df[i,1])
                  & (rownames(df) != i) )
    rows_kept = rows_kept & (!dup_inverse)
  }
  df = df[rows_kept,]
  if ( n_rows_sav != nrow(df) ) {
    if (n_rows_sav - nrow(df) == 1)
      miic_warning (df_name, "1 row is duplicated. Only one instance",
        " of the row will be used.")
    else
      miic_warning (df_name, n_rows_sav - nrow(df), " rows are duplicated.",
        " Only one instance of these rows will be used.")
  }

  if (nrow(df) == 0)
    miic_warning (df_name, "The provided dataframe is empty.")
  return (df)
}

#-------------------------------------------------------------------------------
# Check_parameters
#-------------------------------------------------------------------------------
# Check all input parameters that are not df.
# Params : input_data + all possible parameters of miic method
# Returns: a list with all the parameters, eventually modified or initialized
#-------------------------------------------------------------------------------
check_parameters <- function (input_data, n_threads, cplx,
  orientation, ori_proba_ratio, ori_consensus_ratio, propagation, latent,
  n_eff, n_shuffles, conf_threshold, sample_weights, test_mar,
  consistent, max_iteration, consensus_threshold, negative_info, verbose) {
  list_ret = list()
  if ( is.null (n_threads) || (!is.numeric(n_threads))
       || (round(n_threads,0) != n_threads) || (n_threads < 1) ) {
    miic_warning ("parameters", "supplied value ", n_threads,
      " for the n_threads parameter is invalid. It must be an integer >= 1.",
      " The default value (1) will be used.")
    n_threads = 1
  }
  list_ret$n_threads = n_threads

  if (  is.null (cplx)
     || ( ! is.character(cplx) )
     || ( ! (cplx %in% c("nml", "mdl") ) ) ) {
    miic_warning ("parameters", "supplied value ", cplx,
      " for the complexity parameter is invalid. It must be 'nml' or 'mdl'.",
      " The default value ('nml') will be used.")
    cplx = "nml"
  }
  list_ret$cplx = cplx

  if ( is.null (orientation) || (!is.logical(orientation)) ) {
    miic_warning ("parameters", "supplied value ", orientation,
      " for the orientation parameter is invalid. It must be TRUE/FALSE.",
      " The default value (TRUE) will be used.")
    orientation = TRUE
  }
  list_ret$orientation = orientation

  if (  is.null (ori_proba_ratio) || (!is.numeric(ori_proba_ratio))
     || (ori_proba_ratio < 0) || (ori_proba_ratio > 1) ) {
    miic_warning ("parameters", "supplied value ", ori_proba_ratio,
      " for the orientation probabilty ratio parameter is invalid.",
      " It must be a floating point between 0 and 1.",
      " The default value (1) will be used.")
    ori_proba_ratio = 1
  }
  list_ret$ori_proba_ratio = ori_proba_ratio

  if ( is.null (ori_consensus_ratio) )
    ori_consensus_ratio = ori_proba_ratio
  else if (  (!is.numeric (ori_consensus_ratio))
          || (ori_consensus_ratio < 0) || (ori_consensus_ratio > 1) ) {
    miic_warning ("parameters", "supplied value ", ori_consensus_ratio,
      " for the orientation concensus ratio parameter is invalid.",
      " It must be a floating point between 0 and 1.",
      " The default value (same as orientation probabilty ratio: ",
      ori_proba_ratio, ") will be used.")
    ori_consensus_ratio = ori_proba_ratio
  }
  list_ret$ori_consensus_ratio = ori_consensus_ratio

  if ( is.null (propagation) || (!is.logical(propagation)) ) {
    miic_warning ("parameters", "supplied value ", propagation,
      " for the propagation parameter is invalid. It must be TRUE/FALSE.",
      " The default value (FALSE) will be used.")
    propagation = FALSE
  }
  list_ret$propagation = propagation

  if ( is.null(latent) || (!(latent %in% MIIC_VALID_LATENT)) ) {
    msg_str = paste (MIIC_VALID_LATENT, collapse=", ")
    miic_warning ("parameters", "supplied value ", latent,
      " for the latent parameter is invalid. Possible values are: ", msg_str,
      ". The default value (", MIIC_VALID_LATENT[[1]], ") will be used.")
    latent = MIIC_VALID_LATENT[[1]]
  }
  list_ret$latent = latent

  if (  is.null(n_eff) || (!is.numeric(n_eff)) || (round(n_eff,0) != n_eff)
     || ((n_eff <= 0) && (n_eff != -1)) || (n_eff > nrow(input_data)) ) {
    miic_warning ("parameters", "supplied value ", n_eff,
      " for the number of effective samples is invalid.",
      " The number of effective samples must be an integer that can be -1",
      " for an automatic assignment or a positive number less or",
      " equal to the number of samples. The default value (-1) will be used.")
    n_eff = -1
  }
  list_ret$n_eff = n_eff

  if (  is.null(n_shuffles) || (!is.numeric(n_shuffles))
     || (round(n_shuffles,0) != n_shuffles) || (n_shuffles < 0) ) {
    miic_warning ("parameters", "supplied value ", n_shuffles,
      " for the number of shufflings is invalid. It must be an integer >= 0.",
      " The default value (0) will be  used.")
    n_shuffles = 0
  }

  if (n_shuffles == 0) {
    if ( (!is.null (conf_threshold)) && (conf_threshold != 0) )
      miic_warning ("parameters", "supplied value ", conf_threshold,
        " for the confidence threshold parameter will be ignored",
        " as the number of shufflings is set to 0.",
        " To activate the confidencence cut, both the number of shufflings",
        " and conf_threshold = 0.01).")
    conf_threshold = 0
  } else {
    if (  is.null (conf_threshold)
       || (!is.numeric(conf_threshold))
       || (conf_threshold < 0) ) {
      miic_warning ("parameters", "supplied value ", conf_threshold,
        " for the confidence threshold parameter is invalid.",
        " When confidence cut is activated (when n_shuffles > 0),",
        " the confidence threshold must be a floating point > 0. The",
        " confidence cut will be desactivated and default values will be used",
        " for the number of shufflings (0) and the confidence threshold (0).")
      n_shuffles = 0
      conf_threshold = 0
    } else if (conf_threshold == 0) {
      miic_warning ("parameters", "the confidence threshold parameter is 0",
        " but it must be > 0 when confidence cut is activated",
        " (when n_shuffles > 0). The confidence cut will be desactivated.",
        " To activate the confidencence cut, both the number of shufflings",
        " and the confidence threshold must be > 0 (i.e.: n_shuffles = 100",
        " and conf_threshold = 0.01).")
      n_shuffles = 0
    }
  }
  list_ret$n_shuffles = n_shuffles
  list_ret$conf_threshold = conf_threshold

  if (  ( ! is.null (sample_weights) )
     && (  (length(sample_weights) != nrow(input_data))
        || (any(!is.numeric(sample_weights)))
        || (any(sample_weights < 0))
        || (any(sample_weights > 1)) ) ) {
    miic_warning ("parameters", "supplied value for the sample_weights parameter",
      " is invalid. It must be a vector of the same size as the number of",
      " samples in the input data and all weights must be floating points",
      " in the [0,1] range. The parameter will be ignored.")
    sample_weights = NULL
  }
  list_ret$sample_weights = sample_weights

  if ( is.null (test_mar) || (!is.logical(test_mar)) ) {
    miic_warning ("parameters", "Supplied value ", test_mar,
      " for the missing at random test parameter is invalid.",
      " It must be TRUE/FALSE. The default value (TRUE) will be used.")
    test_mar = TRUE
  }
  list_ret$test_mar = test_mar

  if ( is.null(consistent) || (!(consistent %in% MIIC_VALID_CONSISTENT)) ) {
    msg_str = paste (MIIC_VALID_CONSISTENT, collapse=", ")
    miic_warning ("parameters", "supplied value ", consistent,
      " for the consistent parameter is invalid. Possible values are: ", msg_str,
      ". The default value (", MIIC_VALID_CONSISTENT[[1]], ") will be used.")
    consistent = MIIC_VALID_CONSISTENT[[1]]
  }
  list_ret$consistent = consistent

  if (consistent == "no") {
    if ( (!is.null (max_iteration)) && (max_iteration != 100) )
      miic_warning ("parameters", "supplied value ", max_iteration,
        " for the maximum iteration parameter will not be used",
        " as consistency is off.")
    max_iteration = 100
    if ( (!is.null (consensus_threshold)) && (consensus_threshold != 0.8) )
      miic_warning ("parameters", "Supplied value ", consensus_threshold,
        " for the consensus threshold parameter will not be used",
        " as consistency is off.")
    consensus_threshold = 0.8
  } else { # Consistency on
    if (  is.null (max_iteration) || (!is.numeric(max_iteration))
       || (round(max_iteration,0) != max_iteration) || (max_iteration <= 0) ) {
      miic_warning ("parameters", "supplied value ", max_iteration,
        " for the maximum iteration parameter is invalid.",
        " It must be a stricly positive integer when consistency is activated.",
        " The default value (100) will be used.")
      max_iteration = 100
    }
    if (  is.null (consensus_threshold) || (!is.numeric(consensus_threshold))
       || (consensus_threshold < 0.5) || (consensus_threshold > 1) ) {
      miic_warning ("parameters", "supplied value ", consensus_threshold,
        " for the consensus threshold parameter is invalid.",
        " It must be a floating point between 0.5 and 1 when consistency is",
        " activated. The default value (0.8) will be used.")
      consensus_threshold = 0.8
    }
  }
  list_ret$max_iteration = max_iteration
  list_ret$consensus_threshold = consensus_threshold

  if ( is.null (negative_info) || (!is.logical(negative_info)) ) {
    miic_warning ("parameters", "supplied value ", negative_info, " for",
      " parameter allowing/disallowing negative shifted mutual information is",
      " invalid. It must be TRUE/FALSE. The default value (FALSE) will be used.")
    negative_info = FALSE
  }
  list_ret$negative_info = negative_info

  if ( is.null (verbose) || (!is.logical(verbose)) ) {
    miic_warning ("parameters", "supplied value ", verbose,
      " for the verbose parameter is invalid. It must be TRUE/FALSE.",
      " The default value (FALSE) will be used.")
    verbose = FALSE
  }
  list_ret$verbose = verbose
  return (list_ret)
}
