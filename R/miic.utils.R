#*******************************************************************************
# Filename   : miic.utils.R
#
# Description: various utilities functions and constants for miic
#*******************************************************************************

#===============================================================================
# CONSTANTS
#===============================================================================
MIIC_TEMPORAL_MODES = c ("TS", "TNS")
MIIC_VALID_MODES = c ("S", MIIC_TEMPORAL_MODES, "L")

MIIC_VALID_LATENT <- c ("orientation", "yes", "no")
MIIC_VALID_CONSISTENT <- c ("no", "orientation", "skeleton")

MIIC_CONTINUOUS_TRESHOLD <- 5

STATE_ORDER_STANDARD_VALID_COLUMS <- c ("var_names", "var_type",
    "levels_increasing_order", "is_contextual", "is_consequence",
    "group", "group_color")
STATE_ORDER_TEMPORAL_VALID_COLUMNS = c (STATE_ORDER_STANDARD_VALID_COLUMS,
                                        "n_layers", "delta_t", "mov_avg")

#===============================================================================
# FUNCTIONS
#===============================================================================
# list_to_str
#-------------------------------------------------------------------------------
# Utility function to transform the first n_max items of a list into a
# displayable (comma + space separated) string: "item1, item2, item3, ..."
# Parameters:
# - list: a list
# - n_max: int, optional, NULL by default, maximum of items used. If NULL,
#   all items are used. If the list has more than n_max items, ", ..." is added
# Return: string, the displayable string
#-------------------------------------------------------------------------------
list_to_str <- function (a_list, n_max=NULL)
  {
  if (is.null (a_list) )
    return ("NULL")
  if (length (a_list) == 0)
    return ("")
  if ( (! is.null (n_max)) && (length (a_list) > n_max) )
    a_list <- c (a_list[1:n_max], "...")
  ret <- paste (unlist(a_list), collapse=", ")
  return (ret)
  }

#-------------------------------------------------------------------------------
# miic_error
#-------------------------------------------------------------------------------
# Utility function to raise an error and stop
#-------------------------------------------------------------------------------
miic_error <- function (context, ...)
  {
  stop (paste0 ("Error in ", context, ": ", ...), call.=FALSE)
  }

#-------------------------------------------------------------------------------
# miic_warning
#-------------------------------------------------------------------------------
# Utility function to raise a warning
#-------------------------------------------------------------------------------
miic_warning <- function (context, ...)
  {
  warning (paste0 ("Warning in ", context, ": ", ...), call.=FALSE)
  }

#-------------------------------------------------------------------------------
# miic_msg
#-------------------------------------------------------------------------------
# Utility function to display a message
#-------------------------------------------------------------------------------
miic_msg <- function (...)
  {
  cat (paste0 (..., "\n") )
  }

#-------------------------------------------------------------------------------
# Check input data
#-------------------------------------------------------------------------------
# Check applied;
# - data is a non empty data frame
# - columns full of NA are discarded
# - rows full of NA are discarded in non temporal modes
# - warn about constants variables
#
# Parameters:
# - input_data: a data frame with variables as columns and rows as samples
# - mode : MIIC mode
# Return:
# - input_data: the input data data frame, eventually without full of NAs rows
#-------------------------------------------------------------------------------
check_input_data <- function (input_data, mode)
  {
  if ( is.null(input_data) )
    miic_error ("input data", "The input data is required.")
  if ( ! is.data.frame (input_data) )
    miic_error ("input data", "The input data must be a dataframe.")
  #
  # Ensure we have a true data frame (i.e.: not a tibble) with data
  #
  input_data = as.data.frame (input_data)
  if (nrow (input_data) == 0)
    miic_error ("input data", "The input data is empty.")
  if (ncol (input_data) == 0)
    miic_error ("input data", "The input data has no variable.")
  if ( (mode %in% MIIC_TEMPORAL_MODES) && (ncol (input_data) == 1) )
    miic_error ("input data", "The input data has only one column (assumed to be the time steps), no variable.")
  #
  # Check variables full of NAs
  #
  cols_only_na <- colSums (is.na (input_data)) == nrow (input_data)
  input_data <- input_data [, !cols_only_na, drop=F]
  if (ncol (input_data) == 0)
    miic_error ("input data", "The input data contains only NAs.")
  if ( (mode %in% MIIC_TEMPORAL_MODES) && (ncol (input_data) == 1) )
    miic_error ("input data", "The input data has only one column (assumed to be the timesteps) not full of NAs.")
  if ( any (cols_only_na) )
    miic_warning ("input data", "the input data contains ", sum(cols_only_na),
             " variables(s) with only NAs. These variables(s) will be removed.")
  #
  # Check about rows with only NAs (only if mode is not temporal, as NAs check
  # in temporal mode will be done after data lagging)
  #
  if ( ! (mode %in% MIIC_TEMPORAL_MODES) )
    {
    rows_only_na <- rowSums (is.na (input_data)) == ncol (input_data)
    input_data <- input_data [!rows_only_na, , F]
    if ( any (rows_only_na) )
      miic_warning ("input data", "the input data contains ", sum(rows_only_na),
               " row(s) with only NAs. These row(s) will be removed.")
    }
  #
  # Check constant variables
  #
  if (mode %in% MIIC_TEMPORAL_MODES)
    {
    col_names <- colnames (input_data[, 2:ncol(input_data), drop=F])
    n_unique_vals <- unlist (lapply (input_data[, 2:ncol(input_data), drop=F],
                          function (x) { length (unique (x[!is.na(x)] ) ) } ) )
    }
  else
    {
    col_names <- colnames (input_data)
    n_unique_vals <- unlist (lapply (input_data,
                          function (x) { length (unique (x[!is.na(x)] ) ) } ) )
    }
  vars_constant <- (n_unique_vals <= 1)
  if ( any (vars_constant) )
    {
    msg_str <- list_to_str (col_names[vars_constant], n_max=10)
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
# Basically, ensure, if a state_order is supplied, that all the values are valid
# and set the wrong ones to default values (or ignore them if it can't be fixed)
# Missing information will be added by examining the input_data.
# The returned state_order will have all expected columns for the mode chosen
# and a var_names column with all variables ordered in the same way as the data.
# Few errors lead to a stop : if all variables are contextual or consequences,
# other issues will raise a warning and be "fixed" by using default values.
#
# Please look at the comments inside the function for the detail of checks done
#
# Parameters:
# - input_data: a data frame with variables as columns and rows as samples
# - state_order: a data frame, can be NULL.
#   possible/expected columns are:
#   * var_names: the list of all variables in the input data
#   * var_type: 0=discrete, 1=continuous (default deduced from the input data,
#     note that numerical variables with less than x unique values will be
#     discrete)
#   * levels_increasing_order: NA for continuous variables. For discrete,
#     can be NA or the full ordered list of the unique values. (default NA)
#   * is_contextual: 0=not contextual, 1=contextual (default 0)
#   * is_consequence: 0=not consequence, 1=consequence (default 0)
#   additional possible columns in temporal mode are:
#   * n_layers: the number of layers in the time unfolded graph
#   * delta_t: the number of time steps between layers
#   * mov_avg: if a moving average must applied on some variables
#   NB: is_consequence is not allowed in temporal mode
# - mode: the MIIC mode
# Return: the checked and eventually generated or completed state order
#-------------------------------------------------------------------------------
check_state_order <- function (input_data, state_order, mode)
  {
  if (mode %in% MIIC_TEMPORAL_MODES)
    input_data <- input_data[, 2:ncol(input_data), drop=F]
  data_var_names <- colnames (input_data)
  n_vars <- length (data_var_names)
  #
  # Basic checks
  #
  if ( is.null (state_order) )
    state_order <- data.frame ("var_names"=data_var_names, stringsAsFactors=F)
  if ( ! is.data.frame(state_order) )
    {
    miic_warning ("state order",
      "the supplied state_order is not a data frame and will be ignored.")
    state_order <- data.frame ("var_names"=data_var_names, stringsAsFactors=F)
    }
  else
    #
    # Ensure we have a true data frame (i.e.: not a tibble)
    #
    state_order <- as.data.frame (state_order)
  #
  # Factors lead to wrong test results
  #
  factor_cols <- which (unlist (lapply (state_order, is.factor) ) )
  for (i in factor_cols)
    state_order[,i] <- as.character (state_order[,i])
  #
  # Check content
  #
  if ( ! ("var_names" %in% colnames (state_order)) )
    {
    miic_warning ("state order", "the column var_names is missing,",
                  " the supplied state_order will be ignored.")
    state_order <- data.frame ("var_names"=data_var_names, stringsAsFactors=F)
    }
  #
  # Check if the state_order columns are valid
  #
  if (mode %in% MIIC_TEMPORAL_MODES)
    valid_cols <- STATE_ORDER_TEMPORAL_VALID_COLUMNS
  else
    valid_cols <- STATE_ORDER_STANDARD_VALID_COLUMS
  mismatch <- is.na (match (colnames (state_order), valid_cols))
  if ( any (mismatch) )
    {
    msg_str <- list_to_str (colnames (state_order)[mismatch], n_max=10)
    if (sum (mismatch) == 1)
      miic_warning ("state order", "the column ", msg_str,
        " is not valid and will be ignored.")
    else
      miic_warning ("state order", sum (mismatch), " columns (", msg_str,
        ") are not valid and will be ignored.")
    state_order <- state_order[, !mismatch, drop=F]
    }
  #
  # We ensure that the var_names column is the first
  #
  idx_var_names <- which (colnames(state_order) == "var_names")
  idx_others <- 1:ncol (state_order)
  idx_others <- idx_others[idx_others != idx_var_names]
  state_order <- state_order[, c(idx_var_names, idx_others), drop=F]
  #
  # Check variables in state_order not in data
  #
  mismatch <- is.na (match (state_order$var_names, data_var_names))
  if ( any (mismatch) )
    {
    msg_str <- list_to_str (state_order$var_names[mismatch], n_max=10)
    if (sum (mismatch) == 1)
      miic_warning ("state order", "the variable ", msg_str,
        " does not match any name in input data and will be ignored.")
    else
      miic_warning ("state order", sum (mismatch), " variables (", msg_str,
        ") do not match any name in input data and will be ignored.")
    state_order <- state_order[!mismatch, , drop=F]
    }
  #
  # Before checking variables in data not in the state_order
  # if var_type, is_contextual or is_consequence are present, we flag NA
  # in these columns as "NA" string (that is not a NA value)
  # to be able to display correct warnings later.
  # The same applies for the specific columns of the temporal modes.
  # ( !! this changes the column type to character, even if no NA is detected !! )
  #
  if ("var_type" %in% colnames (state_order) )
    state_order$var_type[ is.na (state_order$var_type) ] <- "NA"
  if ("is_contextual" %in% colnames (state_order) )
    state_order$is_contextual[ is.na (state_order$is_contextual) ] <- "NA"
  if ("is_consequence" %in% colnames (state_order) )
    state_order$is_consequence[ is.na (state_order$is_consequence) ] <- "NA"
  if (mode %in% MIIC_TEMPORAL_MODES)
    {
    if ("n_layers" %in% colnames (state_order) )
      state_order$n_layers[ is.na (state_order$n_layers) ] <- "NA"
    if ("delta_t" %in% colnames (state_order) )
      state_order$delta_t[ is.na (state_order$delta_t) ] <- "NA"
    if ("mov_avg" %in% colnames (state_order) )
      state_order$mov_avg[ is.na (state_order$mov_avg) ] <- "NA"
    }
  #
  # Check variables in data not in the state_order
  #
  not_found <- is.na (match (data_var_names, state_order$var_names))
  if ( any (not_found) )
    {
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
  state_order <- state_order [order (match(state_order$var_names,
                                           data_var_names)), , drop=F]
  rownames (state_order) <- NULL
  #
  # var_type (0=discrete / 1=continuous)
  #
  data_is_num <- unlist (lapply (input_data, is.numeric) )
  n_unique_vals <- unlist (lapply (input_data, function (x) {
    length (unique (x[!is.na(x)] ) ) } ) )
  var_type_specified <- rep (F, n_vars)
  if ( ! ("var_type" %in% colnames (state_order) ) )
    {
    state_order$var_type <- as.integer (data_is_num)
    #
    # Continuous Variables with less than MIIC_CONTINUOUS_TRESHOLD are
    # considered as discrete
    #
    state_order$var_type [ n_unique_vals < MIIC_CONTINUOUS_TRESHOLD ] = 0
    state_order$var_type = as.integer(state_order$var_type)
    }
  else
    {
    var_type_specified <- rep (T, n_vars)
    #
    # Exclude NAs from the warning (NA = row added because var name missing)
    #
    non_valid <- ( ( ! (is.na (state_order$var_type) ) )
                 & ( ! (state_order$var_type %in% c(0,1)) ) )
    if ( any (non_valid) )
      {
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
    # All non 0 or 1 need to be fixed (same as if var_type column is not supplied)
    #
    non_valid <- ! (state_order$var_type %in% c(0,1))
    if ( any (non_valid) )
      {
      state_order$var_type[non_valid] <- as.integer(data_is_num)[non_valid]
      state_order$var_type [ non_valid & (n_unique_vals < MIIC_CONTINUOUS_TRESHOLD) ] = 0
      var_type_specified[non_valid] <- F
      }
    #
    # Ensure the type of var_type is integer
    # (because when looking for NAs present before, the column type has been
    # shifted to character. Now, we are sure that we have only O and 1
    #
    state_order$var_type = as.integer(state_order$var_type)
    #
    # Check var_type against data
    #
    pb_continuous <- (state_order$var_type == 1) & (!data_is_num)
    if ( any (pb_continuous) )
      {
      msg_str <- list_to_str (state_order$var_names[pb_continuous], n_max=10)
      if ( sum (pb_continuous) == 1)
        miic_warning ("state order", "the variable ", msg_str,
          " is declared continuous in the var_type column but is not numeric.",
          " This variable will be considered as discrete.")
      else
        miic_warning ("state order", sum (pb_continuous), " variables (", msg_str,
          ") are declared continuous in the var_type column but these variables",
          " are not numeric. These variables will be considered as discrete.")
      state_order$var_type[pb_continuous] <- as.integer (0)
      var_type_specified[pb_continuous] <- F
      }
    }
  #
  # In temporal mode, we store if var_type was specified by the user
  # It will be used when assessing the type of lagged variables to decide
  # if a warning is raised (when user had specified the type)
  # or not (automatic type determination)
  #
  if (mode %in% MIIC_TEMPORAL_MODES)
    state_order$var_type_specified <- var_type_specified
  #
  # Check the number of unique values versus var_type
  #
  for (i in 1:n_vars)
    {
    if (state_order[i, "var_type"] == 1) # Continuous
      {
      # 2 unique values or less does not make sense for a continuous variable
      #
      if (n_unique_vals[[i]] <= 2)
        {
        if (var_type_specified[[i]])
          miic_warning ("state order", "variable ", data_var_names[[i]],
              " specified as continuous has only ", n_unique_vals[[i]],
              " non-NA unique values. It will be processed as discrete.")
        state_order$var_type[[i]] <- as.integer(0)
        if (mode %in% MIIC_TEMPORAL_MODES)
          state_order$var_type_specified[[i]] <- F
        }
      #
      # Less than MIIC_CONTINUOUS_TRESHOLD unique variables can be discretized
      # but may not be truly continuous
      #
      else if (n_unique_vals[[i]] < MIIC_CONTINUOUS_TRESHOLD)
        miic_warning ("state order", "numerical variable ", data_var_names[[i]],
                 " is treated as continuous but has only ", n_unique_vals[[i]],
                 " non-NA unique values.")
      }
    else # discrete var
      {
      if ( data_is_num[[i]] && (n_unique_vals[[i]] >= MIIC_CONTINUOUS_TRESHOLD * 2) )
        miic_warning ("state order", "numerical variable ", data_var_names[[i]],
          " is treated as discrete but has ", n_unique_vals[[i]], " levels.")
      }
    }
  #
  # is_contextual
  #
  if ( ! ("is_contextual" %in% colnames (state_order) ) )
    state_order$is_contextual <- as.integer ( rep (0, n_vars) )
  else
    {
    # Exclude NAs from the warning (NA = row added because var name missing)
    #
    non_valid <- ( ( ! (is.na (state_order$is_contextual) ) )
                 & ( ! (state_order$is_contextual %in% c(0,1)) ) )
    if (any (non_valid))
      {
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
    state_order$is_contextual <- as.integer(state_order$is_contextual)
    #
    # Stop if all variables are contextual
    #
    if (all (state_order$is_contextual == 1))
      miic_error ("state order", "All variables have been defined as",
        " contextual. No network can be infered with these settings.")
    }
  #
  # is_consequence
  #
  if ( ! ("is_consequence" %in% colnames (state_order) ) )
    state_order$is_consequence <- as.integer ( rep (0, n_vars) )
  else
    {
    if (mode %in% MIIC_TEMPORAL_MODES)
      {
      # Exclude NAs from warnings (NA = row added because var name missing)
      # => Look of anything not NA and != 0
      #
      conseq_def <- ( ( ! (is.na (state_order$is_consequence) ) )
                    & ( is.null (state_order$is_consequence)
                      | (state_order$is_consequence != 0) ) )
      if (any (conseq_def))
        {
        msg_str <- list_to_str (state_order$var_names[conseq_def], n_max=10)
        if (sum (conseq_def) == 1)
          miic_warning ("state order", "the variable ", msg_str,
            " is defined as consequence but consequence prior is not compatible",
            " with temporal mode. This variable will be considered as not",
            " consequence.")
        else
          miic_warning ("state order", sum (conseq_def), " variables (", msg_str,
            ") are defined as consequence but consequence prior is not compatible",
            " with temporal mode. These variables will be considered as not",
            " consequence.")
        }
      state_order$is_consequence <- as.integer ( rep (0, n_vars) )
      }
    else # Not temporal mode
      {
      # Exclude NAs from warnings (NA = row added because var name missing)
      #
      non_valid <- ( ( ! (is.na (state_order$is_consequence) ) )
                   & ( ! (state_order$is_consequence %in% c(0,1)) ) )
      if (any (non_valid))
        {
        msg_str <- list_to_str (state_order$var_names[non_valid], n_max=10)
        if (sum (non_valid) == 1)
          miic_warning ("state order", "the variable ", msg_str,
            " does not have a valid value in the is_consequence column,",
            " this variable will be considered as not consequence.")
        else
          miic_warning ("state order", sum (non_valid), " variables (", msg_str,
            ") do not have a valid value in the is_consequence column,",
            " these variables will be considered as not consequence.")
        }
      #
      # All non 0 or 1 are not valid => set to not consequence
      #
      non_valid <- ! (state_order$is_consequence %in% c(0,1))
      if (any (non_valid))
        state_order$is_consequence[non_valid] <- 0
      #
      # Ensure the type of is_consequence is integer
      # (because when looking for NAs present before, the column type has been
      # shifted to character. Now, we are sure that we have only int O and 1
      #
      state_order$is_consequence = as.integer (state_order$is_consequence)
      #
      # Stop if all variables are consequences
      #
      if (all (state_order$is_consequence == 1))
        miic_error ("state order", "All variables have been defined as",
          " consequences. No network can be infered with these settings.")
      }
    }
  #
  # levels_increasing_order
  #
  if ( ! ("levels_increasing_order" %in% colnames (state_order) ) )
    state_order$levels_increasing_order <- NA_character_
  else
    {
    state_order$levels_increasing_order <- as.character (state_order$levels_increasing_order)
    for (i in 1:n_vars)
      {
      order_str <- state_order[i, "levels_increasing_order"]
      if ( is.na (order_str) )
        next
      if (order_str == "")
        {
        state_order[i, "levels_increasing_order"] <- NA_character_
        next
        }
      if (state_order[i, "var_type"] == 1)
        {
        miic_warning ("state order", "variable ", state_order[i, "var_names"],
                      " is considered as a continuous variable,",
                      " the provided levels order will be ignored.")
        state_order[i, "levels_increasing_order"] <- NA_character_
        next
        }
      #
      # The format is a comma separated list
      #
      if (!grepl (",", order_str, fixed=T))
        {
        miic_warning ("state order", "variable ", state_order[i, "var_names"],
                      " has an invalid levels order ", as.character (order_str),
                      ", it should be a comma separated list,",
                      " the provided levels order will be ignored.")
        state_order[i, "levels_increasing_order"] <- NA_character_
        next
        }
      #
      # Discrete var, check the match of unique values in data
      # and values in levels_increasing_order
      #
      orders <- trimws (unlist (strsplit (as.character (order_str), ",") ) )
      values <- unique (input_data[!is.na(input_data[,i]),i])
      #
      # Remove "NA" from levels_increasing_order
      # NB : we test here only "NA" string as only "NA" is converted as NA in R
      # by default when using read.table or read.csv. So a field containing
      # "NA,1,2,3", when split, will contain [ "NA", "1", "2", "3" ]
      #
      if ("NA" %in% orders)
        {
        miic_warning ("state order", "variable ", state_order[i, "var_names"],
          " has NA in the provided levels order. NA can not be used to",
          " order levels and should not be included in the provided levels order.")
        orders <- orders[ orders != "NA"]
        if ( length (orders) == 0 )
          {
          state_order[i, "levels_increasing_order"] <- NA_character_
          next
          }
        }
      #
      # Values in levels_increasing_order must be unique
      #
      duplicate_orders = unique( orders[duplicated(orders)] )
      orders = unique (orders)
      if (length (duplicate_orders) > 0)
        miic_warning ("state order", "variable ", state_order[i, "var_names"],
          " has duplicated values (", list_to_str( duplicate_orders ),
          "). These duplicates will be ignored.")
      #
      # A valid levels order contains at least 2 non NA values
      #
      if (length (orders) <= 1)
        {
        miic_warning ("state order", "the provided levels order ",
                      as.character (order_str), " for variable ",
                      state_order[i, "var_names"],
                      " should contain at least 2 non NA values, it will be ignored.")
        state_order[i, "levels_increasing_order"] <- NA_character_
        next
        }
      #
      # Convert values in state order using the same type as data.
      # It will avoid issues when comparing TRUE/FALSE with T/F or 1.0 with 1
      # If the values coming from the state_order can not be converted,
      # leave the value unchanged to display a warning about the value later
      # i.e.: ["T","A","FALSE"] for a logical column in the data will be
      # converted into [TRUE, NA, FALSE] and then to ["TRUE","A","FALSE"]
      # => "T" will not raise a warning and "A" will
      #
      if (is.logical (values))
        {
        suppressWarnings ( { orders_log <- as.logical(orders) } )
        orders[!is.na (orders_log)] <- orders_log[!is.na (orders_log)]
        }
      else if (is.integer (values))
        {
        suppressWarnings ( { orders_int <- as.integer(orders) } )
        orders[!is.na (orders_int)] <- orders_int[!is.na (orders_int)]
        }
      else if (is.numeric (values))
        {
        suppressWarnings ( { orders_num <- as.numeric(orders) } )
        orders[!is.na (orders_num)] <- orders_num[!is.na (orders_num)]
        }
      orders <- as.character (orders)
      values <- as.character (values)
      #
      # Check if some provided levels are not in the data
      #
      not_in_data <- is.na (match (orders, values) )
      if ( any (not_in_data) )
        {
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
        if ( length (orders) <= 1 )
          {
          state_order[i, "levels_increasing_order"] <- NA_character_
          next
          }
        }
      #
      # Check if missing levels compared to data
      #
      absent <- is.na (match (values, orders) )
      if ( any (absent) )
        {
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
        state_order[i, "levels_increasing_order"] <- NA_character_
        next
        }
      #
      # If the levels_increasing_order was not turned into NA,
      # update the levels_increasing_order to have a clean string without
      # leading or trailing blanks and same type format between data and state
      # order (e.g. if data column is logical (TRUE/FALSE), values (T/F) in the
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
  if (any (ctx_and_csq))
    {
    msg_str <- list_to_str (state_order$var_names[ctx_and_csq], n_max=10)
    if (sum (ctx_and_csq) == 1)
      miic_warning ("state order", "the variable ", msg_str,
        " can not be defined as both contextual and consequence. This variable",
        " will be considered as neither contextual nor consequence.")
    else
      miic_warning ("state order", sum (ctx_and_csq), " variables (", msg_str,
        ") can not be defined as both contextual and consequence. These",
        " variables will be considered as neither contextual nor consequence.")
    state_order$is_contextual[ctx_and_csq] = as.integer (0)
    state_order$is_consequence[ctx_and_csq] = as.integer (0)
    }
  return (state_order)
  }

#-------------------------------------------------------------------------------
# check_other_df
#-------------------------------------------------------------------------------
# Check optional data frame like black box or true edges
# Parameters:
# - input_data: a data frame with variables as columns and rows as samples
# - state_order: the data frame returned by check_state_order. As it has been
#   checked, the state_order should never been NULL or empty. It should also
#   contain an is_contextual column filled with a valid value for each variable
# - df: the data fame to check, expected to be a 2 columns data frame in
#   standard mode and 3 columns data frame in temporal mode.
#   All values in 2 first columns of the data frame are expected to be variables
#   names, and in temporal mode, the 3rd column is expected to contain lags.
#   An invalid data frame will be ignored, Invalid rows will be discarded
# - df_name: the data fame name (i.e. :"black box", "true edges")
#   This value is used only to display messages
# - mode: the MIIC mode
# Return: the data frame checked
#-------------------------------------------------------------------------------
# TODO split in 3 fonctions std, TS, TNS
check_other_df <- function (input_data, state_order, df, df_name, mode)
  {
  if ( is.null(df) )
    return (NULL)
  #
  # Basic checks
  #
  if ( ! is.data.frame(df) )
    {
    miic_warning (df_name, "The ", df_name, " parameter, if provided,",
      " must be a data frame. The ", df_name, " will be ignored.")
    return (NULL)
    }
  #
  # Ensure we have a true data frame (i.e.: not a tibble)
  #
  df <- as.data.frame (df)
  #
  # Check number of cols
  #
  cols_of_variables <- c(1,2)
  if (mode %in% MIIC_TEMPORAL_MODES)
    {
    input_data = input_data[, 2:ncol(input_data), drop=F]
    if (mode == "TS")
      n_cols <- 3
    else
      {
      cols_of_variables = c(1,3)
      n_cols <- 4
      }
    }
  else
    n_cols <- 2
  if (ncol(df) != n_cols)
    {
    miic_warning (df_name, "The expected dataframe must have ", n_cols,
      " columns but the provided one has ", ncol(df), " and will be ignored.")
    return (NULL)
    }
  #
  # Check not empty
  #
  if (nrow(df) == 0)
    {
    miic_warning (df_name, "The provided dataframe is empty.")
    return (NULL)
    }
  #
  # Factors lead to wrong test results
  #
  factor_cols <- which (unlist (lapply (df, is.factor) ) )
  for (i in factor_cols)
    df[,i] <- as.character (df[,i])
  #
  # Check variables names
  #
  data_var_names <- colnames (input_data)
  rows_ok <- rep (T, nrow(df))
  for ( row_idx in 1:nrow(df) )
    {
    for (col_idx in cols_of_variables)
      {
      one_var_name <- df[row_idx, col_idx]
      if (! (one_var_name %in% data_var_names) )
        {
        miic_warning (df_name, "The variable ", one_var_name,
          " is not present in the input data. The row ", row_idx, " will be ignored.")
        rows_ok[[row_idx]] <- F
        }
      }
    if ( ! rows_ok[[row_idx]] )
      next
    #
    # Self loops are not authorized except in temporal modes
    #
    if ( ( ! (mode %in% MIIC_TEMPORAL_MODES) )
       && (df[row_idx, 1] == df[row_idx, 2]) )
      {
      miic_warning (df_name, "the variables must be different for each row (found ",
        df[row_idx, 1], " two times at row ", row_idx, "). This row will be ignored.")
      rows_ok[[row_idx]] <- F
      }
    }
  df <- df [rows_ok, , drop=F]
  if (nrow(df) == 0)
    {
    miic_warning (df_name, "The provided dataframe is empty.")
    return (NULL)
    }
  #
  # In temporal mode, check that the lag column(s)
  #
  if (mode %in% MIIC_TEMPORAL_MODES)
    {
    if (mode == "TS")
      {
      # Temporal stationnary, we expect a 3col with integer >= 0
      #
      test_wrong_lag <- function(x)
        {
        if ( is.null (x) )                                     # NULL: KO
          return (TRUE)
        if ( is.na (x) )                                       # NA: OK for now
          return (FALSE)
        else if ( is.na ( suppressWarnings (as.numeric(x)) ) ) # Not num: KO
          return (TRUE)
        else if ( round(as.numeric(x),0) != as.numeric(x) )    # Not int: KO
          return (TRUE)
        else if ( (as.numeric(x) < 0) )                        # <0: KO
          return (TRUE)
        else
          return (FALSE)                                       # OK
        }
      wrong_lags <- unlist (lapply (df[,3], FUN=test_wrong_lag) )
      }
    else
      {
      # Non stationary, we expect 2nd and 4th columns with integers
      #
      test_wrong_lag <-  function(x)
        {
        if ( is.null (x) )                                     # NULL: KO
         return (TRUE)
        if ( is.na (x) )                                       # NA: KO in non stat
         return (TRUE)
        else if ( is.na ( suppressWarnings (as.numeric(x)) ) ) # Not num: KO
         return (TRUE)
        else if ( round(as.numeric(x),0) != as.numeric(x) )    # Not int: KO
         return (TRUE)
        # else if ( (as.numeric(x) < 0) )                        # <0 is OK
        # return (TRUE)
        else
         return (FALSE)                                       # OK
        }
      wrong_lags2 <- unlist (lapply (df[,2], FUN=test_wrong_lag) )
      wrong_lags4 <- unlist (lapply (df[,4], FUN=test_wrong_lag) )
      wrong_lags <- wrong_lags2 | wrong_lags4
      }
    if ( any (wrong_lags) )
      {
      msg_str <- list_to_str (which(wrong_lags), n_max=10)
      if (sum (wrong_lags) == 1)
        miic_warning (df_name, "lag is incorrect at row ", msg_str,
          ", this line will be ignored.")
      else
        miic_warning (df_name, "lag is incorrect for multiple rows (", msg_str,
          "), these rows will be ignored.")
      df <- df [!wrong_lags, , drop=F]
      }
    if (nrow(df) == 0)
      {
      miic_warning (df_name, "The provided dataframe is empty.")
      return (NULL)
      }
    #
    # Check that contextual lag are NA and no self loop on a contextual var
    # in temporal stationary mode
    #
    wrongs_ctx <- rep ( F, nrow(df) )
    wrongs_lagged <- rep ( F, nrow(df) )
    wrongs_selfs_ctx <- rep ( F, nrow(df) )

    if (mode == "TS")
      {
      contextuals <- unlist ( apply ( df, MARGIN=1, FUN=function (x) {
        # The state_order should have been checked before.
        # So, it should never been NULL or empty. It should also contain an
        # is_contextual column filled with a valid value for each variable
        orig_idx <- which (state_order$var_names == x[[1]])
        dest_idx <- which (state_order$var_names == x[[2]])
        return (  (state_order[orig_idx, "is_contextual"] == 1)
               || (state_order[dest_idx, "is_contextual"] == 1) ) } ) )
      wrongs_ctx <- ( contextuals & ( ! is.na (df[,3]) ) )
      if ( any (wrongs_ctx) )
        {
        if (sum (wrongs_ctx) == 1)
          miic_warning (df_name, "lags for contextual variables must be NA.",
            " The line ", df[wrongs_ctx, 1], " - ", df[wrongs_ctx, 2], " lag ",
            df[wrongs_ctx, 3], " will be ignored.")
        else
          miic_warning (df_name, "lags for contextual variables must be NAs. ",
            sum (wrongs_ctx), " wrong lines will be ignored.")
        }
      #
      # Check that lag >= 0 if not contextual
      #
      wrongs_lagged <- ( (!contextuals) & is.na (df[,3]) )
      if ( any (wrongs_lagged) )
        {
        if (sum (wrongs_lagged) == 1)
          miic_warning (df_name, "lag for non contextual variables must be >= 0.",
            " The line ", df[wrongs_lagged, 1], " - ", df[wrongs_lagged, 2], " lag ",
            df[wrongs_lagged, 3], " will be ignored.")
        else
          miic_warning (df_name, "lags for non contextual variables must be >= 0. ",
            sum (wrongs_lagged), " wrong lines will be ignored.")
        }
      #
      # Check that no self loop with a contextual var
      #
      wrongs_selfs_ctx <- ( contextuals & (df[,1] == df[,2]) )
      wrongs_selfs_ctx[wrongs_ctx | wrongs_lagged] <- F
      if ( any (wrongs_selfs_ctx) )
        {
        if (sum (wrongs_selfs_ctx) == 1)
          miic_warning (df_name, "no self loop possible for a contextual variable.",
            " The line ", df[wrongs_selfs_ctx, 1], " - ", df[wrongs_selfs_ctx, 2], " lag ",
            df[wrongs_selfs_ctx, 3], " will be ignored.")
        else
          miic_warning (df_name, "no self loops possible for contextual variables. ",
            sum (wrongs_selfs_ctx), " wrong lines will be ignored.")
        }
      }
    #
    # The self loops need a lag > 0
    #
    if (mode == "TS")
      wrongs_selfs <- ( (!contextuals) & (df[,1] == df[,2]) & (df[,3] == 0) )
    else
      wrongs_selfs <- ( (df[,1] == df[,3]) & (df[,2] == df[,4]) )
    wrongs_selfs [wrongs_ctx | wrongs_lagged | wrongs_selfs_ctx] <- F
    if ( any (wrongs_selfs) )
      {
      if (sum (wrongs_selfs) == 1)
        miic_warning (df_name, "lag for self loops must be > 0.",
          " The line ", df[wrongs_selfs, 1], " - ", df[wrongs_selfs, 2], " lag ",
          df[wrongs_selfs, 3], " will be ignored.")
      else
        miic_warning (df_name, "lags for self loops must be > 0. ",
          sum (wrongs_selfs), " wrong lines will be ignored.")
      }

    df <- df [ (!wrongs_ctx) & (!wrongs_lagged) & (!wrongs_selfs_ctx) & (!wrongs_selfs), , drop=F]
    if (nrow(df) == 0)
      {
      miic_warning (df_name, "The provided dataframe is empty.")
      return (NULL)
      }
    #
    # Ensure the lags are integer types
    #
    if (mode == "TS")
      df[,3] = as.integer (df[,3])
    else # TNS
      {
      df[,2] = as.integer (df[,2])
      df[,4] = as.integer (df[,4])
      }
    }
  #
  # Remove duplicate row
  #
  n_rows_sav = nrow(df)
  #
  # Equal rows
  #
  df = unique (df)
  rownames(df) = NULL
  #
  # We remove equal rows but with variable names swapped
  # as edges in black box are not oriented and, for true edges,
  # the post-processing will not be able to process opposite edges
  #
  rows_kept = rep (T, nrow(df))
  for (i in 1:nrow(df))
    {
    if ( ! rows_kept[[i]] )
      next
    if (mode == "TS")
      {
      # In temporal stationary, lag != 0 with variable swapped are not duplicate
      #
      if ( (!is.na(df[i,3])) && (df[i,3] != 0) )
        next
      dup_inverse = ( (df[,1] == df[i,2])
                    & (df[,2] == df[i,1])
                    & (rownames(df) != i)
                    & (is.na(df[,3]) | (df[,3] == 0)) )
      }
    else if (mode == "TNS")
      dup_inverse = ( (df[,1] == df[i,3])
                    & (df[,2] == df[i,4])
                    & (df[,3] == df[i,1])
                    & (df[,4] == df[i,2])
                    & (rownames(df) != i) )
    else # non temporal
      dup_inverse = ( (df[,1] == df[i,2])
                    & (df[,2] == df[i,1])
                    & (rownames(df) != i) )
    rows_kept = rows_kept & (!dup_inverse)
    }
  df <- df[rows_kept, , drop=F]
  if ( n_rows_sav != nrow(df) )
    {
    if (df_name == "true edges")
      {
      miic_warning (df_name, "the implementation of truth edges",
        " is not designed to handle opposite edges.",
        " Only one direction will be considered for the opposite edge(s).")
      }
    else
      {
      if (n_rows_sav - nrow(df) == 1)
        miic_warning (df_name, "1 row is duplicated. Only one instance",
          " of the row will be used.")
      else
        miic_warning (df_name, n_rows_sav - nrow(df), " rows are duplicated.",
          " Only one instance of these rows will be used.")
      }
    }

  if (nrow(df) == 0)
    {
    miic_warning (df_name, "The provided data frame is empty.")
    return (NULL)
    }
  return (df)
  }

#-------------------------------------------------------------------------------
# check_param_string
#-------------------------------------------------------------------------------
# Utility function to check a string parameter
# Parameters:
# - value: the parameter to check
# - name: the name of the parameter
# - list: the possible values, 1st value is the default
# Return: the checked parameter, eventually reset to its default value
#-------------------------------------------------------------------------------
check_param_string <- function (value, name, possibles)
  {
  if (  is.null(value)
     || (length (value) != 1)
     || is.na(value)
     || (!is.character(value))
     || (!(value %in% possibles)) )
    {
    msg_str = paste0 (paste0 ("'", possibles, "'"), collapse=", ")
    if ( is.null (value) )
      val_str = "NULL"
    else
      val_str = list_to_str (value)
    miic_warning ("parameters", "supplied value '", val_str,
      "' for the ", name, " parameter is invalid. Possible values are: ",
      msg_str, ". The default value ('", possibles[[1]], "') will be used.")
    value <- possibles[[1]]
    }
  return (value)
  }

#-------------------------------------------------------------------------------
# check_param_logical
#-------------------------------------------------------------------------------
# Utility function to check a logical parameter
# Parameters:
# - value: the parameter to check
# - name: the name of the parameter
# - default: the default value
# Return: the checked parameter, eventually reset to its default value
#-------------------------------------------------------------------------------
check_param_logical <- function (value, name, default)
  {
  if ( is.null (value)
    || (length (value) != 1)
    || is.na (value)
    || (!is.logical(value)) )
    {
    if ( is.null (value) )
      val_str = "NULL"
    else
      val_str = list_to_str (value)
    miic_warning ("parameters", "supplied value ", val_str,
      " for the ", name, " parameter is invalid. It must be TRUE/FALSE.",
      " The default value (", default, ") will be used.")
    value <- default
    }
  return (value)
  }

#-------------------------------------------------------------------------------
# test_param_wrong_int
#-------------------------------------------------------------------------------
# Utility function to detect a wrong int
# Parameters:
# - value: the parameter to check
# - min: the min value, NA if none
# - max: the max values, NA if none
# Return: TRUE if the value is not an int or not in the range, FALSE otherwise
#-------------------------------------------------------------------------------
test_param_wrong_int <- function (value, min=NA, max=NA)
  {
  return (  is.null (value)
         || (length (value) != 1)
         || is.na (value)
         || (!is.numeric(value))
         || (round(value,0) != value)
         || ((!is.na (min)) && (value < min))
         || ((!is.na (max)) && (value > max)) )
  }

#-------------------------------------------------------------------------------
# check_param_int
#-------------------------------------------------------------------------------
# Utility function to check an int parameter
# Parameters:
# - value: the parameter to check
# - name: the name of the parameter
# - min_max: a tuple with min and max values. NA if no min and/or no max
# - default: the default value
# Return: the checked parameter, eventually reset to its default value
#-------------------------------------------------------------------------------
check_param_int <- function (value, name, default, min=NA, max=NA)
  {
  if ( test_param_wrong_int (value, min, max) )
    {
    msg_str = " It must be an integer."
    if ( (!is.na(min)) && (!is.na(max)) )
      msg_str = paste0 (" It must be an integer in the range [",
                        min, ", ", max, "].")
    else if ( ! is.na (min) )
      msg_str = paste0 (" It must be an integer >= ", min, ".")
    else if ( ! is.na (max) )
      msg_str = paste0 (" It must be an integer <= ", max, ".")
    if ( is.null (value) )
      val_str = "NULL"
    else
      val_str = list_to_str (value)
    miic_warning ("parameters", "supplied value ", val_str,
      " for the ", name, " parameter is invalid." , msg_str,
      " The default value (", default, ") will be used.")
    value = default
    }
  return (value)
  }

#-------------------------------------------------------------------------------
# test_param_wrong_float
#-------------------------------------------------------------------------------
# Utility function to detect a wrong float
# Parameters:
# - value: the parameter to check
# - min: the min value, NA if none
# - max: the max values, NA if none
# - strict_min: if TRUE, value must be < min. if FALSE, value must be <= min
# TODO Revoir strict_min pas clair
# Return: TRUE if the value is not a float or not in the range, FALSE otherwise
#-------------------------------------------------------------------------------
test_param_wrong_float <- function (value, min=NA, strict_min=T, max=NA)
  {
  return (  is.null (value)
         || (length (value) != 1)
         || is.na (value)
         || (!is.numeric(value))
         || ((!is.na (min)) && strict_min && (value < min))
         || ((!is.na (min)) && (!strict_min) && (value <= min))
         || ((!is.na (max)) && (value > max)) )
  }

#-------------------------------------------------------------------------------
# check_parameters
#-------------------------------------------------------------------------------
# Check all input parameters that are not data frames and not specific to a mode
# Parameters:
# - input_data: a data frame with input data
# - all possible parameters not specific to a mode of MIIC method
# Return: a list with all the parameters, eventually modified or initialized
#-------------------------------------------------------------------------------
# REVIEW strict_min
check_parameters <- function (input_data, n_threads, cplx,
  orientation, ort_proba_ratio, ort_consensus_ratio, propagation, latent,
  n_eff, n_shuffles, conf_threshold, sample_weights, test_mar,
  consistent, max_iteration, consensus_threshold,
  negative_info, mode, verbose)
  {
  list_ret <- list ("mode" = mode)
  list_ret$n_threads <- check_param_int (n_threads, "n_threads", 1, min=1, max=NA)
  list_ret$cplx <- check_param_string (cplx, "complexity", c("nml", "bic"))
  list_ret$orientation <- check_param_logical (orientation, "orientation", TRUE)

  if ( test_param_wrong_float (ort_proba_ratio, min=0, strict_min=F, max=1) )
    {
    miic_warning ("parameters", "supplied value ", list_to_str (ort_proba_ratio),
      " for the orientation probability ratio parameter is invalid.",
      " It must be a floating point between 0 and 1.",
      # TODO 2.0.3 >=0 vs lmiic > 0 " It must be a floating point in the range ]0,1].",
      " The default value (1) will be used.")
    ort_proba_ratio = 1
    }
  list_ret$ort_proba_ratio <- ort_proba_ratio

  if ( is.null (ort_consensus_ratio) )
    ort_consensus_ratio <-list_ret$ort_proba_ratio
  else if ( test_param_wrong_float (ort_consensus_ratio, min=0, strict_min=F, max=1) )
    {
    miic_warning ("parameters", "supplied value ", list_to_str (ort_consensus_ratio),
      " for the orientation concensus ratio parameter is invalid.",
      " It must be a floating point between 0 and 1.",
      # TODO check if >0 " It must be a floating point in the range ]0,1].",
      " The default value (same as orientation probabilty ratio: ",
      ort_proba_ratio, ") will be used.")
    ort_consensus_ratio <- list_ret$ort_proba_ratio
    }
  list_ret$ort_consensus_ratio <- ort_consensus_ratio

  list_ret$propagation <- check_param_logical (propagation, "propagation", FALSE)
  list_ret$latent <- check_param_string (latent, "latent", MIIC_VALID_LATENT)

  if (  test_param_wrong_int (n_eff, min=-1, max=nrow(input_data) )
     || (n_eff == 0)  )
    {
    miic_warning ("parameters", "supplied value ", list_to_str (n_eff),
      " for the number of effective samples is invalid.",
      " The number of effective samples must be an integer that can be -1",
      " for an automatic assignment or a positive number less or",
      " equal to the number of samples. The default value (-1) will be used.")
    n_eff <- -1
    }
  list_ret$n_eff <- n_eff

  n_shuffles <- check_param_int (n_shuffles, "number of shufflings", 0, min=0, max=NA)
  if (n_shuffles == 0)
    {
    if (  (length(conf_threshold) > 1)
       || (  (!is.null (conf_threshold))
          && (!is.na (conf_threshold))
          && (conf_threshold != 0) ) )
      miic_warning ("parameters", "supplied value ", list_to_str (conf_threshold),
        " for the confidence threshold parameter will be ignored",
        " as the number of shufflings is set to 0.",
        " To activate the confidencence cut, both the number of shufflings",
        " and the confidence threshold must be > 0 (i.e.: n_shuffles = 100",
        " and conf_threshold = 0.01).")
    conf_threshold <- 0
    }
  else
    {
    if ( test_param_wrong_float (conf_threshold, min=0, max=NA) )
      {
      miic_warning ("parameters", "supplied value ", list_to_str (conf_threshold),
        " for the confidence threshold parameter is invalid.",
        " When confidence cut is activated (when n_shuffles > 0),",
        " the confidence threshold must be a floating point > 0. The",
        " confidence cut will be desactivated and default values will be used",
        " for the number of shufflings (0) and the confidence threshold (0).")
      n_shuffles <- 0
      conf_threshold <- 0
      }
    else if (conf_threshold == 0)
      {
      miic_warning ("parameters", "the confidence threshold parameter is 0",
        " but it must be > 0 when confidence cut is activated",
        " (when n_shuffles > 0). The confidence cut will be desactivated.",
        " To activate the confidencence cut, both the number of shufflings",
        " and the confidence threshold must be > 0 (i.e.: n_shuffles = 100",
        " and conf_threshold = 0.01).")
      n_shuffles <- 0
      }
    }
  list_ret$n_shuffles <- n_shuffles
  list_ret$conf_threshold <- conf_threshold

  if (  ( ! is.null (sample_weights) )
     && (  (length(sample_weights) != nrow(input_data))
        || (any(is.na (sample_weights)))
        || (any(!is.numeric(sample_weights)))
        || (any(sample_weights < 0))
        || (any(sample_weights > 1)) ) )
    {
    miic_warning ("parameters", "supplied value for the sample_weights parameter",
      " is invalid. It must be a vector of the same size as the number of",
      " samples in the input data and all weights must be floating points",
      " in the [0,1] range. The parameter will be ignored.")
    sample_weights <- NULL
    }
  list_ret$sample_weights <- sample_weights

  list_ret$test_mar <- check_param_logical (test_mar, "missing at random test", TRUE)
  list_ret$consistent <- check_param_string (consistent, "consistent", MIIC_VALID_CONSISTENT)

  if (list_ret$consistent == "no")
    {
    if ( (length (max_iteration) > 1)
       || (  (!is.null (max_iteration))
          && (!is.na (max_iteration))
          && (max_iteration != 100) ) )
      miic_warning ("parameters", "supplied value ", list_to_str(max_iteration),
        " for the maximum iteration parameter will not be used",
        " as consistency is off.")
    max_iteration <- 100
    if (  (length (consensus_threshold) > 1)
       || (  (!is.null (consensus_threshold))
          && (!is.na (consensus_threshold))
          && (consensus_threshold != 0.8) ) )
      miic_warning ("parameters", "Supplied value ", list_to_str(consensus_threshold),
        " for the consensus threshold parameter will not be used",
        " as consistency is off.")
    consensus_threshold <- 0.8
    }
  else # Consistency on
    {
    if ( test_param_wrong_int (max_iteration, min=1, max=NA) )
      {
      miic_warning ("parameters", "supplied value ", list_to_str (max_iteration),
        " for the maximum iteration parameter is invalid.",
        " It must be a stricly positive integer when consistency is activated.",
        " The default value (100) will be used.")
      max_iteration <- 100
      }
    if ( test_param_wrong_float (consensus_threshold, min=0.5, max=1) )
      {
      miic_warning ("parameters", "supplied value ", list_to_str (consensus_threshold),
        " for the consensus threshold parameter is invalid.",
        " It must be a floating point between 0.5 and 1 when consistency is",
        " activated. The default value (0.8) will be used.")
      consensus_threshold <- 0.8
      }
    }
  list_ret$max_iteration <- max_iteration
  list_ret$consensus_threshold <- consensus_threshold

  list_ret$negative_info <- check_param_logical (negative_info,
    "allowing/disallowing negative shifted mutual information", FALSE)
  list_ret$verbose <- check_param_logical (verbose, "verbose", FALSE)

  return (list_ret)
  }

#-------------------------------------------------------------------------------
# check_cross_inputs
#-------------------------------------------------------------------------------
# Extra checks needing several inputs
#
# Params:
# - the list constructed by check_inputs
# Return:
# - the list constructed by check_inputs
#   (unmodified for now as no check modifies the list)
#-------------------------------------------------------------------------------
check_cross_inputs <- function (list_in)
  {
  if ( ! (list_in$params$mode %in% MIIC_TEMPORAL_MODES) )
    {
    # Warning on discrete columns with level = rows
    # (NB: for temporal modes, it will be done after data lagging)
    #
    for (one_col in list_in$state_order[list_in$state_order$var_type == 0,
                                        "var_names"])
      {
      vals_not_na = list_in$input_data [!is.na(list_in$input_data[,one_col]), one_col]
      # NB: No warning on constant vars, is done in check_input_data
      if (  (length (unique (vals_not_na)) > 1)
         && (length (vals_not_na) == length (unique (vals_not_na)) ) )
        miic_warning ("input data", "the discrete variable ", one_col,
                      " has only one occurence for each level.")
      }
    }
  return (list_in)
  }

#-------------------------------------------------------------------------------
# prepare_inputs
#-------------------------------------------------------------------------------
# This function checks all the inputs that will be, if needed, corrected.
# If not supplied, parameters will be set to their default values.
# It will return a list with all the items needed by the c++ reconstruction part
# => In temporal mode, as the list must be usable directly by the C++ part, the
# input_data, state_order, black_box and true_edges items are the lagged ones
#
# Parameters:
# All possible parameters of the miic method. Default values are the same as the
# miic function to allow a call to this function only to check/fix parameters
#
# Return: A list that contains:
# - input_data: the checked input data, lagged if in temporal mode
# - params: the checked parameters
# - state_order: the checked state_order, lagged if in temporal mode
# - black_box: if supplied, the checked black box, lagged if in temporal mode
# - true_edges: if supplied, the checked true edges, lagged if in temporal mode
# - non_lagged: only in temporal mode, nested list containing the inputs checked
#   but not lagged:
#   * input_data: the checked non lagged input data
#   * state_order: the checked non lagged state_order
#   * black_box: if a black box is supplied, the checked non lagged black box
#   * true_edges: if true edges are supplied, the checked non lagged true edges
#-------------------------------------------------------------------------------
prepare_inputs <- function (input_data,
                            state_order = NULL,
                            true_edges = NULL,
                            black_box = NULL,
                            n_threads = 1,
                            cplx = "nml",
                            orientation = TRUE,
                            ori_proba_ratio = 1,
                            ori_consensus_ratio = NULL,
                            propagation = FALSE,
                            latent = "orientation",
                            n_eff = -1,
                            n_shuffles = 0,
                            conf_threshold = 0,
                            sample_weights = NULL,
                            test_mar = TRUE,
                            consistent = "no",
                            max_iteration = 100,
                            consensus_threshold = 0.8,
                            negative_info = FALSE,
                            mode = "S",
                            n_layers = NULL,
                            delta_t = NULL,
                            movavg = NULL,
                            keep_max_data = FALSE,
                            max_nodes = 50,
                            var_interest = NULL,
                            var_interest_condition = NULL,
                            window_position = "start",
                            verbose = FALSE)
  {
  if ( is.null(mode) || ( ! (mode %in% MIIC_VALID_MODES) ) )
    miic_error ("parameters check", "invalid mode ", mode,
      ". Possible modes are S (Standard), TS (Temporal Stationnary),",
      " TNS (Temporal Non Stationnary).")

  if ( ! is.null(verbose) )
    {
    if ( is.logical(verbose) )
      {
      verbose <- as.integer ( check_param_logical (verbose, "verbose", FALSE) ) + 1
      miic_warning ("parameters", " type expected for verbose is now integer.")
      }
    else
      verbose <- check_param_int (verbose, "verbose", 1, min=0, max=2)
    }
  if (verbose >= 1)
    {
    if (mode == "TS")
      miic_msg ("Start MIIC in temporal stationary mode")
    else if (mode == "TNS")
      miic_msg ("Start MIIC in temporal non stationary mode")
    else if (mode == "L")
      miic_msg ("Start MIIC in layered mode")
    else
      miic_msg ("Start MIIC...")
    }
  #
  # Basic checks applicable with few differences between the different modes.
  # Once these functions done, we will have
  # - a basic cleaned dataset
  # - all parameters not specific to a mode checked
  # - a correct state order will all columns expected for the chosen mode
  #   (but the content of columns for non standard modes is not checked)
  # - other optional data frames (black box and true edges) checked against the
  #   variables in the input data.
  #
  list_ret <- list()
  list_ret$input_data <- check_input_data (input_data, mode)
  list_ret$params <- check_parameters (input_data = list_ret$input_data,
                                       n_threads = n_threads,
                                       cplx = cplx,
                                       orientation = orientation,
                                       ori_proba_ratio = ori_proba_ratio,
                                       ori_consensus_ratio = ori_consensus_ratio,
                                       propagation = propagation,
                                       latent = latent,
                                       n_eff = n_eff,
                                       n_shuffles = n_shuffles,
                                       conf_threshold = conf_threshold,
                                       sample_weights = sample_weights,
                                       test_mar = test_mar,
                                       consistent = consistent,
                                       max_iteration = max_iteration,
                                       consensus_threshold = consensus_threshold,
                                       negative_info = negative_info,
                                       mode = mode,
                                       verbose = verbose)
  list_ret$state_order <- check_state_order (
    list_ret$input_data, state_order, list_ret$params$mode)
  list_ret$black_box <- check_other_df (list_ret$input_data,
    list_ret$state_order, black_box, "black box", list_ret$params$mode)
  list_ret$true_edges <- check_other_df (list_ret$input_data,
    list_ret$state_order, true_edges, "true edges", list_ret$params$mode)
  #
  # Extra checks needing several inputs
  #
  list_ret = check_cross_inputs (list_ret)
  #
  # Other steps depending on the mode
  #
  # NB: these functions are always called whatever the mode. If some parameters
  # specific to a mode are set but for the wrong mode, it will raise a warning.
  # i.e.: number of layers, specific to temporal mode supplied in standard mode
  #
  list_ret <- tmiic_prepare_inputs (list_in = list_ret,
                                    n_layers = n_layers,
                                    delta_t = delta_t,
                                    movavg = movavg,
                                    keep_max_data = keep_max_data,
                                    max_nodes = max_nodes,
                                    var_interest = var_interest,
                                    var_interest_condition = var_interest_condition,
                                    window_position = window_position)
  #
  # Convert TNS into layers
  #
  if (list_ret$params$mode == "TNS")
    list_ret = tmiic_tns_to_ml (list_ret)
  return (list_ret)
  }

#-------------------------------------------------------------------------------
# tmiic_tns_to_ml
#-------------------------------------------------------------------------------
tmiic_tns_to_ml = function(list_in)
  {
  return (list_in)
  }
