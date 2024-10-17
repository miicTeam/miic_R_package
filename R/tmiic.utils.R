#*******************************************************************************
# Filename   : tmiic.utils.R                         Creation date: 10 may 2023
#
# Description: Utility functions for temporal MIIC
#
# Author     : Franck SIMON
#*******************************************************************************

#===============================================================================
# FUNCTIONS
#===============================================================================
# tmiic_check_parameters_not_temporal
#-------------------------------------------------------------------------------
# Raise warnings if some temporal parameters are used in non temporal mode
# Parameters: all possible parameters specific to  temporal modes
# Returns: none
#-------------------------------------------------------------------------------
tmiic_check_parameters_not_temporal <- function (n_layers, delta_t, mov_avg,
  keep_max_data, max_nodes, var_interest, var_interest_condition, window_position)
  {
  if ( ! is.null (n_layers) )
    miic_warning ("parameters", "the n_layers parameter ", n_layers,
      " will not be used as mode is not temporal.")
  if ( ! is.null (delta_t) )
    miic_warning ("parameters", "the delta_t parameter ", delta_t,
      " will not be used as mode is not temporal.")
  if ( ! is.null (mov_avg) )
    miic_warning ("parameters", "the mov_avg parameter ", mov_avg,
      " will not be used as mode is not temporal.")
  if (  (!is.null (keep_max_data))
     && (keep_max_data != FALSE) )
    miic_warning ("parameters", "the keep_max_data parameter ", keep_max_data,
      " will not be used as mode is not temporal.")
  if (  (!is.null (max_nodes))
     && (max_nodes != 50) )
    miic_warning ("parameters", "the max_nodes parameter ", max_nodes,
      " will not be used as mode is not temporal.")
  if ( ! is.null (var_interest) )
    miic_warning ("parameters", "the var_interest parameter ", var_interest,
      " will not be used as mode is not temporal.")
  if ( ! is.null (var_interest_condition) )
    miic_warning ("parameters", "the var_interest_condition parameter ",
      var_interest_condition, " will not be used as mode is not temporal.")
  if (  (!is.null (window_position))
     && (window_position != "start") )
    miic_warning ("parameters", "the window_position parameter ",
      window_position, " will not be used as mode is not temporal.")
  }

#-------------------------------------------------------------------------------
# tmiic_check_parameters
#-------------------------------------------------------------------------------
# Checks on parameters for temporal mode
#
# The "simple" temporal parameters, having one value not tuned per variable,
# are added in the list of parameters. For the parameters that can be different
# per variable: n_layers, delta_t and mov_avg, if these variables are supplied
# as parameters (!= NULL) and not already in the state_order, (otherwise the
# parameters are ignored), these parameters are moved into the state_order.
# The move into the state_order is tuned for each variable, i.e.: mov_avg is
# applied on continuous variables but not on discrete, or, in stationary mode,
# the contextual variables have n_layers=1, delta_t=0.
#
# Parameters:
# - list_in: the list of inputs prepared by the tmiic_check_inputs.
#   The list has 2 items: params and non_lagged.
#   non_lagged is a nested list with non lagged inputs : input_data,
#   state_order and eventually black box, true edges
# - all possible parameters specific to  temporal modes
# Returns: the updated list of inputs, 2 items can be modified:
# - non_lagged$state_order: n_layers, delta_t, mov_avg can be added as new
#   columns in the non lagged state_order:
# - params: the list of parameters will include checked temporal parameters
#-------------------------------------------------------------------------------
tmiic_check_parameters <- function (list_in, n_layers, delta_t, mov_avg,
  keep_max_data, max_nodes, var_interest, var_interest_condition, window_position)
  {
  # Start with parameters not moved to the state order
  #
  if (list_in$params$mode == "TS")
    list_in$params$keep_max_data <- check_param_logical (
      keep_max_data, "keep_max_data", FALSE)
  else if (  (!is.null (keep_max_data))
          && (!is.na (keep_max_data))
          && (keep_max_data != FALSE) )
    miic_warning ("parameters", "the keep_max_data parameter ", keep_max_data,
      " will not be used as mode is temporal non stationnary.")
  #
  # We check only max_nodes > numbers of variables
  # as the number of layers can be unknown at this point
  #
  n_vars <- ncol (list_in$non_lagged$input_data) - 1 # -1 for time steps column
  list_in$params$max_nodes <- check_param_int (
    max_nodes, "maximum number of lagged nodes",
    default=max (50, n_vars + 1), min=n_vars + 1)
  #
  # Check specific parameters of stationary or non stationary mode
  #
  if (list_in$params$mode == "TS")
    {
    if (!is.null (var_interest))
      miic_warning ("parameters", "the variable of interest parameter ", var_interest,
        " will not be used as mode is temporal stationnary.")
    if (!is.null (var_interest_condition))
      miic_warning ("parameters", "the variable of interest condition ", var_interest_condition,
        " will not be used as mode is temporal stationnary.")
    if (  (!is.null (window_position))
       && (!is.na (window_position))
       && (window_position != "start") )
      miic_warning ("parameters", "the window position parameter ",
        window_position, " will not be used as mode is temporal stationnary.")
    }
  else # Not stationary
    {
    if (!is.null (var_interest))
      {
      if (  (length (var_interest) != 1)
         || is.na(var_interest)
         || (!is.character(var_interest))
         || (!(var_interest %in% colnames(list_in$non_lagged$input_data)[-1])) )
        miic_warning ("parameters", "the variable of interest parameter ",
          list_to_str (var_interest), " is invalid and will be ignored.")
      else
        list_in$params$var_interest <- var_interest
      }

    if (is.null (list_in$params$var_interest))
      {
      if (!is.null (var_interest_condition))
        miic_warning ("parameters", "the variable of interest condition ",
          var_interest_condition,
          " will not be used as variable of interest is not defined.")
      }
    else
      {
      list_in$params$var_interest_condition <- tmiic_non_stat_check_var_interest_condition (
        list_in$non_lagged$input_data, list_in$params$var_interest, var_interest_condition)
      if ( is.null (list_in$params$var_interest_condition) )
        list_in$params$var_interest <- NULL
      }

    list_in$params$window_position <- tmiic_non_stat_check_window_position (
      list_in$params$var_interest, window_position)
    }
  #
  # Check number of layers parameter
  #
  if ( ! is.null (n_layers) )
    {
    if ( "n_layers" %in% colnames(list_in$non_lagged$state_order) )
      miic_warning ("parameters", "supplied value ", list_to_str (n_layers),
        " for n_layers parameter will be ignored as the number of layers",
        " is also provided in the state_order.")
    else
      {
      # The n_layers parameter is used, so we need to check it. We have 2 cases:
      # - in non stationary mode around a variable of interest,
      #   n_layers can be of the form n,m (n layers before the position of the
      #   variable of interest and m layers after)
      # - otherwise (in stationary or non stationary + window_position 'start'
      #   or 'end"), it must be an int >= 2
      #
      if (  (list_in$params$mode == "TNS")
         && (list_in$params$window_position == "around") # => var interest used
         && (length (n_layers) == 1)
         && (!is.na(n_layers))
         && is.character(n_layers)
         && (grepl (",", n_layers, fixed=T)) )
        {
        # n_layers can be and is of the form n,m => two int check
        # the minimum acceptable is 2 int values >= 0 and one of the two >= 1
        #
        two_values = strsplit( n_layers, split=",", fixed)[[1]]
        two_values = suppressWarnings( as.numeric(two_values) )
        if (   ( length(two_values) != 2 )
            # || test_param_wrong_int (two_values[[1]], min=0, max=NA)
            # || test_param_wrong_int (two_values[[2]], min=0, max=NA)
            || (  test_param_wrong_int (two_values[[1]], min=1, max=NA)
               && test_param_wrong_int (two_values[[2]], min=1, max=NA) ) )
          #
          # Not two values or no values >= 1, not good at all, ignore the param
          #
          miic_warning ("parameters", "supplied value ", list_to_str (n_layers),
            " for the number of layers is invalid.",
            " If not NULL, it must be an integer >= 2",
            " or a couple of integers 'n,m' such as n >= 1 and m >= 1.",
            " The number of layers will be determined from the data.")
        else
          {
          # Even if not completely ok, we should been able to use the n_layers
          #
          if (test_param_wrong_int (two_values[[1]], min=1, max=NA) )
            {
            miic_warning ("parameters", "first value of ", list_to_str (n_layers),
              " for the number of layers is invalid.",
              " It should be an integer >= 1. This value will be ignored and",
              " the window position will be 'start'.")
            list_in$params$window_position = "start"
            list_in$non_lagged$state_order$n_layers <- as.integer (two_values[[2]]) + 1
            }
          else if (test_param_wrong_int (two_values[[2]], min=1, max=NA) )
            {
            miic_warning ("parameters", "second value of ", list_to_str (n_layers),
              " for the number of layers is invalid.",
              " It should be an integer >= 1. This value will be ignored and",
              " the window position will be 'end'.")
            list_in$params$window_position = "end"
            list_in$non_lagged$state_order$n_layers <- as.integer (two_values[[1]]) + 1
            }
          else # two values are ok
            list_in$non_lagged$state_order$n_layers <- paste0 (two_values, collapse=",")
          }
        }
      else
        {
        # n_layers is not or should not be of the form n,m => check one int >= 2
        #
        if ( (list_in$params$mode == "TNS")
           && (length (n_layers) == 1)
           && (!is.na(n_layers))
           && is.character(n_layers)
           && (grepl (",", n_layers, fixed=T)) )
          miic_warning ("parameters", "supplied value ", list_to_str (n_layers),
            " for the number of layers can not be of the form 'n,m' for",
            " the window position '", list_in$params$window_position,
            "'. If not NULL, it must be an integer >= 2.",
            " The number of layers will be determined from the data.")
        else if ( test_param_wrong_int (n_layers, min=2, max=NA) )
          miic_warning ("parameters", "supplied value ", list_to_str (n_layers),
            " for the number of layers is invalid.",
            " If not NULL, it must be an integer >= 2.",
            " The number of layers will be determined from the data.")
        else # valid n_layers
          {
          list_in$non_lagged$state_order$n_layers <- as.integer (n_layers)
          if (list_in$params$mode == "TS")
            list_in$non_lagged$state_order$n_layers[list_in$non_lagged$state_order$is_contextual == 1] <- as.integer (1)
          }
        }
      }
    }
  #
  # Check delta_t
  #
  if ( ! is.null (delta_t) )
    {
    if ( "delta_t" %in% colnames(list_in$non_lagged$state_order) )
      miic_warning ("parameters", "supplied value ", list_to_str (delta_t),
        " for the delta_t parameter will be ignored as the delta t",
        " is also provided in the state_order.")
    else
      {
      if ( test_param_wrong_int (delta_t, min=1, max=NA) )
        miic_warning ("parameters", "supplied value ", list_to_str (delta_t),
          " for the delta t parameter is invalid.",
          " If not NULL, it must be an integer >= 1.",
          " The delta t will be determined from the data.")
      else # valid delta_t
        {
        list_in$non_lagged$state_order$delta_t <- delta_t
        if (list_in$params$mode == "TS")
          list_in$non_lagged$state_order$delta_t[list_in$non_lagged$state_order$is_contextual == 1] <- 0
        }
      }
    }
  #
  # Check mov_avg
  #
  if ( ! is.null (mov_avg) )
    {
    if ( "mov_avg" %in% colnames(list_in$non_lagged$state_order) )
      miic_warning ("parameters", "supplied value ", list_to_str (mov_avg),
        " for the moving average parameter will be ignored as the moving average",
        " is also provided in the state_order.")
    else
      {
      if (   test_param_wrong_int (mov_avg, min=0, max=NA)
         || (mov_avg == 1) )
        miic_warning ("parameters", "supplied value ", list_to_str (mov_avg),
          " for the moving average parameter is invalid.",
          " If not NULL or 0, it must be an integer >= 2.",
          " The moving average parameter will be ignored.")
      else
        {
        # Valid mov_avg
        #
        list_in$non_lagged$state_order$mov_avg <- mov_avg
        #
        # No mov_avg on discrete
        #
        list_in$non_lagged$state_order$mov_avg[list_in$non_lagged$state_order$var_type == 0] <- 0
        #
        # No mov_avg on contextual vars (constant in stationary
        # and expected to be not 'averageable' in non stationary,
        # i.e.: addition of treatment, cell division, temperature threshold
        # NB: in non stationary, if the user wants a moving average on a
        # contextual variable, it is possible by specifying the mov_avg in
        # the state_order
        #
        list_in$non_lagged$state_order$mov_avg[list_in$non_lagged$state_order$is_contextual == 1] <- 0
        }
      }
    }

  return (list_in)
  }

# TODO review
#-------------------------------------------------------------------------------
# tmiic_check_state_order
#-------------------------------------------------------------------------------
# This function performs the first part checks of the state order columns
# specific to temporal mode: n_layers, delta_t and mov_avg.
#
# In most cases, these columns will not be present at this stage,
# as these information will be likely provided as parameters
# (cf tmiic_check_parameters to see how the n_layers, delta_t and mov_avg
# parameters are moved into the state_order).
#
# Parameters:
# - state_order: a data frame, the state order returned by check_state_order
# - mode : the temporal mode ("TS" or "TNS")
# Return:
# - state_order: the state_order, with temporal parameters eventually modified
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# tmiic_check_state_order
#-------------------------------------------------------------------------------
# Check and prepare the state order for temporal modes.
#
# This function is designed to be called once the check_parameters_temporal
# function has moved (if needed) the n_layers, delta_t and mov_avg parameters
# into the state_order.
# This function check n_layers, delta_t and mov_avg in term of possible types
# and values. Then, it will try to fill possible missing values and will
# check/fix the values against the var_type and is_contextual information.
#
# Params:
# - state_order: a data frame, the state order returned by tmiic_check_parameters
# - params: the list of parameters (needed to know the mode)
# Returns: a list with 2 items:
# - state_order: the state_order, with temporal parameters eventually modified
# - params: the list of parameters eventually modified
#-------------------------------------------------------------------------------
tmiic_check_state_order <- function (list_in)
  {
  n_vars = nrow (list_in$non_lagged$state_order)
  #
  # In stationary mode, contextual are not lagged
  # => we can set n_layers to 1 and delta_t to 0
  # In non stationary, contextual are lagged (unless specified by the user)
  # => we consider them as normal vars
  #
  if (list_in$params$mode == "TS")
    are_contextual <- (list_in$non_lagged$state_order$is_contextual == 1)
  else
    are_contextual <- rep( F, n_vars )
  #
  # Check state order n_layers
  #
  flag_two_int = F
  if (  ( ! ("n_layers" %in% colnames(list_in$non_lagged$state_order)) )
     || all(is.na(list_in$non_lagged$state_order$n_layers)) )
    {
    # If n_layers column is not in state_order (or full of NA in the unlikely
    # but possible case that all variables have been added by check_state_order)
    # then add an integer column n_layers full of NA => automatic estimation
    # and, in stationary mode, set 1 for contextual vars
    #
    list_in$non_lagged$state_order$n_layers <- NA
    list_in$non_lagged$state_order$n_layers[are_contextual] <- 1
    }
  else
    {
    # n_layers is in state_order and some values != NA, check values
    # We have two cases: n_layers can be in the form 'n' or  'n,m'
    #
    wrongs <- unlist (lapply (list_in$non_lagged$state_order$n_layers,
                              FUN=function(x) {
      #
      # NA: OK (missing var added by check_state_order)
      #
      if  (is.na (x))
        return (FALSE)
      #
      # In non stationary mode and around a var of interest, n,m is valid
      #
      is_row_with_two = grepl (",", x, fixed=T)
      if (  (list_in$params$mode == "TNS")
         && (list_in$params$window_position == "around") # => var_interest not null
         && is_row_with_two)
        {
        two_values = strsplit( x, split=",", fixed)[[1]]
        if ( length(two_values) != 2 )                # Not "x,y" => KO
          return (TRUE)
        two_values = trimws (two_values)
        two_values = suppressWarnings( as.numeric(two_values) )
        if (any (is.na (two_values)))                 # Not 2 numeric => KO
          return (TRUE)
        if (any (round(two_values, 0) != two_values)) # Not 2 int => KO
          return (TRUE)
        if (any (two_values < 0))                     # Not >= 0 => KO
          return (TRUE)
        return (FALSE)                                # OK, 2 int >= 0
        }
      #
      # Else only an unique int is acceptable
      #
      if ( is.na ( suppressWarnings (as.numeric(x)) ) ) # Not num: KO
        return (TRUE)
      if ( round(as.numeric(x),0) != as.numeric(x) )    # Not int: KO
        return (TRUE)
      if (as.numeric(x) < 1)                            # Not >= 1: KO
        return (TRUE)
      return (FALSE)                                    # OK, one int >= 1
      } ) )

    if ( any (wrongs) )
      {
      msg_str <- list_to_str (list_in$non_lagged$state_order$var_names[wrongs],
                              n_max=10)
      if (sum (wrongs) == 1)
        miic_warning ("state order", "the number of layers ",
                      list_in$non_lagged$state_order$n_layers[wrongs],
                      " is incorrect for the variable ", msg_str,
                      ", this value will be ignored.")
      else
        miic_warning ("state order", "the number of layers are incorrect for",
                      " several variables (", msg_str,
                      "), these values will be ignored.")
      list_in$non_lagged$state_order$n_layers[wrongs] <- NA
      }
    #
    # Test if we have a mix between 'n' and 'n,m'
    # (as we don't known how to manage a mix)
    #
    rows_with_two = grepl (",", list_in$non_lagged$state_order$n_layers, fixed=T)
    nb_two = sum (rows_with_two)
    nb_na = sum ( is.na(list_in$non_lagged$state_order$n_layers) )
    nb_one = n_vars - nb_two - nb_na
    if (nb_two > nb_one)
      {
      flag_two_int = T
      if (nb_one > 0)
        {
        # More "n,m" than n => drop the n
        #
        miic_warning ("state order", "the number of layers",
          " contains a mix of 'n' and 'n,m'. The 'n' values for variables ",
          list_to_str( list_in$non_lagged$state_order$var_names [!rows_with_two],
                       n_max=10),
          " will be ignored.")
        list_in$non_lagged$state_order$n_layers [!rows_with_two] <- NA
        }
      }
    if ( (nb_one >= nb_two) && (nb_two > 0) )
      {
      # More 'n' than 'n,m' => drop the 'n,m'
      #
      miic_warning ("state order", "the number of layers",
        " contains a mix of 'n' and 'n,m'. The 'n,m' values for variables ",
        list_to_str( list_in$non_lagged$state_order$var_names [rows_with_two],
                     n_max=10),
        " will be ignored.")
      list_in$non_lagged$state_order$n_layers [rows_with_two] <- NA
      }
    #
    # If 2 int, transform into a clean string
    #
    na_in_so <- is.na (list_in$non_lagged$state_order$n_layers)
    if (flag_two_int)
      list_in$non_lagged$state_order$n_layers[!na_in_so] = unlist (lapply (
        list_in$non_lagged$state_order$n_layers[!na_in_so],
        FUN = function(x) {
          two_values = unlist (trimws (strsplit( x, split=",", fixed)[[1]]) )
          return( paste0( two_values[[1]], ",", two_values[[2]] ) )
          }) )
    #
    # Invalid values have been replaced by NA, fill NA by a value if possible
    #
    # For contextual in stationary mode, replace NAs with 1
    # (NB: in non stationary, are_contextual was set to FALSE for all variables)
    #
    if ( any (na_in_so & are_contextual) )
      {
      msg_str <- list_to_str (
        list_in$non_lagged$state_order$var_names[na_in_so & are_contextual],
        n_max=10)
      miic_warning ("state order", "the missing number of layers have been",
        " set to 1 for contextual variables (", msg_str, ").")
      list_in$non_lagged$state_order$n_layers[na_in_so & are_contextual] <- 1
      na_in_so <- is.na (list_in$non_lagged$state_order$n_layers)
      }
    #
    # Looks for remaining NAs
    #
    if ( any (na_in_so) )
      {
      # For remaining vars with n_layers equal to NA:
      # - if no other var has a n_layers, go for automatic estimate
      # - if there is an unique n_layers, apply this value to all
      # - if there is multiple n_layers values,
      #   * In stationary, can't decide, stop
      #   * In non stationary, go for automatic estimate
      #
      uniq_vals <- unique ( list_in$non_lagged$state_order$n_layers[
        (!na_in_so) & (!are_contextual) ] )
      if (  ( (list_in$params$mode == "TS") && (length (uniq_vals) == 0) )
         || ( (list_in$params$mode == "TNS") && (length (uniq_vals) != 1) ) )
        {
        msg_str <- list_to_str (list_in$non_lagged$state_order$var_names[na_in_so],
                                n_max=10)
        miic_warning ("state order", "the missing number of layers will be",
          " determined automatically for variables ", msg_str, ".")
        }
      else if (length (uniq_vals) > 1) # => stationary mode
        {
        miic_error ("state order",
          "some number of layers are missing and they can not be completed",
          " automatically as multiple values are already present.")
        }
      else # length (uniq_vals) == 1, affect the unique value to all NA
        {
        msg_str <- list_to_str (list_in$non_lagged$state_order$var_names[na_in_so], n_max=10)
        miic_warning ("state order", "the missing number of layers will be",
          " set to ", uniq_vals, " for variables ", msg_str, ".")
        list_in$non_lagged$state_order$n_layers[na_in_so & (!are_contextual)] <- uniq_vals
        }
      }
    #
    # Extra check in stationary mode, for contextual vars, n_layers must be 1
    #
    if (list_in$params$mode == "TS")
      {
      wrongs <- ( ( ! is.na (list_in$non_lagged$state_order$n_layers) )
                & (list_in$non_lagged$state_order$n_layers != 1)
                & are_contextual )
      if ( any (wrongs) )
        {
        msg_str <- list_to_str (list_in$non_lagged$state_order$var_names[wrongs],
                                n_max=10)
        if (sum (wrongs) == 1)
          miic_warning ("state order", "the contextual variable ", msg_str,
            " has an invalid number of layers. It will be set to 1.")
        else
          miic_warning ("state order", "several contextual variables (", msg_str,
            ") have an invalid number of layers. They will be set to 1.")
        list_in$non_lagged$state_order$n_layers[wrongs] <- 1
        }
      }
    #
    # Warning if multiple values of n_layers (excluding contextual in stationary)
    #
    uniq_vals <- unique ( list_in$non_lagged$state_order$n_layers[
        (!is.na (list_in$non_lagged$state_order$n_layers))
      & (!are_contextual) ] )
    if (length (uniq_vals) > 1)
      {
      msg_str <- list_to_str (uniq_vals)
      if (list_in$params$mode == "TS")
        miic_warning ("state order", "different values (", msg_str,
          ") have been defined for the number of layers.",
          " Such setting should be avoided unless \"specific\" reason",
          " as the result will likely not be accurate.")
      else
        miic_warning ("state order", "different values (", msg_str,
          ") have been defined for the number of layers.")
      }
    #
    # Test about "0,x" and "x,0"
    #
    if (flag_two_int)
      {
      na_in_so <- is.na (list_in$non_lagged$state_order$n_layers)
      first_vals = unlist( lapply( list_in$non_lagged$state_order$n_layers[!na_in_so],
        FUN=function (x) { strsplit(x, split=",", fixed)[[1]][[1]] } ) )
      second_vals = unlist( lapply( list_in$non_lagged$state_order$n_layers[!na_in_so],
        FUN=function (x) { strsplit(x, split=",", fixed)[[1]][[2]] } ) )
      if ( all (first_vals == "0") )
        {
        miic_warning ("state order", "all the first values of 'n,m' layers are 0.",
          " The window position will be changed to 'start'.")
        list_in$params$window_position = "start"
        list_in$non_lagged$state_order$n_layers[!na_in_so] = as.integer(second_vals) + 1
        flag_two_int = F
        }
      if ( all (second_vals == "0") )
        {
        miic_warning ("state order", "all the second values of 'n,m' layers are 0.",
          " The window position will be changed to 'end'.")
        list_in$params$window_position = "end"
        list_in$non_lagged$state_order$n_layers[!na_in_so] = as.integer(first_vals) + 1
        flag_two_int = F
        }
      }
    #
    # Stop if all nb layers == 1
    #
    if (  (!any (is.na (list_in$non_lagged$state_order$n_layers) ) )
       && ( all (list_in$non_lagged$state_order$n_layers == "0,0") # => Non stationary
          || all (list_in$non_lagged$state_order$n_layers <= 1) ) )
      miic_error ("state order", "there must be one variable",
        " at least with a number of layers > 1.")
    }
  #
  # Check state order delta_t
  #
  if (  ( ! ("delta_t" %in% colnames(list_in$non_lagged$state_order)) )
     || all(is.na(list_in$non_lagged$state_order$delta_t)) )
    {
    # If delta_t column is not in the state_order (or full of NA in the unlikely
    # but possible case that all variables have been added by check_state_order)
    # then add an integer column delta_t full of NA => automatic estimation
    # and, in stationary mode, set 0 for contextual vars
    #
    list_in$non_lagged$state_order$delta_t <- NA
    list_in$non_lagged$state_order$delta_t[are_contextual] <- 0
    }
  else
    {
    # delta_t is in the state_order and some values != NA, check values
    #
    # NB: delta_t has been turned into character and initial NA have been
    # turned into "NA" by the check_state_order function
    # => "NA" must raise warning, as initial NA is not a valid delta_t
    # => true NAs don't raise warnings as they are rows added because the
    #    varible name was missing in the state_order and these true NAs
    #    indicates "will be automatically estimated or determined"
    #
    wrongs <- unlist (lapply (list_in$non_lagged$state_order$delta_t,
                              FUN=function(x) {
      if  (is.na (x))                   # NA: OK (missing row added before)
        return (FALSE)
      else if ( is.na ( suppressWarnings (as.numeric(x)) ) ) # Not num: KO
        return (TRUE)
      else if ( round(as.numeric(x),0) != as.numeric(x) )    # Not int: KO
        return (TRUE)
      else if (as.numeric(x) < 0)                            # Not >= 1: KO
        return (TRUE)
      else
        return (FALSE)                                       # OK
      } ) )
    if ( any (wrongs) )
      {
      msg_str <- list_to_str (list_in$non_lagged$state_order$var_names[wrongs],
                              n_max=10)
      if (sum (wrongs) == 1)
        miic_warning ("state order", "the delta t ",
                      list_in$non_lagged$state_order$delta_t[wrongs],
                      " is incorrect for the variable ", msg_str,
                      ", this value will be ignored.")
      else
        miic_warning ("state order", "the delta t are incorrect for",
                      " several variables (", msg_str,
                      "), these values will be ignored.")
      list_in$non_lagged$state_order$delta_t[wrongs] <- NA
      }
    list_in$non_lagged$state_order$delta_t <- as.integer (
      list_in$non_lagged$state_order$delta_t )
    #
    # Replace NA vals if possible:
    #
    na_in_so <- is.na (list_in$non_lagged$state_order$delta_t)
    if ( any (na_in_so & are_contextual) )
      {
      # For contextual, in stationary mode, replace NAs with 1
      #
      msg_str <- list_to_str (
        list_in$non_lagged$state_order$var_names[na_in_so & are_contextual],
        n_max=10)
      miic_warning ("state order", "the missing delta t have been",
        " set to 0 for contextual variables (", msg_str, ").")
      list_in$non_lagged$state_order$delta_t[na_in_so & are_contextual] <- 0
      }
    #
    # Look if still remaining NAs
    #
    na_in_so <- is.na (list_in$non_lagged$state_order$delta_t)
    if ( any (na_in_so) )
      {
      # For remaining vars with delta_t equal to NA:
      # - if no other var has a delta_t, go for automatic estimate
      # - if there is an unique delta_t, apply this value to all
      # - if there is multiple delta_t values, in stationary, stop
      #
      uniq_vals <- unique ( list_in$non_lagged$state_order$delta_t[
                              (!na_in_so) & (!are_contextual) ] )
      if (  ( (list_in$params$mode == "TS") && (length (uniq_vals) == 0) )
         || ( (list_in$params$mode == "TNS") && (length (uniq_vals) != 1) ) )
        {
        msg_str <- list_to_str (list_in$non_lagged$state_order$var_names[na_in_so],
                                n_max=10)
        miic_warning ("state order", "the missing delta t will be",
          " determined automatically for variables ", msg_str, ".")
        }
      else if (length (uniq_vals) > 1) # => stationary mode
        {
        miic_error ("state order", "the state order contains NAs",
          " for the delta t and it can not be completed",
          " automatically as multiple values are already present.")
        }
      else # affect the unique value to all NA
        {
        msg_str <- list_to_str (list_in$non_lagged$state_order$var_names[na_in_so],
                                n_max=10)
        miic_warning ("state order", "the missing delta t will be",
          " set to ", uniq_vals, " for variables ", msg_str, ".")
        list_in$non_lagged$state_order$delta_t [na_in_so & (!are_contextual)] <- uniq_vals
        }
      }
    #
    # Check/fix invalid values:
    # in stationary mode, for contextual vars, delta_t must be 0
    #
    if (list_in$params$mode == "TS")
      {
      wrongs <- ( (!is.na (list_in$non_lagged$state_order$delta_t))
                & (list_in$non_lagged$state_order$delta_t != 0)
                & are_contextual )
      if ( any (wrongs) )
        {
        msg_str <- list_to_str (list_in$non_lagged$state_order$var_names[wrongs],
                                n_max=10)
        if (sum (wrongs) == 1)
          miic_warning ("state order", "the contextual variable ", msg_str,
            " has an invalid delta t. It will be set to 0.")
        else
          miic_warning ("state order", "several contextual variables (", msg_str,
            ") have an invalid delta t. They will be set to 0.")
        list_in$non_lagged$state_order$delta_t[wrongs] <- 0
        }
      }
    #
    # Warning if multiple values of delta_t (excluding contextual in stationary)
    #
    uniq_vals <- unique ( list_in$non_lagged$state_order$delta_t[
      (!is.na (list_in$non_lagged$state_order$delta_t)) & (!are_contextual) ] )
    if (length (uniq_vals) > 1)
      {
      msg_str <- list_to_str (uniq_vals)
      if (list_in$params$mode == "TS")
        miic_warning ("state order", "different values (", msg_str,
          ") have been defined for the delta t.",
          " Such setting should be avoided unless \"specific\" reason",
          " as the result will likely not be accurate.")
      else
        miic_warning ("state order", "different values (", msg_str,
          ") have been defined for the delta t.")
      }
    #
    # Stop if all delta t == 0
    #
    if (  (!any (is.na (list_in$non_lagged$state_order$delta_t)))
       && (all (list_in$non_lagged$state_order$delta_t <= 0)) )
      miic_error ("state order",
                  "there must be one variable at least with a delta t > 0.")
    }
  #
  # Check state order mov_avg
  #
  if (  ( ! ("mov_avg" %in% colnames(list_in$non_lagged$state_order)) )
     || all(is.na(list_in$non_lagged$state_order$mov_avg)) )
    {
    # Add mov_avg column with 0 for all vars
    #
    list_in$non_lagged$state_order$mov_avg <- 0
    }
  else
    {
    # mov_avg supplied, check values
    #
    wrongs <- unlist (lapply (list_in$non_lagged$state_order$mov_avg, FUN=function(x) {
      if  (is.na (x))                   # NA: OK (missing row added before)
        return (FALSE)
      else if ( is.na ( suppressWarnings (as.numeric(x)) ) ) # Not num: KO
        return (TRUE)
      else if ( round(as.numeric(x),0) != as.numeric(x) )    # Not int: KO
        return (TRUE)
      else if ( (as.numeric(x) != 0) && (as.numeric(x) < 2) ) # <2 and !=0: KO
        return (TRUE)
      else
        return (FALSE)                                       # OK
      } ) )
    if ( any (wrongs) )
      {
      msg_str <- list_to_str (list_in$non_lagged$state_order$var_names[wrongs],
                              n_max=10)
      if (sum (wrongs) == 1)
        miic_warning ("state order", "the moving average is incorrect for",
          " the variable ", msg_str, ", this value will be ignored.")
      else
        miic_warning ("state order", "the moving average are incorrect for",
          " several variables (", msg_str, "), these values will be ignored.")
      list_in$non_lagged$state_order$mov_avg[wrongs] <- NA
      }
    list_in$non_lagged$state_order$mov_avg <- as.integer (list_in$non_lagged$state_order$mov_avg)
    #
    # Replace NA vals by 0
    #
    na_in_so <- is.na (list_in$non_lagged$state_order$mov_avg)
    if ( any (na_in_so) )
      {
      msg_str <- list_to_str (list_in$non_lagged$state_order$var_names[na_in_so],
                              n_max=10)
      miic_warning ("state order", "the missing moving average have been set",
        " to 0 for variables ", msg_str, ".")
      list_in$non_lagged$state_order$mov_avg[na_in_so] <- 0
      }
    #
    # Check/fix invalid values: for discrete vars, no moving average
    #
    wrongs <- ( (list_in$non_lagged$state_order$mov_avg != 0)
              & (list_in$non_lagged$state_order$var_type == 0) )
    if ( any (wrongs) )
      {
      msg_str <- list_to_str (list_in$non_lagged$state_order$var_names[wrongs],
                              n_max=10)
      if (sum (wrongs) == 1)
        miic_warning ("state order", "a moving average can not be applied",
          " on the discrete variable ", msg_str, ".")
      else
        miic_warning ("state order", "moving average operations can not",
        " be applied on discrete variables (", msg_str, ").")
      list_in$non_lagged$state_order$mov_avg[wrongs] <- 0
      }
    #
    # In stationary mode, for contextual vars, no moving average
    #
    wrongs <- ( (list_in$non_lagged$state_order$mov_avg != 0) & are_contextual )
    if ( any (wrongs) )
      {
      msg_str <- list_to_str (list_in$non_lagged$state_order$var_names[wrongs],
                              n_max=10)
      if (sum (wrongs) == 1)
        miic_warning ("state order", "a moving average can not be applied",
          " on the contextual variable ", msg_str, ".")
      else
        miic_warning ("state order", "moving average operations can not",
        " be applied on contextual variables (", msg_str, ").")
      list_in$non_lagged$state_order$mov_avg[wrongs] <- 0
      }
    #
    # In non stationary mode, warning if moving average on contextual vars
    # but don't change the mov_avg value (it is the user choice)
    #
    if (list_in$params$mode == "TNS")
      {
      wrongs <- ( (list_in$non_lagged$state_order$mov_avg != 0)
                & list_in$non_lagged$state_order$is_contextual )
      if ( any (wrongs) )
        {
        msg_str <- list_to_str (list_in$non_lagged$state_order$var_names[wrongs],
                                n_max=10)
        if (sum (wrongs) == 1)
          miic_warning ("state order", "a moving average will be applied",
            " on the contextual variable ", msg_str, ".")
        else
          miic_warning ("state order", "moving average operations will",
          " be applied on contextual variables (", msg_str, ").")
        }
      }
    #
    # Warning if multiple values of moving average excluding 0
    #
    uniq_vals <- unique ( list_in$non_lagged$state_order$mov_avg[
        (list_in$non_lagged$state_order$mov_avg != 0) ] )
    if (length (uniq_vals) > 1)
      {
      msg_str <- list_to_str (uniq_vals)
      miic_warning ("state order", "different values (", msg_str,
        ") have been defined for the moving averages.")
      }
    }
  #
  # Enforce the columns type
  #
  if (flag_two_int)
    list_in$non_lagged$state_order$n_layers = as.character (list_in$non_lagged$state_order$n_layers)
  else
    list_in$non_lagged$state_order$n_layers = as.integer(list_in$non_lagged$state_order$n_layers)
  list_in$non_lagged$state_order$delta_t = as.integer (list_in$non_lagged$state_order$delta_t)
  list_in$non_lagged$state_order$mov_avg = as.integer (list_in$non_lagged$state_order$mov_avg)
  #
  # Cross checks
  #
  if (  (!any (is.na (list_in$non_lagged$state_order$n_layers)))
     && (!any (is.na (list_in$non_lagged$state_order$delta_t))) )
    {
    if (flag_two_int)
      {
      n_layers = unlist (lapply( list_in$non_lagged$state_order$n_layers, FUN=function (x) {
        sum( as.integer( strsplit( x, split=",", fixed )[[1]] ) ) + 1
        }) )
      t_max <- max (n_layers * list_in$non_lagged$state_order$delta_t)
      }
    else
      t_max <- max (  list_in$non_lagged$state_order$n_layers
                    * list_in$non_lagged$state_order$delta_t)
    if (t_max < 2)
      miic_error ("state order", "there must be one variable",
                  " at least with 2 layers and a delta t >= 1.")
    }
  return (list_in)
  }

#-------------------------------------------------------------------------------
# tmiic_mov_avg_onecol
#-------------------------------------------------------------------------------
# Utility function to a apply a moving average over a list
# params:
# - x: the list
# - w: the length of the window
# This moving average is centered, so the first (w-1) %/% 2 and the last
# (w-1) - low_shift items will be filled with NA_real_
#-------------------------------------------------------------------------------
tmiic_mov_avg_onecol <- function (x, w)
  {
  low_shift <- (w-1) %/% 2
  high_shift <- (w-1) - low_shift
  ret <- rep(-1, length(x))
  ret[1:low_shift] <- NA_real_
  # print (head (ret))

  start_idx <- low_shift+1
  end_idx <- length(x) - high_shift
  # i <- start_idx + 1
  # i <- start_idx
  for (i in start_idx:end_idx)
    {
    idx_low <- i - low_shift
    idx_high <- i + high_shift
    # TODO VOIR 2.0.3: ret[i] <- mean (x[idx_low:idx_high], na.rm=TRUE)
    ret[i] <- mean (x[idx_low:idx_high], na.action=na.omit)
    }
  # print (head (ret))
  # print (tail (ret))
  ret[(end_idx+1):length(ret)] <- NA_real_
  # print (tail (ret))
  return (ret)
  }

#-------------------------------------------------------------------------------
# tmiic_mov_avg
#-------------------------------------------------------------------------------
# Apply moving averages on data
# - list_traj: a list of data frames, each item representing a trajectory.
#   Each data frame must contain the time step information in the 1st column
#   and the variables in the other columns.
# - mov_avg: the list of moving average to be applied, optional, NULL by defaut.
#   The length of the mov_avg list is the number of columns of the data frames - 1
#   (because the 1st column in data frames is the time step).
#   When the mov_avg item value is >= 2, a moving average using this value as
#   window size is applied on the corresponding column: mov_avg item 1 is
#   applied to data column 2, mov_avg item 2 to data column 3, ...
# - keep_max_data: boolean flag, optional, FALSE by default
#   When FALSE, the rows containing NA introduced by the moving average(s)
#   are deleted, otherwise when TRUE, the rows are kept
# - verbose_level: integer in the range [0,2], 1 by default. The level of
#   verbosity: 0 = no display, 1 = summary display, 2 = maximum display.
# Returns:
# - list_traj: the list trajectories with moving averages applied
#-------------------------------------------------------------------------------
tmiic_mov_avg <- function (list_traj, mov_avg=NULL, keep_max_data=F, verbose_level=0)
  {
  if ( is.null (mov_avg) || all (mov_avg < 2) )
    return (list_traj)
  if (verbose_level >= 1)
    miic_msg ("Applying moving averages...")
  #
  # Apply mov_avg on each trajectory and variable of the dataset
  #
  n_vars <- ncol(list_traj[[1]])-1
  var_names <- colnames (list_traj[[1]])[-1]
  for (i in 1:length(list_traj) )
    for (j in 1:n_vars)
      if (mov_avg[[j]] >= 2)
        {
        # print (paste0 (j, " => mov_avg = ", mov_avg[[j]]))
        list_traj[[i]][,j+1] <- tmiic_mov_avg_onecol (list_traj[[i]][,j+1], mov_avg[[j]])
        if (verbose_level >= 2)
          miic_msg ("- ", var_names[[j]], ": moving average of window size ",
                    mov_avg[[j]], " applied")
        }
  #
  # Remove starting and ending rows where NAs were introduced
  #
  if (!keep_max_data)
    {
    mov_avg_max <- max(mov_avg)
    low_shift <- (mov_avg_max-1) %/% 2
    high_shift <- (mov_avg_max-1) - low_shift
    start_idx <- 1
    if (low_shift > 0)
      start_idx <- start_idx + low_shift
    # i <- 1
    for (i in 1:length(list_traj) )
      {
      end_idx <- nrow(list_traj[[i]]) - high_shift
      list_traj[[i]] <- list_traj[[i]][start_idx:end_idx,]
      }
    }
  return (list_traj)
  }

#-------------------------------------------------------------------------------
# tmiic_check_data_after_lagging
#-------------------------------------------------------------------------------
# Check variables and rows full of NAs after lagging, delete rows full of NAs.
# Raise a warning if discrete variables have only 1 occurrence per level
#
# Params:
# - lagged_data: a data frame, the lagged input data
# - lagged_so: a data frame, the lagged state order
# Returns:
# - a data frame: the lagged input data verified
#-------------------------------------------------------------------------------
tmiic_check_data_after_lagging <- function (lagged_data, lagged_so)
  {
  # Warning on columns full of NAs
  #
  cols_only_na <- colSums (is.na (lagged_data)) == nrow (lagged_data)
  if ( any (cols_only_na) )
    {
    if ( sum (cols_only_na) == 1)
      miic_warning ("lagged data", "the variable ", colnames(lagged_data)[cols_only_na],
                    " contains only NAs after lagging.")
    else
      miic_warning ("lagged data",  sum(cols_only_na), " variables (",
                    list_to_str (colnames(lagged_data)[cols_only_na], n_max=10),
                    ") contains only NAs after lagging.")
    }
  #
  # Warning on discrete columns with level = rows
  #
  for (one_col in lagged_so$var_names[lagged_so$var_type == 0])
    {
    vals_not_na = lagged_data [!is.na(lagged_data[,one_col]), one_col]
    # NB: No warning on constant vars, will be done in check_lagged_state_order
    if (  (length (unique (vals_not_na)) > 1)
       && (length (vals_not_na) == length (unique (vals_not_na)) ) )
      miic_warning ("lagged data", "the discrete variable ", one_col,
                    " has only one occurence for each level.")
    }
  #
  # Warning on rows full of NAs + remove rows
  #
  rows_only_na <- rowSums (is.na (lagged_data)) == ncol (lagged_data)
  if ( any (rows_only_na) )
    {
    if ( sum (rows_only_na) == 1)
      miic_warning ("lagged data", "one row contains only NAs after lagging.",
                    " It will be discarded.")
    else
      miic_warning ("lagged data",  sum(rows_only_na),
        " rows contain only NAs after lagging. They will be discarded.")
    lagged_data = lagged_data[!rows_only_na, , F]
    }
  return (lagged_data)
  }

#-------------------------------------------------------------------------------
# tmiic_check_state_order_after_lagging
#-------------------------------------------------------------------------------
# The var_type in the state_order may need to be re-evaluated after lagging
# as some numerical lagged variables can have less unique values
# and will/can not be any more considered as continuous
#
# Params :
# - lagged_data: a data frame, the lagged input data
# - lagged_so: a data frame, the lagged state order
# Returns:
# - a data frame: the lagged state_order, modified if needed
#-------------------------------------------------------------------------------
tmiic_check_state_order_after_lagging <- function (lagged_data, lagged_so)
  {
  n_unique_vals <- unlist (lapply (lagged_data, function (x) {
    length (unique (x[!is.na(x)] ) ) } ) )
  for (i in 1:nrow (lagged_so))
    {
    if (n_unique_vals[[i]] == 1)
      {
      miic_warning ("lagged data", "the variable ", lagged_so[i, "var_names"],
        " is constant after lagging.")
      lagged_so[i, "var_type"] <- as.integer(0)
      next
      }
    if (lagged_so[i, "var_type"] == 0)
      next
    if (n_unique_vals[[i]] == 2)
      {
      if (  ("var_type_specified" %in% colnames(lagged_so) )
         && (lagged_so[i, "var_type_specified"]) )
        miic_warning ("lagged data", "the variable ", lagged_so[i, "var_names"],
          " was specified as continuous but contains only two values",
          " after lagging. It will be considered as discrete.")
      lagged_so[i, "var_type"] <- as.integer(0)
      next
      }
    if (n_unique_vals[[i]] < MIIC_CONTINUOUS_TRESHOLD)
      {
      if (  ("var_type_specified" %in% colnames(lagged_so) )
         && (lagged_so[i, "var_type_specified"]) )
        miic_warning ("lagged data", "the variable ", lagged_so[i, "var_names"],
          " was specified as continuous but contains only ", n_unique_vals[[i]],
          " values after lagging.")
      else
        lagged_so[i, "var_type"] <- as.integer(0)
      }
    }
  lagged_so$var_type_specified <- NULL
  return (lagged_so)
  }

#-------------------------------------------------------------------------------
# tmiic_check_other_df_after_lagging
#-------------------------------------------------------------------------------
# Check the optional data frame true edges or black box after lagging.
# It complements the previous checks by verifying that lagged edges in these
# data frames exist in the lagged data
#
# Params :
# - var_names: a list, the list of lagged variables names
# - lagged_df: a data frame, the lagged true edges or black box
# - df_name: the data frame name, "true edges" or "black box"
# Returns:
# - a data frame: the lagged data frame, eventually modified
#-------------------------------------------------------------------------------
tmiic_check_other_df_after_lagging <- function (var_names, lagged_df, df_name)
  {
  if ( is.null(lagged_df) )
    return (lagged_df)
  all_varnames_in_df <- unique (c (lagged_df[,1], lagged_df[,2]) )
  vars_absent <- ( ! (all_varnames_in_df %in% var_names) )
  if (any (vars_absent))
    {
    if (sum (vars_absent) == 1)
      miic_warning (df_name, "the variable ", all_varnames_in_df[vars_absent],
                    " is not present in the lagged data.",
                    " Row(s) with this variable will be ignored.")
    else
      miic_warning (df_name, "several variables (",
                    list_to_str (all_varnames_in_df[vars_absent], n_max=10),
                    ") are not present in the lagged data.",
                    " Row(s) with these variables will be ignored.")
    }
  rows_ok <- ( (lagged_df[,1] %in% var_names)
             & (lagged_df[,2] %in% var_names) )
  lagged_df <- lagged_df[rows_ok, , drop=F]
  return (lagged_df)
  }

#-------------------------------------------------------------------------------
# tmiic_prepare_inputs
#-------------------------------------------------------------------------------
# Extra steps to check and prepare the inputs in temporal mode
#
# Params :
# - list_in: the list of inputs, after non temporal checks done
# - all temporal parameters
#
# Return: a list.
# - if called with a non temporal mode, it is the input list unchanged
# - if called with a temporal mode, the list will contain the following items:
#   * input_data: the input data checked and lagged
#   * params: the checked parameters
#   * state_order: the state_order checked and lagged
#   * black_box: if supplied, the black box, checked and lagged
#   * true_edges: if supplied, the true edges, checked and lagged
#   * non_lagged: a nested list containing the inputs checked but not lagged:
#     ° input_data: the checked non lagged input data
#     ° state_order: the checked non lagged state_order
#     ° black_box: if a black box is supplied, the checked non lagged black box
#     ° true_edges: if true edges are supplied, the checked non lagged true edges
#-------------------------------------------------------------------------------
tmiic_prepare_inputs <- function (list_in,
                                  n_layers,
                                  delta_t,
                                  mov_avg,
                                  keep_max_data,
                                  max_nodes,
                                  var_interest,
                                  var_interest_condition,
                                  window_position)
  {
  # First of all, we check temporal parameters when not in temporal mode:
  # raise warnings if some temporal parameters are specified
  #
  if ( ! (list_in$params$mode %in% MIIC_TEMPORAL_MODES) )
    {
    tmiic_check_parameters_not_temporal (
      n_layers = n_layers,
      delta_t = delta_t,
      mov_avg = mov_avg,
      keep_max_data = keep_max_data,
      max_nodes = max_nodes,
      var_interest = var_interest,
      var_interest_condition = var_interest_condition,
      window_position = window_position)
    return (list_in)
    }
  #
  # We are in temporal mode, then we need to return at the end of the function
  # the list of inputs with a structure usable to call the reconstruct function.
  # The expected structure of the list is 'params', (lagged) 'input_data',
  # (lagged) 'state_order' and, if supplied, (lagged) 'black_box' and
  # (lagged) 'true_edges'
  # => The first step is to move the non lagged inputs into a nested list
  # 'non_lagged', the lagged version of the inputs will be added when calling
  # the different functions below
  #
  list_ret <- list ("params" = list_in$params, "non_lagged" = list_in)
  list_ret$non_lagged$params = NULL
  #
  # Check the temporal params, initialize to default value the non supplied
  # NB: if n_layers, delta_t and mov_avg are supplied as parameters, but not
  # also supplied in the state order, they are "moved" in the state order
  #
  list_ret <- tmiic_check_parameters (
    list_in = list_ret,
    n_layers = n_layers,
    delta_t = delta_t,
    mov_avg = mov_avg,
    keep_max_data = keep_max_data,
    max_nodes = max_nodes,
    var_interest = var_interest,
    var_interest_condition = var_interest_condition,
    window_position = window_position)
  #
  # Init / check (/ harmonize if possible) n_layers, delta_t and mov_avg
  #
  list_ret <- tmiic_check_state_order (list_in = list_ret)
  #
  # Prepare trajectories
  #
  list_traj <- tmiic_extract_trajectories (list_ret$non_lagged$input_data)
  verbose_level <- ifelse (list_ret$params$verbose, 2, 1)
  list_traj <- tmiic_mov_avg (list_traj, list_ret$non_lagged$state_order$mov_avg,
    keep_max_data=list_ret$params$keep_max_data, verbose_level=verbose_level)
  #
  # The way we estimate the temporal window and the lagging depend on the
  # stationary or non stationary mode
  #
  if (list_ret$params$mode == "TS") # Stationary
    {
    # Estimate dynamic (if n layers and delta t are not specified by the user)
    #
    list_ret$non_lagged$state_order <- tmiic_stat_estimate_dynamic (list_traj,
      list_ret$non_lagged$state_order, max_nodes=list_ret$params$max_nodes,
      verbose_level=verbose_level)
    #
    # Lag inputs according to n layers and delta t
    #
    list_ret$state_order <- tmiic_stat_lag_state_order (list_ret$non_lagged$state_order)
    list_ret$true_edges <- tmiic_stat_lag_other_df (list_ret$non_lagged$state_order,
                                                    list_ret$non_lagged$true_edges)
    list_ret$black_box <- tmiic_stat_lag_other_df (list_ret$non_lagged$state_order,
                                                   list_ret$non_lagged$black_box)
    list_traj_lagged <- tmiic_stat_lag_input_data (list_traj,
      list_ret$state_order, keep_max_data=list_ret$params$keep_max_data)
    }
  else # "TNS", non stationary
    {
    # Align the trajectories
    #
    list_tmp <- tmiic_non_stat_align_data (list_traj = list_traj,
                                           params = list_ret$params)
    list_traj <- list_tmp$list_traj
    list_ret$params <- list_tmp$params
    pos_interest <- list_tmp$pos_interest
    #
    # Determine (if not specified by the user) the layers and delta t
    #
    list_ret$non_lagged$state_order <- tmiic_non_stat_cover_dynamic (
                            list_traj = list_traj,
                            state_order = list_ret$non_lagged$state_order,
                            params = list_ret$params,
                            pos_interest = pos_interest,
                            verbose = list_ret$params$verbose)
    #
    # The steps could have been skipped, we would be able to create lagged
    # samples without this step. However, as we will return all the inputs
    # before lagging, apply extra step to return the inputs fitting
    # in the best way
    #
    list_tmp <- tmiic_non_stat_tune_inputs (list_traj = list_traj,
                                            state_order = list_ret$non_lagged$state_order,
                                            params = list_ret$params,
                                            pos_interest = pos_interest)
    list_traj <- list_tmp$list_traj
    list_ret$params <- list_tmp$params
    list_ret$non_lagged$state_order <- list_tmp$state_order
    pos_interest <- list_tmp$pos_interest
    #
    # Lag inputs according to n layers and delta t
    #
    list_ret$state_order <- tmiic_non_stat_lag_state_order (
      state_order = list_ret$non_lagged$state_order)
    list_ret$true_edges <- tmiic_non_stat_lag_other_df (
      df = list_ret$non_lagged$true_edges)
    list_ret$black_box <- tmiic_non_stat_lag_other_df (
      df = list_ret$non_lagged$black_box)
    list_traj_lagged <- tmiic_non_stat_lag_input_data (
      list_traj = list_traj,
      lagged_order = list_ret$state_order,
      pos_interest = pos_interest)
    }
  #
  # Even if not used by MIIC, we return the checked inputs before lagging.
  # This can be useful for the user especially for the non stationary mode
  # as the trajectories can have been aligned on a position of interest
  #
  list_ret$non_lagged$input_data <- tmiic_group_trajectories (list_traj)
  #
  # The lagged trajectories will be the input of MIIC reconstruct
  #
  list_ret$input_data <- tmiic_group_trajectories (list_traj_lagged)
  #
  # Checks post lagging:
  #
  # For state_order, check the number of unique values per variable
  # and review discrete/continuous after lagging as some columns may have less
  # number of unique values
  #
  list_ret$state_order <- tmiic_check_state_order_after_lagging (
    list_ret$input_data, list_ret$state_order)
  #
  # Check columns and rows full of NAs and discrete with 1 occurrence per level
  #
  list_ret$input_data <- tmiic_check_data_after_lagging (
    list_ret$input_data, list_ret$state_order)
  #
  # For other df, check that the lagged variables are in the lagged data
  #
  list_ret$true_edges <- tmiic_check_other_df_after_lagging (
    list_ret$state_order$var_names, list_ret$true_edges, "true edges")
  list_ret$black_box <- tmiic_check_other_df_after_lagging (
    list_ret$state_order$var_names, list_ret$black_box, "black box")
  #
  # In stationary, adjust n_eff if delta_t > 1 and no eff supplied by the user
  #
  if ( (list_ret$params$mode == "TS") && (list_ret$params$n_eff == -1) )
    {
    avg_delta_t <- mean (list_ret$state_order$delta_t [list_ret$state_order$is_contextual == 0])
    if (avg_delta_t > 1)
      {
      list_ret$params$n_eff <- round (nrow (list_ret$input_data) / avg_delta_t, 0)
      miic_msg ("Note : the n_eff has been set to ", list_ret$params$n_eff,
                " (nb lagged samples= ", nrow (list_ret$input_data),
                " / delta_t=", round(avg_delta_t, 2), ").")
      }
    }
  #
  # Clean the returned state_order from column used internally
  #
  list_ret$non_lagged$state_order <- list_ret$non_lagged$state_order[,
    colnames(list_ret$non_lagged$state_order) %in% STATE_ORDER_TEMPORAL_VALID_COLUMNS]
  #
  # The non stationary mode has been integrated as part of multi-layers version
  #
  if (list_ret$params$mode == "TNS")
    list_ret <- tmiic_non_stat_to_ml (list_ret)
  print (list_ret$layers)
  print (class (list_ret$layers$layer))
  print (class (list_ret$layers$connected))
  print (class (list_ret$layers$contributors))

  return (list_ret)
  }
