#*****************************************************************************
# Filename   : tmiic.utils.R                   Creation date: 11 may 2023
#
# Description: Utility functions for temporal MIIC
#
# Author     : Franck SIMON
#
# Changes history:
# - 10 may  2023 : initial version
#*****************************************************************************

#-------------------------------------------------------------------------------
# tmiic_check_state_order_part1
#-------------------------------------------------------------------------------
# This function performs the first part checks of the state order columns
# specific to temporal mode: n_layers, delta_t and movavg.
# In most cases, these columns will not be present at this stage,
# as these information will be likely provided as parameters
# (cf tmiic_check_parameters to see how the n_layers, delta_t and movavg
# parameters are moved into the state_order).
# Checks here are basic and cover NULL, integer type and minimal values only.
# NAs are excluded from warnings (NA = row added because var name missing)
#
# Params:
# - input_data: a dataframe with input data
# - state_order: a dataframe, the state order returned by check_state_order
# Returns: a list with 2 items:
# - state_order: the state_order, with temporal parameters eventually modified
#-------------------------------------------------------------------------------
tmiic_check_state_order_part1 <- function (state_order)
  {
  # n_layers check
  #
  if ("n_layers" %in% colnames (state_order) )
    {
    wrongs = unlist (lapply (state_order$n_layers, FUN=function(x) {
      if  (is.na (x))                   # NA: OK (missing row added before)
        return (FALSE)
      else if ( is.na ( suppressWarnings (as.numeric(x)) ) ) # Not num: KO
        return (TRUE)
      else if ( round(as.numeric(x),0) != as.numeric(x) )    # Not int: KO
        return (TRUE)
      else if (as.numeric(x) < 1)                            # Not >= 1: KO
        return (TRUE)
      else
        return (FALSE)                                       # OK
      } ) )
    if ( any (wrongs) )
      {
      msg_str <- list_to_str (state_order$var_names[wrongs], n_max=10)
      if (sum (wrongs) == 1)
        miic_warning ("state order", "the number of layers is incorrect for",
          " the variable ", msg_str, ", this value will be ignored.")
      else
        miic_warning ("state order", "the number of layers are incorrect for",
          " several variables (", msg_str, "), these values will be ignored.")
      state_order$n_layers[wrongs] = NA
      }
    state_order$n_layers = as.integer (state_order$n_layers)
    }
  #
  # delta_t check
  #
  if ("delta_t" %in% colnames (state_order) )
    {
    wrongs = unlist (lapply (state_order$delta_t, FUN=function(x) {
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
      msg_str <- list_to_str (state_order$var_names[wrongs], n_max=10)
      if (sum (wrongs) == 1)
        miic_warning ("state order", "the delta t is incorrect for",
          " the variable ", msg_str, ", this value will be ignored.")
      else
        miic_warning ("state order", "the delta t are incorrect for",
          " several variables (", msg_str, "), these values will be ignored.")
      state_order$delta_t[wrongs] = NA
      }
    state_order$delta_t = as.integer (state_order$delta_t)
    }
  #
  # movavg check
  #
  if ("movavg" %in% colnames (state_order) )
    {
    wrongs = unlist (lapply (state_order$movavg, FUN=function(x) {
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
      msg_str <- list_to_str (state_order$var_names[wrongs], n_max=10)
      if (sum (wrongs) == 1)
        miic_warning ("state order", "the moving average is incorrect for",
          " the variable ", msg_str, ", this value will be ignored.")
      else
        miic_warning ("state order", "the moving average are incorrect for",
          " several variables (", msg_str, "), these values will be ignored.")
      state_order$movavg[wrongs] = NA
      }
    state_order$movavg = as.integer (state_order$movavg)
    }
  return (state_order)
  }

#-------------------------------------------------------------------------------
# tmiic_check_parameters
#-------------------------------------------------------------------------------
# Checks on parameters for temporal mode
#
# As the temporal parameters n_layers, delta_t, movavg need to take different
# values depending on the type (discrete/continuous) or contextual,
# these parameters are moved in the state_order to have a value defined
# for each variable (unless these information are already in the state_order,
# in such case, the parameters are ignored).
# The other temporal parameters, having one value not tuned per variable
# are added in the list of parameters.
#
# Params:
# - state_order: the dataframe  returned by tmiic_check_state_order_part1
# - params: the list of parameters (used only to add temporal parameters)
# - all possible temporal parameters of miic method
# Returns: a list with 2 items:
# - state_order: the state_order, with temporal parameters eventually added
# - params: the list of parameters with temporal parameters added
#-------------------------------------------------------------------------------
tmiic_check_parameters <- function (state_order, params,
  n_layers, delta_t, movavg, keep_max_data, max_nodes)
  {
  # Check number of layers parameter
  #
  if ( ! is.null (n_layers) )
    {
    if ( test_param_wrong_int (n_layers, min=2, max=NA) )
      {
      if ( "n_layers" %in% colnames(state_order) )
        miic_warning ("parameters", "supplied value ", list_to_str (n_layers),
          " for the number of layers is invalid,",
          " if not NULL, it must be an integer >= 2.",
          " This issue has no impact as the number of layers is provided",
          " in the state_order.")
      else
        miic_warning ("parameters", "supplied value ", list_to_str (n_layers),
          " for the number of layers is invalid,",
          " if not NULL, it must be an integer >= 2.",
          " The number of layers will be estimated from the data.")
      }
    else # valid n_layers
      {
      if ( ! ("n_layers" %in% colnames(state_order)) )
        {
        state_order$n_layers = n_layers
        state_order$n_layers[state_order$is_contextual == 1] = 1
        }
      else # n_layers in state_order
        {
        na_in_so = is.na (state_order$n_layers)
        if ( any (na_in_so) )
          {
          miic_warning ("parameters", "the number of layers is supplied both",
            " in the state_order and as parameter. As some values are missing",
            " in the state_order, the parameter will be used to fill these",
            " missing values.")
          state_order$n_layers[na_in_so & (state_order$is_contextual == 0)] = n_layers
          }
        else
          miic_warning ("parameters", "the number of layers is supplied both",
            " in the state_order and as parameter. The parameter will be",
            " ignored.")
        }
      }
    }
  #
  # Check delta_t
  #
  if ( ! is.null (delta_t) )
    {
    if ( test_param_wrong_int (delta_t, min=1, max=NA) )
      {
      if ( "delta_t" %in% colnames(state_order) )
        miic_warning ("parameters", "supplied value ", list_to_str (delta_t),
          " for the delta t parameter is invalid,",
          " if not NULL, it must be an integer >= 1.",
          " This issue has no impact as the delta t is provided",
          " in the state_order.")
      else
        miic_warning ("parameters", "supplied value ", list_to_str (delta_t),
          " for the delta t parameter is invalid,",
          " if not NULL, it must be an integer >= 1.",
          " The delta t will be estimated from the data.")
      }
    else # valid delta_t
      {
      if ( ! ("delta_t" %in% colnames(state_order)) )
        {
        state_order$delta_t = delta_t
        state_order$delta_t[state_order$is_contextual == 1] = 0
        }
      else # delta_t in state_order
        {
        na_in_so = is.na (state_order$delta_t)
        if ( any (na_in_so) )
          {
          miic_warning ("parameters", "the delta t is supplied both",
            " in the state_order and as parameter. As some values are missing",
            " in the state_order, the parameter will be used to fill these",
            " missing values.")
          state_order$delta_t[na_in_so & (state_order$is_contextual == 0)] = delta_t
          }
        else
          miic_warning ("parameters", "the delta t is supplied both",
            " in the state_order and as parameter. The parameter will be",
            " ignored.")
        }
      }
    }
  #
  # Check movavg
  #
  if ( ! is.null (movavg) )
    {
    if (   test_param_wrong_int (movavg, min=0, max=NA)
       || (movavg == 1) )
      {
      if ( "movavg" %in% colnames(state_order) )
        miic_warning ("parameters", "supplied value ", list_to_str (movavg),
          " for the moving average parameter is invalid,",
          " if not NULL or 0, it must be an integer >= 2.",
          " This issue has no impact as the moving average is provided",
          " in the state_order.")
      else
        miic_warning ("parameters", "supplied value ", list_to_str (movavg),
          " for the moving average parameter is invalid,",
          " if not NULL or 0, it must be an integer >= 2.",
          " The moving average parameter will be ignored.")
      }
    else # valid movavg
      {
      if ( ! ("movavg" %in% colnames(state_order)) )
        {
        state_order$movavg = movavg
        # No movavg on discrete or contextual vars
        state_order$movavg[state_order$var_type == 0] = 0
        state_order$movavg[state_order$is_contextual == 1] = 0
        }
      else # movavg in state_order
        {
        na_in_so = is.na (state_order$movavg)
        if ( any (na_in_so) )
          {
          miic_warning ("parameters", "the moving average is supplied both",
            " in the state_order and as parameter. As some values are missing",
            " in the state_order, the parameter will be used to fill these",
            " missing values.")
          state_order$movavg[ na_in_so
                            & (state_order$var_type == 1)
                            & (state_order$is_contextual == 0)] = movavg
          }
        else
          miic_warning ("parameters", "the moving average is supplied both",
            " in the state_order and as parameter. The parameter will be",
            " ignored.")
        }
      }
    }

  params$keep_max_data = check_param_logical (keep_max_data, "keep_max_data", FALSE)
  params$max_nodes = check_param_int (max_nodes, "maximum number of lagged nodes",
                                      default=50, min=nrow(state_order)+1)

  return (list ("params"=params, "state_order"=state_order))
  }

#-------------------------------------------------------------------------------
# tmiic_check_state_order_part2
#-------------------------------------------------------------------------------
# Second part of the check state order for temporal mode.
# This function is designed to be called after the check_parameters_temporal
# function has moved (if needed) the n_layers, delta_t and movavg parameters
# into the state_order.
# This function will try to fill possible missing values and will check/fix
# the temporal settings against the var_type and is_contextual information.
#
# Params :
# - state_order: a dataframe, the state order returned by tmiic_check_parameters
# Returns: a list with 2 items:
# - state_order: the state_order, with temporal parameters eventually modified
#-------------------------------------------------------------------------------
tmiic_check_state_order_part2 <- function (state_order)
  {
  # Check state order n_layers. During the check_state_order function,
  # we already checked NULL, not integer and <1 values, they became NA.
  #
  if ( ! ("n_layers" %in% colnames(state_order)) )
    {
    # Add n_layers column with 1 for is_contextual, NA otherwise
    #
    state_order$n_layers = NA
    state_order$n_layers[ state_order$is_contextual == 1] = 1
    }
  else
    {
    # Replace NA vals if possible:
    #
    na_in_so = is.na (state_order$n_layers)
    are_contextual = (state_order$is_contextual == 1)
    if ( any (na_in_so & are_contextual) )
      {
      # For contextual, replace NAs with 1
      #
      msg_str = list_to_str (state_order$var_names[na_in_so & are_contextual],
                             n_max=10)
      miic_warning ("temporal checks", "the missing number of layers have been",
        " set to 1 for contextual variables (", msg_str, ").")
      state_order$n_layers[na_in_so & are_contextual] = 1
      }
    #
    # Contextual vars done, look if still NAs on not contextual
    #
    na_in_so = is.na (state_order$n_layers)
    if ( any (na_in_so) )
      {
      # For non contextual vars with n_layers equal to NA:
      # - if no other var has a n_layers, go for automatic estimate
      # - if there is an unique n_layers, apply this value to all
      # - if there is multiple n_layers values, stop
      #
      uniq_vals = unique (state_order$n_layers[(!na_in_so) & (!are_contextual)])
      if (length (uniq_vals) == 0)
        {
        msg_str = list_to_str (state_order$var_names[na_in_so], n_max=10)
        miic_warning ("temporal checks", "the missing number of layers will be ",
          " determined from data for variables ", msg_str, ".")
        }
      else if (length (uniq_vals) > 1)
        {
        miic_error ("temporal checks",
          "some number of layers are missing and they can not be completed",
          " automatically as multiple values are already present.")
        }
      else
        {
        msg_str = list_to_str (state_order$var_names[na_in_so], n_max=10)
        miic_warning ("temporal checks", "the missing number of layers will be ",
          " set to ", uniq_vals, " for variables ", msg_str, ".")
        state_order$n_layers[na_in_so & are_contextual] = uniq_vals
        }
      }
    #
    # Check/fix invalid values:  for contextual vars, n_layers must be 1
    #
    wrongs = ( ( ! is.na (state_order$n_layers) )
             & (state_order$n_layers != 1)
             & (state_order$is_contextual == 1) )
    if ( any (wrongs) )
      {
      msg_str = list_to_str (state_order$var_names[wrongs], n_max=10)
      if (sum (wrongs) == 1)
        miic_warning ("temporal checks", "the variable ", msg_str, ", as",
          " contextual, has an invalid number of layers. It will be set to 1.")
      else
        miic_warning ("temporal checks", "several variables (", msg_str, "), as",
          " contextual, have an invalid number of layers. They will be set to 1.")
      state_order$n_layers[wrongs] = 1
      }
    #
    # Warning if multiple values of n_layers excluding contextual
    #
    uniq_vals = unique (state_order$n_layers[ (!is.na (state_order$n_layers))
                                            & (state_order$is_contextual == 0) ])
    if (length (uniq_vals) > 1)
      {
      msg_str = list_to_str (uniq_vals)
      miic_warning ("temporal checks", "different values (", msg_str,
        ") have be defined for the number of layers.",
        " Such setting is experimental and not recommanded.")
      }
    #
    # Stop if all nb layers == 1
    #
    if (  (!any (is.na (state_order$n_layers)))
       && (all (state_order$n_layers <= 1)) )
      miic_error ("temporal checks", "there must be one variable",
        " at least with a number of layers > 1.")
    }
  #
  # Check state order delta_t (idem as n_layers)
  #
  if ( ! ("delta_t" %in% colnames(state_order)) )
    {
    # Add delta_t column with 0 for contextual, NA otherwise
    #
    state_order$delta_t = NA
    state_order$delta_t[ state_order$is_contextual == 1] = 0
    }
  else
    {
    # Replace NA vals if possible:
    #
    na_in_so = is.na (state_order$delta_t)
    are_contextual = (state_order$is_contextual == 1)
    if ( any (na_in_so & are_contextual) )
      {
      # For contextual, replace NAs with 1
      #
      msg_str = list_to_str (state_order$var_names[na_in_so & are_contextual],
                             n_max=10)
      miic_warning ("temporal checks", "the missing delta t have been",
        " set to 0 for contextual variables (", msg_str, ").")
      state_order$delta_t[na_in_so & are_contextual] = 0
      }
    #
    # Contextual vars done, look if still NAs on not contextual
    #
    na_in_so = is.na (state_order$delta_t)
    if ( any (na_in_so) )
      {
      # For non contextual vars with delta_t equal to NA:
      # - if no other var has a delta_t, go for automatic estimate
      # - if there is an unique delta_t, apply this value to all
      # - if there is multiple delta_t values, stop
      #
      uniq_vals = unique (state_order$delta_t[(!na_in_so) & (!are_contextual)])
      if (length (uniq_vals) == 0)
        {
        msg_str = list_to_str (state_order$var_names[na_in_so], n_max=10)
        miic_warning ("state order", "the missing delta t will be ",
          " determined from data for variables ", msg_str, ".")
        }
      else if (length (uniq_vals) > 1)
        {
        miic_error ("state_order", "the state order contains NAs",
          " for the delta t and it can not be completed",
          " automatically as multiple values are already present.")
        }
      else
        {
        msg_str = list_to_str (state_order$var_names[na_in_so], n_max=10)
        miic_warning ("state order", "the missing delta t will be ",
          " set to ", uniq_vals, " for variables ", msg_str, ").")
        state_order$delta_t[na_in_so & are_contextual] = uniq_vals
        }
      }
    #
    # Check/fix invalid values: for contextual vars, delta_t must be 0
    #
    wrongs = ( (!is.na (state_order$delta_t))
             & (state_order$delta_t != 0)
             & (state_order$is_contextual == 1) )
    if ( any (wrongs) )
      {
      msg_str = list_to_str (state_order$var_names[wrongs], n_max=10)
      if (sum (wrongs) == 1)
        miic_warning ("temporal checks", "the variable ", msg_str, ", as",
          " contextual, has an invalid delta t. It will be set to 0.")
      else
        miic_warning ("temporal checks", "several variables (", msg_str, "), as",
          " contextual, have an invalid delta_t. They will be set to 0.")
      state_order$delta_t[wrongs] = 0
      }
    #
    # Warning if multiple values of delta_t excluding contextual
    #
    uniq_vals = unique (state_order$delta_t[ (!is.na (state_order$delta_t))
                                           & (state_order$is_contextual == 0) ])
    if (length (uniq_vals) > 1)
      {
      msg_str = list_to_str (uniq_vals)
      miic_warning ("temporal checks", "different values (", msg_str,
        ") have be defined for the delta t.",
        " Such setting is experimental and not recommanded.")
      }
    #
    # Stop if all delta t == 0
    #
    if (  (!any (is.na (state_order$delta_t)))
       && (all (state_order$delta_t <= 0)) )
      miic_error ("temporal checks",
                  "there must be one variable at least with a delta t > 0.")
    }
  #
  # Check state order movavg
  #
  if ( ! ("movavg" %in% colnames(state_order)) )
    {
    # Add movavg column with 0 for all vars
    #
    state_order$movavg = 0
    }
  else
    {
    # Replace NA vals by 0
    #
    na_in_so = is.na (state_order$movavg)
    if ( any (na_in_so) )
      {
      msg_str = list_to_str (state_order$var_names[na_in_so], n_max=10)
      miic_warning ("state order", "the missing moving average have been",
        " set to 0 for variables ", msg_str)
      state_order$movavg[na_in_so] = 0
      }
    #
    # Check/fix invalid values: for discrete vars, no moving average
    #
    wrongs = ( (state_order$movavg != 0) & (state_order$var_type == 0) )
    if ( any (wrongs) )
      {
      msg_str = list_to_str (state_order$var_names[wrongs], n_max=10)
      if (sum (wrongs) == 1)
        miic_warning ("temporal checks", "a moving average can not be applied",
          " on a discrete variable ", msg_str, ".")
      else
        miic_warning ("temporal checks", "moving average operations can not",
        " be applied on discrete variables (", msg_str, ").")
      state_order$movavg[wrongs] = 0
      }
    #
    # Check/fix invalid values: for contextual vars, no moving average
    #
    wrongs = ( (state_order$movavg != 0) & (state_order$is_contextual == 1) )
    if ( any (wrongs) )
      {
      msg_str = list_to_str (state_order$var_names[wrongs], n_max=10)
      if (sum (wrongs) == 1)
        miic_warning ("temporal checks", "a moving average can not be applied",
          " on the contextual variable ", msg_str, ".")
      else
        miic_warning ("temporal checks", "moving average operations can not",
        " be applied on contextualvariables (", msg_str, ").")
      state_order$movavg[wrongs] = 0
      }
    #
    # Warning if multiple values of moving average excluding discrete and contextual
    #
    uniq_vals = unique (state_order$movavg[ (!is.na (state_order$movavg))
                                          & (state_order$var_type == 1)
                                          & (state_order$is_contextual == 0) ])
    if (length (uniq_vals) > 1)
      {
      msg_str = list_to_str (uniq_vals)
      miic_warning ("temporal checks", "different values (", msg_str,
        ") have be defined for the moving averages.")
      }
    }
  #
  # Cross checks
  #
  if (  ( "n_layers" %in% colnames(state_order) )
     && ( "delta_t" %in% colnames(state_order) ) )
    {
    if (  (!any (is.na (state_order$n_layers)))
       && (!any (is.na (state_order$delta_t))) )
      {
      t_max = state_order$n_layers * state_order$delta_t
      t_max = t_max[state_order$is_contextual == 0]
      if ( all (t_max < 2) )
        miic_error ("temporal checks", "there must be one variable",
                    " at least with 2 layers and a delta t >= 1.")
      }
    }
  return (state_order)
  }

#-------------------------------------------------------------------------------
# tmiic_check_after_lagging
#-------------------------------------------------------------------------------
# Check the data and the lagged state order: in the state_order,  the var_type
# may need to be  re-evaluated after lagging as some numerical lagged variables
# can have less unique values and are no more considered as discrete
#
# Params :
# - lagged_data: a dataframe, the lagged input data
# - lagged_so  : a dataframe, the lagged state order
# Returns:
# - a dataframe: the lagged state_order, eventually modified
#-------------------------------------------------------------------------------
tmiic_check_after_lagging <- function (lagged_data, lagged_so)
  {
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

  n_unique_vals <- unlist (lapply (lagged_data, function (x) {
    length (unique (x[!is.na(x)] ) ) } ) )
  for (i in 1:nrow (lagged_so))
    {
    if (n_unique_vals[[i]] == 1)
      {
      miic_warning ("lagged data", "the variable ", lagged_so[i, "var_names"],
        " is constant after lagging.")
      lagged_so[i, "var_type"] = 0
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
      lagged_so[i, "var_type"] <- 0
      next
      }
    if (n_unique_vals[[i]] < MIIC_CONTINUOUS_TRESHOLD)
      {
      if (  (! ("var_type_specified" %in% colnames(lagged_so) ) )
         || (!lagged_so[i, "var_type_specified"]) )
        lagged_so[i, "var_type"] <- 0
      }
    }
  lagged_so$var_type_specified <- NULL
  return (lagged_so)
  }

#-------------------------------------------------------------------------------
# tmiic_extract_trajectories
#-------------------------------------------------------------------------------
# Extract the trajectories from a dataframe and return them in a list
# - input_data: a dataframe with the timesteps in the 1st column
#   A new trajectory is identified when timestep < previous timestep
# - check: optional, default=T. Emit warnings when:
#   * there is a gap between 2 consecutive timesteps
#   * the timestep value is not incremented  between 2 consecutive rows
#   * the 1st timesteps of a trajectory is not 1
# Returns:
# - a list: the list of trajectories
#   Note the the timestep information in each trajectory is renumbered from 1
#   to number of timesteps of the trajectory (so no gap, no unchanged timestep)
#-------------------------------------------------------------------------------
tmiic_extract_trajectories <- function (input_data, check=T)
  {
  timesteps = input_data[, 1]
  if ( any ( is.na (timesteps) ) )
    miic_error ("trajectories check", "the timestep column (column 1) contains NA(s)")
  if ( ! all (is.numeric (timesteps)) )
    miic_error ("trajectories check", "the timestep column (column 1) is not integer")
  if ( ! all (round (timesteps, 0) == timesteps) )
    miic_error ("trajectories check", "the timestep column (column 1) is not integer")
  timesteps = as.integer(timesteps)
  timesteps_next = c (timesteps[2:length(timesteps)], 0)
  breaks = which (timesteps_next < timesteps)

  list_ts <- list()
  row_prev = 1
  for ( i in 1:length (breaks) )
    {
    row_new = breaks[[i]]
    list_ts[[i]] = input_data[row_prev:row_new,]
    row_prev = row_new + 1
    }

  if (check)
    {
    no_inc = which (timesteps_next == timesteps)
    if (length (no_inc) > 0)
      miic_warning ("check trajectories", "timestep value unchanged at ",
                    length (no_inc), " position(s)")
    gaps = which (timesteps_next > timesteps + 1)
    if (length (gaps) > 0)
      miic_warning ("check trajectories", "gap in timestep values at ",
                    length (gaps), " position(s)")
    wrong_starts = which ( unlist (lapply (list_ts,
                      FUN=function (x) { return (x[1,1] != 1) } ) ) )
    if (length (wrong_starts) > 0)
      miic_warning ("check trajectories", length (wrong_starts),
        " trajectorie(s) don't start with 1 as first timestep value")
    max_nb_ts = max (unlist (lapply (list_ts, FUN=nrow) ) )
    if (max_nb_ts == 1)
      miic_error ("trajectories check",
                  "all trajectories have only 1 timestep.")
    }
  for ( i in 1:length (list_ts) )
    list_ts[[i]][,1] = 1:nrow (list_ts[[i]])
  return (list_ts)
  }

#-------------------------------------------------------------------------------
# tmiic_group_trajectories
#-------------------------------------------------------------------------------
# Merge a list of trajectories into a dataframe
# - list_ts: the list of trajectories
# - drop_timestep: boolean, FALSE by default. Drop the timestep information
#   (the 1st column) in the returned dataframe
# Returns:
# - a dataframe: dataframe with all the trajectories
#-------------------------------------------------------------------------------
tmiic_group_trajectories = function (list_ts, drop_timestep=FALSE)
  {
  # Pre-allocate the dataframe with the same structure as trajectories
  # and the same number of rows as all the trajectories
  # VOIR data = data[1:n_row_tot,]
  #
  df = list_ts[[1]][FALSE,]
  n_row_tot = sum (unlist (lapply(list_ts, nrow)))
  df <- df[seq_len(n_row_tot),]
  rownames(df) <- NULL

  row_idx = 1
  for (i in 1:length(list_ts) )
    {
    df[row_idx:(row_idx-1+nrow(list_ts[[i]])),] = list_ts[[i]]
    row_idx = row_idx + nrow(list_ts[[i]])
    }
  if (drop_timestep)
    df = df[,-1]
  return (df)
  }

#-------------------------------------------------------------------------------
# tmiic_movavg_onecol
#-------------------------------------------------------------------------------
# Utility function to a apply a moving average over a list
# params:
# - x: the list
# - w: the length of the window
# This moving average is centered, so the first (w-1) %/% 2 and the last
# (w-1) - low_shift items will be filled with NA_real_
#-------------------------------------------------------------------------------
tmiic_movavg_onecol = function (x, w)
  {
  low_shift = (w-1) %/% 2
  high_shift = (w-1) - low_shift
  ret = x
  ret = rep(-1, length(x))
  ret[1:low_shift] = NA_real_
  # print (head (ret))

  start_idx = low_shift+1
  end_idx = length(x) - high_shift
  i = start_idx + 1
  i = start_idx
  for (i in start_idx:end_idx)
    {
    idx_low = i - low_shift
    idx_high = i + high_shift
    ret[i] <- mean (x[idx_low:idx_high], na.action=na.omit)
    }
  # print (head (ret))
  # print (tail (ret))
  ret[(end_idx+1):length(ret)] = NA_real_
  # print (tail (ret))
  return (ret)
  }

#-------------------------------------------------------------------------------
# tmiic_movavg
#-------------------------------------------------------------------------------
# Apply moving averages on data
# - list_ts: a list of dataframe, each item representing a trajectory.
#   Each dataframe must contain the timestep information in the 1st column
#   and the variables in the other columns.
# - movavg: the list of moving average to be applied, optional, NULL by defaut.
#   The length of the movavg list is the number of columns of the dataframes - 1
#   (because the 1st column in dataframes is the timestep).
#   When the movavg item value is >= 2, a moving average using this value as
#   window size is applied on the corresponding column:
#   movavg item 1 is applied data column 2, moavg item 2 to data column 3, ...
# - keep_max_data: boolean flag, optional, FALSE by default
#   When FALSE, the rows containing NA introduced by the moving average(s)
#   are deleted, otherwise when TRUE, the rows are kept
# - verbose_level: integer in the range [0,2], 1 by default. The level of
#   verbosity: 0 = no display, 1 = summary display, 2 = maximum display.
# Returns:
# - list_ts: the list trajectories with moving averages applied
#-------------------------------------------------------------------------------
tmiic_movavg = function (list_ts, movavg=NULL, keep_max_data=F, verbose_level=0)
  {
  if ( is.null (movavg) || all (movavg < 2) )
    return (list_ts)
  if (verbose_level >= 1)
    miic_msg ("Applying moving averages...")
  # Apply movavg on each trajectory and variable of the dataset
  #
  n_vars = ncol(list_ts[[1]])-1
  var_names = colnames (list_ts[[1]])[-1]
  for (i in 1:length(list_ts) )
    for (j in 1:n_vars)
      if (movavg[[j]] >= 2)
        {
        # print (paste0 (j, " => movavg = ", movavg[[j]]))
        list_ts[[i]][,j+1] = tmiic_movavg_onecol (list_ts[[i]][,j+1], movavg[[j]])
        if (verbose_level == 2)
          miic_msg ("- ", var_names[[j]], ": moving average of window size ",
                    movavg[[j]], " applied")
        }
  #
  # Remove starting and ending rows where NAs were introduced
  #
  if (!keep_max_data)
    {
    movavg_max = max(movavg)
    low_shift = (movavg_max-1) %/% 2
    high_shift = (movavg_max-1) - low_shift
    start_idx = 1
    if (low_shift > 0)
      start_idx = start_idx + low_shift
    i = 1
    for (i in 1:length(list_ts) )
      {
      end_idx = nrow(list_ts[[i]]) - high_shift
      list_ts[[i]] = list_ts[[i]][start_idx:end_idx,]
      }
    }
  return (list_ts)
  }

#-----------------------------------------------------------------------------
# tmiic_ajust_window_for_nb_samples
#-----------------------------------------------------------------------------
# Reduce the window size (n_layers or delta_t) if the foreseen n_layers and
# delta_t would lead to too few samples after data lagging
# params:
# - list_ts: a list of dataframe, each item representing a trajectory.
# - n_layers: a list, the n_layers in the state_order column
# - delta_t: a list, the delta_t in the state_order column
# - reduced_param: a string, can be "n_layers" or "delta_t". Indicates with
#   parameter will be reduced if the number of samples is too small
# - verbose: boolean, if TRUE, display a message if the window size is reduced
# returns:
# - a list: the n_layers or delta_t, depending of the reduced_param value.
#   The value are possibly decreased to reduce the window size
#-----------------------------------------------------------------------------
tmiic_ajust_window_for_nb_samples <- function (list_ts, n_layers, delta_t,
                                               reduced_param, verbose)
  {
  tau_per_var = (n_layers - 1) * delta_t
  tau_max = max (tau_per_var)
  ts_lengths = unlist ( lapply (list_ts, nrow) )
  tot_ts = sum (ts_lengths)
  nb_samples = sum ( unlist (lapply (ts_lengths, FUN=function (x) {
                                            max (0, x - tau_max) } ) ) )
  target = min (1000, tot_ts / 10)

  if (nb_samples < target)
    {
    # Look for the best value to reach the target recursivly
    # At each iteration, we keep the half part where the best value is
    # until we can not divide in half further
    #
    recurs_eval = function (target, tau_low, tau_high, ts_lengths)
      {
      if (tau_high - tau_low <= 1)
        return (tau_low)
      tau = round ( (tau_low + tau_high) / 2, 0)
      nb_samples = sum ( unlist (lapply (ts_lengths, FUN=function (x) {
                                            max (0, x - tau) } ) ) )
      if (nb_samples >= target)
        tau_ret = recurs_eval (target, tau, tau_high, ts_lengths)
      else
        tau_ret = recurs_eval (target, tau_low, tau, ts_lengths)
      return (tau_ret)
      }
    tau_red = recurs_eval (target, 1, tau_max, ts_lengths)
    #
    # Max time steps back in time found, (try to) reduce n_layers or delta_t
    #
    if (reduced_param == "n_layers")
      n_layers[tau_per_var > tau_red] = max ( 2,
        floor (tau_red / delta_t[tau_per_var > tau_red]) + 1)
    else # reduce delta_t
      delta_t[tau_per_var > tau_red] = max ( 1,
        floor (tau_red / (n_layers[tau_per_var > tau_red] - 1)) )
    #
    # Check the effect of reduction and feed back to user
    #
    if (reduced_param == "n_layers")
      fixed_param = "delta_t"
    else
      fixed_param = "n_layers"
    tau_max_red = max ( (n_layers - 1) * delta_t )
    nb_samples_red = sum ( unlist (lapply (ts_lengths, FUN=function (x) {
                                          max (0, x - tau_max_red) } ) ) )
    if (nb_samples_red <= 0)
      miic_error ("temporal parameters estimation",
        "with the values supplied in ", fixed_param,
        ", no valid ", reduced_param, " can be estimated.")
    else if ( (nb_samples_red < target) && (nb_samples_red == nb_samples) )
      miic_warning ("temporal parameters estimation",
        "with the estimated or supplied temporal parameters",
        ", the number of usable samples will be ", nb_samples,
        ". Consider to specify manually n_layers and delta_t.")
    else if (nb_samples_red < target)
      miic_warning ("temporal parameters estimation",
        "the ", reduced_param, " parameter has been reduced",
        " to increase the number of samples. However,",
        " the number of usable samples will still only be ", nb_samples_red,
        ". Consider to specify manually n_layers and delta_t.")
    else if (verbose)
      miic_msg ("- The ", reduced_param, " parameter has been reduced ",
        " to increase the number of samples.")
    }

  if (reduced_param == "n_layers")
    return (n_layers)
  else
    return (delta_t)
  }

#-------------------------------------------------------------------------------
# tmiic_estimate_dynamic
#-------------------------------------------------------------------------------
# Estimate tau (the number of total timesteps back to cover the dynamic,
# the number of layers and delta t parameters from the data
# - list_ts: list of dataframe, each item representing a trajectory.
#   Each dataframe must contain the timestep information in the 1st column
#   and the variables in the other columns.
# - state_order: the state_order dataframe. This state_order is expected
#   having being checked by the temporal check functions of inputs.
#   The rows in the state_order must be ordered as the columns in the data.
#   It must contain the var_type, is_contextual, n_layers and delta_t columns.
#   There can be NAs in n_layers and delta_t for continuous and non contextual
#   variables.
# - max_nodes: maximum number of nodes in the inferred time unfolded graph,
#   optional, 50 by default
# - verbose_level: integer in the range [0,2], 1 by default. The level of
#   verbosity: 0 = no display, 1 = summary display, 2 = maximum display.
#-----------------------------------------------------------------------------
tmiic_estimate_dynamic <- function (list_ts, state_order, max_nodes=50,
                                    verbose_level=1)
  {
  #
  # If n_layers and delta_t all defined, nothing to do
  #
  if (  (! any (is.na (state_order$n_layers) ) )
     && (! any (is.na (state_order$delta_t ) ) ) )
    return (state_order)
  #
  # We are going to estimate to temporal dynamic because we need to fill out
  # the missing values in n_layers and/or delta_t.
  #
  # After the checks done on the state order, we know that
  # the missing values are only for continuous and non contextual.
  # In addition, we know that all values for continuous and non contextual
  # are NAs (otherwise the checks would have completed the NAs by
  # generalizing the known values)
  #
  if (verbose_level == 2)
    miic_msg ("Estimating the temporal dynamic...")
  n_ts = length (list_ts)
  n_vars_tot = ncol (list_ts[[1]]) - 1
  n_vars_ctx = sum (state_order$is_contextual)
  n_vars_lag = n_vars_tot - n_vars_ctx
  #
  # Remove timestep, contextual and discrete variables
  #
  for (ts_idx in 1:n_ts)
    {
    if (nrow (list_ts[[ts_idx]]) == 1)
      miic_warning ("Dynamic estimation", "trajectory ", ts_idx,
                    " with only 1 timestep is ignored for dynamic estimation")
    list_ts[[ts_idx]] = list_ts[[ts_idx]][, c(F, ( (!state_order$is_contextual)
                                                 & (state_order$var_type) ) )]
    }
  #
  # Compute mean alpha per variable
  #
  n_vars = ncol (list_ts[[1]])
  var_names = colnames(list_ts[[1]])
  length_to_test = min (unlist (lapply (list_ts, FUN=function (x) {
    ifelse ( nrow(x) <= 1, NA, nrow(x) ) } ) ), na.rm=T)
  alphas_per_var = rep (NA, n_vars)
  taus_per_var = rep (NA, n_vars)
  var_idx = 1
  for (var_idx in 1:n_vars)
    {
    alphas_per_ts = rep (NA, n_ts)
    ts_idx = 1
    for (ts_idx in 1:n_ts)
      {
      if (nrow (list_ts[[ts_idx]]) == 1)
        next
      acf_res = acf (list_ts[[ts_idx]][,(var_idx)], na.action=na.pass,
                     lag.max=length_to_test-1, plot=F)
      if ( all (is.na(acf_res$acf) ) )
        next
      acf_vanish = which (acf_res$acf[,1,1] < 0.05)
      if ( length (acf_vanish) == 0 )
        acf_vanish = length_to_test
      lag_vanish = acf_res$lag[min (acf_vanish), 1, 1]
      lag_4_alpha = max ( 1, round (lag_vanish / 2) )
      alphas_per_ts[[ts_idx]] = acf_res$acf[lag_4_alpha+1,1,1] ^ (1/lag_4_alpha)
      }
    alphas_per_var[[var_idx]] = mean (alphas_per_ts, na.rm=T)
    taus_per_var[[var_idx]] = round ( (1+alphas_per_var[[var_idx]])
                                    / (1-alphas_per_var[[var_idx]]) )
    }
  if (verbose_level == 2)
    {
    miic_msg ("Tau per variable:")
    for (i in 1:length (var_names))
      miic_msg ("- ", var_names[[i]], ": ", taus_per_var[[i]])
    }
  #
  # Compute alphas range and deduce taus range
  #
  # print ("alphas_per_var")
  # print (alphas_per_var)
  # print ("taus")
  # print (unlist (lapply (alphas_per_var, FUN=function(x) {
  #   return ( (1+x) / (1-x) ) })))
  # print ("mean taus")
  # print (mean (unlist (lapply (alphas_per_var, FUN=function(x) {
  #   return ( (1+x) / (1-x) ) })), na.rm=T))

  tau_min  = max ( 1, min (taus_per_var, na.rm=T) )
  tau_mean = max ( 1, round (mean (taus_per_var, na.rm=T), 0) )
  tau_max  = min ( length_to_test, max  (taus_per_var, na.rm=T) )
  tau_max_kept = min (length_to_test, tau_max, tau_mean * 2)
  tau = tau_max_kept
  if (verbose_level >= 1)
    miic_msg ("Automatic estimation of parameters:\n",
      "- Relaxation times goes from ", tau_min, " to ", tau_max,
      " with a mean of ", tau_mean, " => tau max considered = ", tau_max_kept)
  #
  # We know tau : the average maximum time steps back in time to use for the
  # temporal discovery. Now estimate the number of layers 'n_layers'
  # and/or number of time steps between two layers 'delta_t'
  #
  if (  all (!is.na (state_order$n_layers)) ) # n_layers known => NAs in delta_t
    {
    state_order$delta_t[is.na(state_order$delta_t)] = max ( 1,
        ceiling (tau / (state_order$n_layers[is.na(state_order$delta_t)] - 1)) )

    state_order$delta_t = tmiic_ajust_window_for_nb_samples (list_ts,
      state_order$n_layers, state_order$delta_t, reduced_param="delta_t",
      verbose=(verbose_level >= 1) )

    uniq_n_layers = unique (state_order$n_layers[ (state_order$is_contextual == 0)
                                                & (state_order$var_type == 1)] )
    uniq_delta_t = unique (state_order$delta_t[ (state_order$is_contextual == 0)
                                              & (state_order$var_type == 1)] )
    if (verbose_level >= 1)
      {
      if (length (uniq_n_layers) == 1)
        miic_msg ("- As the number of layers was defined to ", uniq_n_layers,
          ", the only parameter tuned is the delta t set to ", uniq_delta_t, ".")
      else
        miic_msg ("- As multiple values of layers were present (",
            list_to_str (uniq_n_layers), "), the delta t have been set,",
            " respectively to ", list_to_str (uniq_delta_t), ".")
      }
    }
  else if (  all (!is.na (state_order$delta_t)) ) # delta_t known => NAs in n_layers
    {
    # To determine the layers, we compute the max number of layers considering
    # the maximum number of nodes in the final grpah
    #
    n_layers_max = max (2, floor ( (max_nodes - n_vars_ctx) / n_vars_lag ) )
    #
    # The final number of layers will (tau / delta_t) + 1 unless if greater
    # than the max number of layers
    #
    state_order$n_layers[is.na(state_order$n_layers)] = min (n_layers_max,
      ceiling (tau / state_order$delta_t[is.na(state_order$n_layers)]) + 1)

    state_order$n_layers = tmiic_ajust_window_for_nb_samples (list_ts,
      state_order$n_layers, state_order$delta_t, reduced_param="n_layers",
      verbose=(verbose_level >= 1) )

    uniq_n_layers = unique (state_order$n_layers[ (state_order$is_contextual == 0)
                                                & (state_order$var_type == 1)] )
    uniq_delta_t = unique (state_order$delta_t[ (state_order$is_contextual == 0)
                                              & (state_order$var_type == 1)] )
    if (verbose_level >= 1)
      {
      if (length (uniq_delta_t) == 1)
        miic_msg ("- As the value of delta t was defined to ", uniq_delta_t,
            ", the only parameter tuned is the number of layers set to ",
            uniq_n_layers, ".")
      else
        miic_msg ("- As multiple values of delta t were present (",
            list_to_str (uniq_delta_t), "), the number of layers have been set,",
            " respectively to ", list_to_str (uniq_n_layers), ".")
      }
    }
  else
    {
    # Both n_layers and delta_t need to be estimated automatically
    #
    delta_t = 1
    if ( (tau + 1) * n_vars_lag + n_vars_ctx <= max_nodes)
      {
      # If when using delta_t = 1, the n_layers (= tau + 1) does not lead to
      # a graph with a total number of nodes > max => OK, nothing more to do
      #
      n_layers = tau + 1
      }
    else
      {
      # We need reduce the number of layers to respect the maximum nodes number
      # and increase de delta t to still cover all the dynamic tau.
      # => Compute the max number of layers and deduce the delta t
      #
      n_layers = max (2, floor ( (max_nodes - n_vars_ctx) / n_vars_lag ) )
      if (n_layers > 2)
        {
        delta_t = max (1, ceiling ( tau / (n_layers-1)  ) )
        tau = (n_layers - 1) * delta_t
        }
      else
        delta_t = tau
      }

    state_order$n_layers[is.na(state_order$n_layers)] = n_layers
    state_order$delta_t[is.na(state_order$delta_t)] = delta_t

    state_order$delta_t = tmiic_ajust_window_for_nb_samples (list_ts,
      state_order$n_layers, state_order$delta_t, reduced_param="delta_t",
      verbose=(verbose_level >= 1) )
    delta_t = unique (state_order$delta_t[ (state_order$var_type == 1)
                                         & (state_order$is_contextual == 0) ])

    if (verbose_level >= 1)
      miic_msg ("- For a final graph with a target of ", max_nodes,
        " nodes having ", n_vars_lag, " lagged variables",
        ifelse (n_vars_ctx > 0, paste0 ("\n  and ", n_vars_ctx, " contextual variables"), ""),
        ":\n  ", n_layers,  " layers spaced by ", delta_t, " timesteps",
        ", dynamic covered goes over t, t-", delta_t,
        ifelse (n_layers > 3, ", ...", ""),
        ifelse (n_layers > 2, paste0 (", t-", tau), "") )
    }

  return (state_order)
  }

#-------------------------------------------------------------------------------
# estimateTemporalDynamic
#-------------------------------------------------------------------------------
#' Estimation of the temporal causal discovery parameters
#'
#' @description This function estimates the number of layers and number of
#' time steps between each layer that are needed to cover the dynamic of a
#' temporal dataset when reconstructing a temporal causal graph.
#' Using autocorrelation decay, the function computes the average relaxation
#' time of the variables and, in regard of a maximum number of nodes, deduces
#' the number of layers and number of time steps between each layer to be used.
#'
#' @param input_data [a data frame]
#' A data frame that contains the observational data.\cr
#' The expected dataframe layout is variables as columns and
#' timeseries/timesteps as rows.
#' The timestep information must be supplied in the first column and,
#' for each timeseries, be consecutive (increment of 1) and in ascending order.
#' Multiple trajectories can be provided, the function will consider that a
#' new trajectory starts each time a smaller timestep than the one of the
#' previous row is encountered.
#'
#' @param state_order [a data frame] An optional data frame providing extra
#' information about variables. It must have d rows where d is the number of
#' input variables, excluding the time step one.\cr
#'
#' The following structure (named columns) is expected:\cr
#'
#' "var_names" (required) contains the name of each variable as specified
#' by colnames(input_data), excluding the time steps column.
#'
#' "var_type" (optional) contains a binary value that specifies if each
#' variable is to be considered as discrete (0) or continuous (1).
#' Discrete variables will be excluded from the temporal dynamic estimation.
#'
#' "is_contextual" (optional) contains a binary value that specifies if a
#' variable is to be considered as a contextual variable (1) or not (0).
#' Contextual variables will be excluded from the temporal dynamic estimation.
#'
#' "movavg" (optional) contains an integer value that specifies the size of
#' the moving average window to be applied to the variable.
#' Note that if "movavg" column is present in the \emph{state_order},
#' its values will overwrite the function parameter.
#'
#' @param movavg [an integer] Optional, NULL by default.\cr
#' When an integer>= 2 is supplied, a moving average operation is applied
#' to all the non discrete and not contextual variables. If no \emph{state_order}
#' is provided, the discrete/continuous variables are deduced from the input
#' data. If you want to apply a moving average only on specific columns,
#' consider to use a \emph{movavg} column in the \emph{state_order} parameter.
#'
#' @param max_nodes [a positive integer] The maximum number of nodes in the
#' final temporal causal graph. The more nodes are allowed in the temporal
#' causal discovery, the more precise will be the discovery but at the cost
#' of longer execution time. The default is set to 50 for a fast causal
#' discovery. On recent computers, values up to 200 or 300 nodes are usually
#' possible (depending on the number of trajectories and time steps in the
#' input data).
#'
#' @param verbose_level [an integer value in the range [0,2], 1 by default]
#' The level of verbosity: 0 = no display, 1 = summary display, 2 = full display.
#'
#' @return A list that contains:
#' \itemize{
#'  \item{n_layers:}{the number of layers}
#'  \item{delta_t:}{the number of time steps between the layers}
#' }
#'
#' @export
#-------------------------------------------------------------------------------
estimateTemporalDynamic <- function (input_data, state_order=NULL, movavg=NULL,
                                     max_nodes=50, verbose_level=1)
  {
  input_data = check_input_data (input_data, "TS")
  state_order = check_state_order (input_data, state_order, "TS")
  state_order$n_layers = NULL
  state_order$delta_t = NULL
  state_order = tmiic_check_state_order_part1 (state_order)
  list_ret = tmiic_check_parameters (state_order = state_order,
                                     params = list(),
                                     n_layers = NULL,
                                     delta_t = NULL,
                                     movavg = movavg,
                                     keep_max_data = F,
                                     max_nodes = max_nodes)
  state_order = tmiic_check_state_order_part2 (list_ret$state_order)

  list_ts = tmiic_extract_trajectories (input_data)
  list_ts = tmiic_movavg (list_ts, state_order$movavg, verbose_level=verbose_level)

  state_order = tmiic_estimate_dynamic (list_ts, state_order, max_nodes=max_nodes,
                                        verbose_level=verbose_level)
  n_layers = unique (state_order$n_layers[ (state_order$var_type == 1)
                                         & (state_order$is_contextual == 0) ])
  delta_t = unique (state_order$delta_t[ (state_order$var_type == 1)
                                         & (state_order$is_contextual == 0) ])
  return ( list ("n_layers"=n_layers, "delta_t"=delta_t))
  }

#-------------------------------------------------------------------------------
# tmiic_check_other_df
#-------------------------------------------------------------------------------
# This function is designed to be called after all the temporal checks on
# parameters and state_order have been done.
# Using the state order completely filled for each variable, the checks
# here will verify that true edges or edges in black box are valid.
#
# Params :
# - df: a dataframe, a true edges or black box dataframe, with 3 columns:
#   origin node, destination node, lag)
# - df: the name of the data frame, used to display warnings
# - state_order: a dataframe, the state order comptletely checked with
#   the column var_type, n_layers and delta_tau fully determined
# Returns:
# - a dataframe: the dataframe, eventually with invalid rows filtered out
#-------------------------------------------------------------------------------
tmiic_check_other_df <- function (df, df_name, state_order)
  {
  if ( (is.null (df)) || (nrow (df) <= 0) )
    return (df)

  for (i in 1:nrow (df))
    {
    orig_node_idx <- which (state_order$var_names == df[i, 1])
    if (state_order[orig_node_idx, "is_contextual"] == 1)
      {
      if ( ! is.na (df[i, 3]) )
        {
        miic_warning (paste0 (df_name, " temporal checks"),
          "the edge at row ", i, " (", df[i, 1], " - ", df[i, 2], ", lag",
          df[i, 3], ") does not exist in the lagged graph and will be ignored.")
        df[i, 3] = -1
        next
        }
      }
    else
      {
      possible_lags = seq (from=0, by=state_order[orig_node_idx, "delta_t"],
                            length=state_order[orig_node_idx, "n_layers"])
      if (! (df[i, 3] %in% possible_lags) )
        {
        miic_warning (paste0 (df_name, " temporal checks"),
          "the edge at row ", i, " (", df[i, 1], " - ", df[i, 2], ", lag",
          df[i, 3], ") does not exist in the lagged graph and will be ignored.")
        df[i, 3] = -1
        next
        }
      }

    dest_node_idx <- which (state_order$var_names == df[i, 2])
    if (state_order[dest_node_idx, "is_contextual"] == 1)
      {
      miic_warning (paste0 (df_name, " temporal checks"),
        "the edge at row ", i, " (", df[i, 1], " - ", df[i, 2], ", lag",
        df[i, 3], ") is invalid. A contextual node can not be the",
        " destination node, the row will be ignored.")
      df[i, 3] = -1
      next
      }
    }
  df = df[ !(df$lag == -1) ,]
  return (df)
  }
