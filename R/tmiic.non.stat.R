#*******************************************************************************
# Filename   : tmiic.non.stat.R                Creation date: 11 may 2023
#
# Description: Utility functions for temporal MIIC non stationary
#
# Author     : Franck SIMON
#
# Changes history:
# - 31 jan  2024 : initial version
#*******************************************************************************

#===============================================================================
# FUNCTIONS
#===============================================================================
# tmiic_non_stat_check_str_value_against_data
#-------------------------------------------------------------------------------
# Check if a string value is compatible and can be used as filter for the data.
# Note that a valid string value is expected to be enclosed by ' or "
# (an x int will stored as x in the string variable and
# a x string will be stored 'x' or "x" in the variable)
#
# Parameters:
# - data: a list, the values in data
# - value: a string, the value to be tested
# - operator: operator used to find matching data
# Returns:
# - NULL : the value is not compatible with the data
# - a list of integer : the value is compatible with the data and the list
#   contains the positions matching the filter with the value
#-------------------------------------------------------------------------------
tmiic_non_stat_check_str_value_against_data <- function (
  possible_values, value, operator="==")
  {
  # Value must not be empty
  #
  if ( (nchar (value) <= 0) || (nchar (trimws (value) ) <= 0) )
    return (NULL)
  value <- trimws (value)
  #
  # Value must have the same type as data
  #
  if (is.logical (possible_values) )
    {
    if ( ! value %in% c("T", "F", "TRUE", "FALSE") )
      return (NULL)
    value_casted <- as.logical (value)
    }
  else if (is.numeric (possible_values) )
    {
    value_casted <- suppressWarnings( as.numeric(value) )
    if ( is.na (value_casted) )
      return (NULL)
    }
  else # assumed to be character
    {
    if ( ! is.character (value) )
      return (NULL)
    if ( ! grepl ("^['\"].*['\"]$", value) )
      return (NULL)
    value_casted <- substr (value, 2, nchar(value) - 1)
    }
  #
  # Find the positions of the value in data
  # Possible operators : ">=x", ">x", "<=x", "<x", "!=x", "==x"
  #
  if (operator == ">=")
    positions <- which (possible_values >= value_casted)
  else if (operator == ">")
    positions <- which (possible_values > value_casted)
  else if (operator == "<=")
    positions <- which (possible_values <= value_casted)
  else if (operator == "<")
    positions <- which (possible_values < value_casted)
  else if (operator == "!=")
    positions <- which (possible_values != value_casted)
  else if (operator == "==")
    positions <- which (possible_values == value_casted)
  else # should not occur as operator should have been checked before
    miic_error ("parameters", "wrong operator ", operator, " given to tmiic_non_stat_check_str_value_against_data.")

  return (positions)
  }

#-------------------------------------------------------------------------------
# tmiic_non_stat_check_var_interest_condition
#-------------------------------------------------------------------------------
# Check the condition on the variable of interest. Valid conditions are :
# "x;y" : the values in data changes from x to y
# ">=x", ">x", "<=x", "<x", "!=x", "==x" : the values in data reaches a point
# where value is >= > <= < != == to x
# For logical, only "x;y", "!=x" and "==x" are allowed
#
# Parameters:
# - input_data: the data frame with input data
# - var_interest: the name of the variable of interest
# - var_interest_condition: the condition on the variable of interest
# Returns:
# - a string, the variable of interest if valid. NULL otherwise
#-------------------------------------------------------------------------------
tmiic_non_stat_check_var_interest_condition <- function (
  input_data, var_interest, var_interest_condition)
  {
  if (  is.null(var_interest_condition)
     || (length (var_interest_condition) != 1)
     || is.na(var_interest_condition)
     || (!is.character(var_interest_condition)) )
    {
    miic_warning ("parameters", "the variable of interest condition ",
      list_to_str (var_interest_condition), " is invalid.",
      " Alignement on ", var_interest, " will be ignored.")
    return (NULL)
    }
  #
  # Ensure we can determine the data type of the variable
  #
  if (  ( ! is.logical (input_data[,var_interest]) )
     && ( ! is.numeric (input_data[,var_interest]) )
     && ( ! is.character (input_data[,var_interest]) ) )
    {
    miic_warning ("parameters", "Unable to check the data type for the variable of interest condition ",
      list_to_str (var_interest_condition), ".",
      " Alignement on ", var_interest, " will be ignored.")
    return (NULL)
    }
  #
  # "x,y" case, condition is when the value in data changes from x to y
  #
  if ( grepl (",", var_interest_condition, fixed=T) )
    {
    values <- trimws (strsplit (var_interest_condition, ",", fixed=T)[[1]])
    if (length (values) != 2)
      {
      miic_warning ("parameters", "the variable of interest condition ",
        list_to_str (var_interest_condition), " is invalid.",
        " Alignement on ", var_interest, " will be ignored.")
      return (NULL)
      }
    if (length (values) > 2)
      {
      miic_warning ("parameters", "the variable of interest condition ",
                    list_to_str (var_interest_condition), " is invalid.",
                    " Only the two first values will be kept.")
      values <- values[1:2]
      }
    pos1 <- tmiic_non_stat_check_str_value_against_data (
      possible_values=input_data[,var_interest], value=values[[1]])
    pos2 <- tmiic_non_stat_check_str_value_against_data (
      possible_values=input_data[,var_interest], value=values[[2]])
    if ( is.null (pos1) || is.null (pos2) )
      {
      miic_warning ("parameters", "the variable of interest condition ",
        list_to_str (var_interest_condition), " is invalid.",
        " Alignement on ", var_interest, " will be ignored.")
      return (NULL)
      }
    if ( (length (pos1) <= 0) || (length (pos2) <= 0) )
      {
      miic_warning ("parameters", "the variable of interest condition ",
        list_to_str (var_interest_condition), " is not in data.",
        " Alignement on ", var_interest, " will be ignored.")
      return (NULL)
      }
    transition_found <- unlist (lapply (pos1, FUN=function(x) { return ( (x+1) %in% pos2) } ) )
    if ( ! any (transition_found) )
      {
      miic_warning ("parameters", "the variable of interest condition ",
        list_to_str (var_interest_condition), " has no match in data.",
        " Alignement on ", var_interest, " will be ignored.")
      return (NULL)
      }
    #
    # Condition is valid
    #
    return (var_interest_condition)
    }
  #
  # Other cases : ">=x", ">x", "<=x", "<x", "!=x", "==x" :
  #
  if (  ( grepl ("^>=", var_interest_condition) )
     || ( grepl ("^<=", var_interest_condition) )
     || ( grepl ("^!=", var_interest_condition) )
     || ( grepl ("^==", var_interest_condition) ) )
    {
    operator <- substr (var_interest_condition, 1, 2)
    one_value <- substr (var_interest_condition, 3, nchar(var_interest_condition))
    }
  else if (  ( grepl ("^>", var_interest_condition) )
     || ( grepl ("^<", var_interest_condition) ) )
    {
    operator <- substr (var_interest_condition, 1, 1)
    one_value <- substr (var_interest_condition, 2, nchar(var_interest_condition))
    }
  else
    {
    miic_warning ("parameters", "the variable of interest condition ",
      list_to_str (var_interest_condition), " is invalid.",
      " Alignement on ", var_interest, " will be ignored.")
    return (NULL)
    }
  #
  # For logical, only != and == operator are allowed
  #
  if ( is.logical (input_data[,var_interest]) )
    {
    if (! (operator %in% c("!=", "==") ))
      {
      miic_warning ("parameters", "the variable of interest condition ",
        list_to_str (var_interest_condition), " is invalid for a logical variable ",
        "(only == and != are allowed).",
        " Alignement on ", var_interest, " will be ignored.")
      return (NULL)
      }
    }
  #
  # value must not be empty
  #
  positions <- tmiic_non_stat_check_str_value_against_data (
    possible_values=input_data[,var_interest], value=one_value, operator=operator)
  if ( is.null (positions) )
    {
    miic_warning ("parameters", "the variable of interest condition ",
      list_to_str (var_interest_condition), " is invalid.",
      " Alignement on ", var_interest, " will be ignored.")
    return (NULL)
    }
  if (length (positions) <= 0)
    {
    miic_warning ("parameters", "the variable of interest condition ",
      list_to_str (var_interest_condition), " is not in data.",
      " Alignement on ", var_interest, " will be ignored.")
    return (NULL)
    }
  #
  # Condition is valid
  #
  return (var_interest_condition)
  }

#-------------------------------------------------------------------------------
# tmiic_non_stat_check_window_position
#-------------------------------------------------------------------------------
# Check the positioning around the variable of interest. Valid positions are :
# "start", "end" and, in addition if we use a variables of interest, "around".
#
# Parameters:
# - var_interest: optional, a string, the variable of interest, if any.
# - window_position: a string, the positioning of the temporal window
# Returns:
# - a string, the window positioning if valid, NULL otherwise
#-------------------------------------------------------------------------------
tmiic_non_stat_check_window_position <- function (var_interest, window_position)
  {
  if ( is.null (var_interest) )
    possibles <- c("start", "end")
  else
    possibles <- c("around", "start", "end")

  if (  is.null (window_position)
     || (length (window_position) != 1)
     || is.na(window_position)
     || (!is.character(window_position))
     || (!(window_position %in% possibles)) )
    {
    msg_str <- paste0 ("the window position parameter ",
        list_to_str (window_position),  " is invalid. Possible values are: ",
        "'start', 'end' and, only if a variable of interest is specified, 'around'.",
        " The default value 'start' will be used.")
    miic_warning ("parameters", msg_str)
    return ("start")
    }
  return (window_position)
  }

#-------------------------------------------------------------------------------
# tmiic_find_var_interest_position
#-------------------------------------------------------------------------------
# Find for each trajectory the position related to the variable of interest
# inputs:
# - list_traj: the list of trajectories
# - params: a list, the list of MIIC params.
# returns;
# - list with  the first matching position, for each trajectory, with the
#   condition on the variable of interest (NA if no match)
#-------------------------------------------------------------------------------
tmiic_find_var_interest_position <- function (list_traj, params)
  {
  # Identify the test to apply
  #
  operator <- NULL
  if ( grepl(",", params$var_interest_condition, fixed=T) )
    two_values <- strsplit (params$var_interest_condition, ",", fixed=T)[[1]]
  else
    {
    if (  ( grepl("^>=", params$var_interest_condition) )
       || ( grepl("^<=", params$var_interest_condition) )
       || ( grepl("^!=", params$var_interest_condition) )
       || ( grepl("^==", params$var_interest_condition) ) )
      {
      operator <- substr (params$var_interest_condition, 1, 2)
      one_value <- substr (params$var_interest_condition, 3, nchar(params$var_interest_condition))
      }
    else if (  ( grepl ("^>", params$var_interest_condition) )
            || ( grepl ("^<", params$var_interest_condition) ) )
      {
      operator <- substr (params$var_interest_condition, 1, 1)
      one_value <- substr (params$var_interest_condition, 2, nchar(params$var_interest_condition))
      }
    else
      miic_error ("trajectories alignemnt", "wrong condition ", params$var_interest_condition, ".")
    }
  #
  # Identify in each trajectory where the condition is met
  #
  n_traj <- length(list_traj)
  pos_interests <- rep (NA_integer_, n_traj) # NA_integer as we are looking for positions
  for (i in 1:n_traj)
    {
    values <- list_traj[[i]][,params$var_interest]
    if ( is.null(operator) )
      {
      pos1 <- tmiic_non_stat_check_str_value_against_data (
        possible_values=values, value=two_values[[1]])
      pos2 <- tmiic_non_stat_check_str_value_against_data (
        possible_values=values, value=two_values[[2]])
      # NB: use pos2 as we want to return the position after the value change
      criteria_found <- pos2[ pos2 %in% (pos1 + 1) ]
      }
    else
      criteria_found <- tmiic_non_stat_check_str_value_against_data (
        possible_values=values, value=one_value, operator=operator)

    if (length (criteria_found) <= 0)
      miic_warning ("trajectories alignement", "the trajectory ", i,
                    " can not be aligned on the variable of interest and will be discarded.")
    else
      pos_interests[[i]] <- criteria_found[[1]]
    }
  return (pos_interests)
  }

#-------------------------------------------------------------------------------
# tmiic_non_stat_align_data
#-------------------------------------------------------------------------------
# Prepare the input_data for the non stationary causal discovery
# The trajectories after alignment will all have the same number of time steps
# and be aligned on start, end or a value change of the variable of interest.
# Note that the window position can be modified by the function is it is not
# possible (i.e.: window position is 'around' but no time step is available
# after the condition => the window position will be turned into 'end').
#
# Params:
# - list_traj: the list of trajectories
# - params: a list, the list of MIIC params.
# Returns: a list with 3 items
# - list_traj: the list of trajectories matching the criteria, if any,
#   aligned, cut if needed to have an equal number of time steps
# - params: the list of MIIC params (window_position can be modified)
# - pos_interest: the index of the position of interest
#   (i.e.: 1 for alignment on 'start', last row for alignment on 'end')
#-------------------------------------------------------------------------------
tmiic_non_stat_align_data <- function (list_traj, params)
  {
  # Delete trajectories with one time step
  #
  ts_lengths <- unlist (lapply (list_traj, nrow))
  ts_ok <- (ts_lengths > 1)
  if (any (!ts_ok) )
    {
    list_traj <- list_traj[ts_ok]
    ts_lengths <- ts_lengths[ts_ok]
    if ( sum (!ts_ok) == 1)
      miic_warning ("trajectories alignement", sum(!ts_ok),
          " trajectory has only 1 time step, it has been discarded.")
    else
      miic_warning ("trajectories alignement", sum(!ts_ok),
          " trajectories have only 1 time step, they have been discarded.")
    }
  if ( length(list_traj) <= 0 )
    miic_error ("trajectories alignement", "no trajectory available to perform the alignment.")
  min_length <- min (ts_lengths)
  #
  # Alignment depends if a var_interest is supplied
  #
  if ( is.null (params$var_interest) )
    {
    # If no alignment on a variable of interest is requested, align the
    # trajectories on start or end with the same (minimal) length
    #
    if (params$window_position == "start")
      {
      if ( min_length != max (ts_lengths) )
        list_traj <- lapply (list_traj, FUN=function (x) {
          return (x[1:min_length,]) })
      aligned_pos_interest <- 1
      }
    else # window_position == "end"
      {
      if ( min_length != max (ts_lengths) )
        list_traj <- lapply (list_traj, FUN=function (x) {
          return( x[(nrow(x)-min_length+1):nrow(x),]) } )
      aligned_pos_interest <- nrow ( list_traj[[1]] )
      }
    }
  else # var_interest != NULL
    {
    # If alignment on a variable of interest is requested
    # use the variable of interest condition to determine positions
    #
    pos_interests <- tmiic_find_var_interest_position (list_traj, params)
    #
    # Filter out trajectories not matching the condition
    #
    ts_ok <- ( ! is.na(pos_interests) )
    list_traj <- list_traj[ts_ok]
    ts_lengths <- ts_lengths[ts_ok]
    pos_interests <- pos_interests[ts_ok]
    if ( length(list_traj) <= 0 )
      miic_error ("trajectories alignement",
        "no trajectory matches the condition ", params$var_interest_condition,
        " on the variable of interest ", params$var_interest, ".")
    #
    # We have the position to align each trajectory. Look how much time steps
    # we have before and after for all trajectories matching the condition
    #
    ts_after <- ts_lengths - pos_interests
    max_before <- max (pos_interests) - 1
    max_after <- max (ts_after)
    if ( (max_after <= 0) && (max_before <= 0) )
      miic_error ("trajectories alignement", "impossible to align on variable ",
        params$var_interest, " around condition ", params$var_interest_condition,
        ": there is no trajectory with time steps before or after the condition.")
    #
    # If window_position around the positions of interest was requested,
    # look if possible : at least one position before and after
    # If not, change the window_position value
    #
    if (params$window_position == "around")
      {
      if (max_before <= 0)
        {
        miic_warning ("trajectories alignement",
          "the requested alignment on the variable ", params$var_interest,
          " was around condition ", params$var_interest_condition,
          " but there is no trajectory with time steps before the condition.",
          " The window position will be turned into 'start'.")
        params$window_position <- "start"
        }
      if (max_after <= 0)
        {
        miic_warning ("trajectories alignement",
          "the requested alignment on the variable ", params$var_interest,
          " was around condition ", params$var_interest_condition,
          " but there is no trajectory with time steps after the condition.",
          " The window position will be turned into 'end'.")
        params$window_position <- "end"
        }
      }
    #
    # If alignment is still around, there is a least one trajectory with
    # one time step before and after the position of interest
    #
    if (params$window_position == "around")
      {
      # As we need to set the window around the position of interest
      # (=> so at least 1 time step before and after), filter out the
      # trajectories with no time step before or after the position of interest
      #
      trajs_ko <- ( ((pos_interests - 1) <= 0) | (ts_after <= 0) )
      if ( any(trajs_ko) )
        {
        if ( sum(trajs_ko) == 1 )
          miic_warning ("trajectories alignement",
            "one trajectory can not be aligned on variable ", params$var_interest,
            " around condition ", params$var_interest_condition,
            " as there is no time step before or after the condition.",
            " This trajectory will be discarded.")
        else
          miic_warning ("trajectories alignement", sum (trajs_ko),
            " trajectories can not be aligned on variable ", params$var_interest,
            " around condition ", params$var_interest_condition,
            " as they have no time steps before or after the condition.",
            " These trajectories will be discarded.")
        list_traj <- list_traj[!trajs_ko]
        if ( length(list_traj) <= 0)
          miic_error ("trajectories alignement",
            "no trajectory available after alignment on variable ", params$var_interest,
            " around condition ", params$var_interest_condition, ".")
        pos_interests <- pos_interests[!trajs_ko]
        ts_lengths <- ts_lengths[!trajs_ko]
        ts_after <- ts_lengths - pos_interests
        }
      #
      # All trajectories have at least one time steps before and after the
      # position of interest, keep the min before and after to set the window
      #
      min_before <- min (pos_interests) - 1
      min_after <- min (ts_after)
      #
      # Align the trajectories with min_before, pos_interest, min_after
      #
      for ( ts_idx in 1:length(list_traj) )
        list_traj[[ts_idx]] <- list_traj[[ts_idx]][
          (pos_interests[[ts_idx]]-min_before):(pos_interests[[ts_idx]]+min_after), , F]
      aligned_pos_interest <- min_before + 1
      }
    #
    # If window_position is "start", keep only trajectories having time steps
    # after the condition, drop time steps before the positions of interest
    # and harmonize to the same length
    #
    else if (params$window_position == "start")
      {
      if (max_after <= 0)
        miic_error ("trajectories alignement",
          "impossible to align on variable ", params$var_interest,
          " from start using condition ", params$var_interest_condition,
          ": there is no trajectory with time steps after the condition.")

      trajs_ko <- (ts_after <= 0)
      if ( any(trajs_ko) )
        {
        if ( sum(trajs_ko) == 1 )
          miic_warning ("trajectories alignement",
            "one trajectory can not be aligned on variable ", params$var_interest,
            " from start using condition ", params$var_interest_condition,
            " as there is no time steps after the condition.",
            " This trajectory will be discarded.")
        else
          miic_warning ("trajectories alignement", sum (trajs_ko),
            " trajectories can not be aligned on variable ", params$var_interest,
            " from start using condition ", params$var_interest_condition,
            " as they have no time steps after the condition.",
            " These trajectories will be discarded.")
        list_traj <- list_traj[!trajs_ko]

        pos_interests <- pos_interests[!trajs_ko]
        ts_lengths <- ts_lengths[!trajs_ko]
        ts_after <- ts_lengths - pos_interests
        }
      min_after <- min (ts_after)
      for ( ts_idx in 1:length(list_traj) )
        list_traj[[ts_idx]] <- list_traj[[ts_idx]][
          pos_interests[[ts_idx]]:(pos_interests[[ts_idx]]+min_after), , F]
      aligned_pos_interest <- 1
      }
    #
    # If window_position is "end", keep only trajectories having time steps
    # before the condition, drop time steps after the positions of interest
    # and harmonize to the same length
    #
    else if (params$window_position == "end")
      {
      if (max_before <= 0)
        miic_error ("trajectories alignement",
          "impossible to align on variable ", params$var_interest,
          " from end using condition ", params$var_interest_condition,
          ": there is no trajectory with time steps before the condition.")

      trajs_ko <- ( (pos_interests - 1) <= 0)
      if ( any(trajs_ko) )
        {
        if ( sum(trajs_ko) == 1 )
          miic_warning ("trajectories alignement",
            "one trajectory can not be aligned on variable ", params$var_interest,
            " from end using condition ", params$var_interest_condition,
            " as there is no time steps before the condition.",
            " This trajectory will be discarded.")
        else
          miic_warning ("trajectories alignement", sum (trajs_ko),
            " trajectories can not be aligned on variable ", params$var_interest,
            " from end using condition ", params$var_interest_condition,
            " as they have no time steps before the condition.",
            " These trajectories will be discarded.")
        list_traj <- list_traj[!trajs_ko]
        pos_interests <- pos_interests[!trajs_ko]
        }

      min_before <- min (pos_interests) - 1
      for ( ts_idx in 1:length(list_traj) )
        list_traj[[ts_idx]] <- list_traj[[ts_idx]][
          (pos_interests[[ts_idx]]-min_before):pos_interests[[ts_idx]], , F]
      aligned_pos_interest <- nrow( list_traj[[1]] )
      }
    else
      {
      # Should never been raised
      miic_error ("trajectories alignement", params$window_position,
                  " is not a valid window position on the variable of interest.")
      }
    }
  return( list( "list_traj"=list_traj, "params"=params, "pos_interest"=aligned_pos_interest ) )
  }

#-------------------------------------------------------------------------------
# tmiic_non_stat_cover_dynamic
#-------------------------------------------------------------------------------
# In temporal non stationary mode, we don't try to estimate the temporal
# dynamic as we expect to have trend or seasonal effect in the data
# So the rule, if the number of layers and/or delta_t is not fixed by user,
# is to cover in the best possible way the length of the trajectory
# - list_traj: a list of data frames, each data frame is a trajectory
# - params: the list of MIIC params
# - pos_interest : position of the time step of interest in aligned trajectories
# - verbose: boolean, FALSE by default. Display some info it TRUE
#-------------------------------------------------------------------------------
tmiic_non_stat_cover_dynamic <- function (list_traj, state_order, params,
                                          pos_interest, verbose=F)
  {
  # NB : if one n_layers is NA in non stationary => all n_layers are NA and
  # same for delta_t : if delta_t is NA in non stationary => all delta_t are NA
  #
  if (  ( any( is.na(state_order$n_layers) ) )
     && (!all( is.na(state_order$n_layers) ) ) )
    miic_error( "temporal window assesment", "mix on undefined and defined number of layers.",
                " Please report this issue on github as it should have been detected earlier.")
  if (  ( any( is.na(state_order$delta_t) ) )
     && (!all( is.na(state_order$delta_t) ) ) )
    miic_error( "temporal window assesment", "mix on undefined and defined delta t.",
                " Please report this issue on github as it should have been detected earlier.")
  #
  # As trajectories are aligned, look only at first
  #
  ts_length <- nrow (list_traj[[1]])
  n_vars <- ncol (list_traj[[1]]) - 1
  nb_timesteps_before <- pos_interest - 1
  nb_timesteps_after <- nrow(list_traj[[1]]) - pos_interest
  #
  # Add n_layers before and after the position of interest
  #
  if ( all( is.na(state_order$n_layers) ) )
    {
    state_order$n_layers_before <- NA_integer_
    state_order$n_layers_after <- NA_integer_
    }
  else
    {
    if (any( grepl(",", state_order$n_layers, fixed=T) ) ) # => window "around"
      {
      state_order$n_layers_before <- unlist (lapply (state_order$n_layers,
        FUN = function(x) { as.integer (strsplit (x, ",", fixed=T)[[1]][[1]]) } ) )
      state_order$n_layers_after <- unlist (lapply (state_order$n_layers,
        FUN = function(x) { as.integer (strsplit (x, ",", fixed=T)[[1]][[2]]) } ) )
      }
    else if (params$window_position == "around") # n_layers is an int => centered window
      {
      state_order$n_layers_before <- (state_order$n_layers - 1) %/% 2
      state_order$n_layers_after <- state_order$n_layers - 1 - state_order$n_layers_before
      }
    else if (params$window_position == "start")
      {
      state_order$n_layers_before <- 0
      state_order$n_layers_after <- state_order$n_layers - 1
      }
    else # params$window_position == "end"
      {
      state_order$n_layers_before <- state_order$n_layers - 1
      state_order$n_layers_after <- 0
      }
    }
  #
  # Do check when both n_layers and delta_t are supplied first
  #
  if (  any( !is.na(state_order$n_layers) )
     && any( !is.na(state_order$delta_t) ) )
    {
    if (verbose)
      miic_msg ("Checking the layers and delta t to cover the trajectories ...")
    #
    # Check if these numbers of layers can fit before and after
    #
    for ( i in 1:nrow(state_order) )
      {
      # We ensure first that delta_t < time steps available
      #
      estim_delta_t <- state_order[i, "delta_t"]
      if (  (nb_timesteps_before > 0)
         && (state_order[i, "n_layers_before"] > 0)
         && (estim_delta_t > nb_timesteps_before) )
        estim_delta_t = nb_timesteps_before
      if (  (nb_timesteps_after > 0)
         && (state_order[i, "n_layers_after"] > 0)
         && (estim_delta_t > nb_timesteps_after) )
        estim_delta_t = nb_timesteps_after

      if (estim_delta_t != state_order[i, "delta_t"])
        {
        miic_warning( "temporal window assesment",
          "the number of time steps between layers for the variable ",
           state_order[i, "var_names"]," has been reduced from ",
           state_order[i, "delta_t"], " to ", estim_delta_t,
          " as the trajectories are too small.")
        state_order[i, "delta_t"] <- estim_delta_t
        }
      #
      # We can have at least one layer before and/or after
      #
      if ( grepl(",", state_order[i, "n_layers"], fixed=T) )
        {
        # We respect the layers defined by the user
        #
        estim_layers_before <- min (nb_timesteps_before %/% estim_delta_t,
                                    state_order[i, "n_layers_before"])
        estim_layers_after <- min (nb_timesteps_after %/% estim_delta_t,
                                   state_order[i, "n_layers_after"])
        updated_val <- paste0 (estim_layers_before, ",", estim_layers_after)
        }
      else
        {
        # We respect the total number of layers => the window can be not centered
        #
        estim_layers_before <- min (nb_timesteps_before %/% estim_delta_t,
                                    state_order[i, "n_layers"])
        estim_layers_after <- min (nb_timesteps_after %/% estim_delta_t,
                                   state_order[i, "n_layers"])
        while (  ( 1 + (estim_layers_before + estim_layers_after) * estim_delta_t >
                   1 + nb_timesteps_before + nb_timesteps_after)
              || (estim_layers_before + 1 + estim_layers_after >
                  state_order[i, "n_layers"]) )
          {
          if (estim_layers_before >= estim_layers_after)
            estim_layers_before <- estim_layers_before - 1
          else
            estim_layers_after <- estim_layers_after - 1
          }
        updated_val <- estim_layers_before + 1 + estim_layers_after
        }

      if ( (updated_val != state_order[i, "n_layers"])
         || (estim_layers_before != state_order[i, "n_layers_before"])
         || (estim_layers_after != state_order[i, "n_layers_after"]) )
        {
        if (updated_val != state_order[i, "n_layers"])
          miic_warning( "temporal window assesment",
            "the number of layers for the variable ",
             state_order[i, "var_names"]," has been reduced from ",
             state_order[i, "n_layers"], " to ", updated_val,
            " as the trajectories are too small.")
        state_order[i, c("n_layers", "n_layers_before", "n_layers_after")] <-
          c(updated_val, estim_layers_before, estim_layers_after)
        }
      }
    }
  else if ( ! any (is.na (state_order$n_layers) ) )
    {
    # Case where n_layers known and delta_t not defined
    #
    miic_msg ("Estimating delta t to cover the trajectories ...")
    for ( i in 1:nrow(state_order) )
      {
      # Even if n_layers are supplied by user, they can not exceed the available
      # time steps
      #
      estim_layers_before <- min (nb_timesteps_before,
                                  state_order[i, "n_layers_before"])
      estim_layers_after <- min (nb_timesteps_after,
                                  state_order[i, "n_layers_after"])
      #
      # Estimate delta_t
      #
      estim_delta_t_before <- NA_integer_
      if (estim_layers_before > 0)
        estim_delta_t_before <- max( 1, nb_timesteps_before %/% estim_layers_before)
      estim_delta_t_after <- NA_integer_
      if (estim_layers_after > 0)
        estim_delta_t_after <- max( 1, nb_timesteps_after %/% estim_layers_after)
      estim_delta_t <- max (estim_delta_t_before, estim_delta_t_after, na.rm=T)
      state_order[i, "delta_t"] <- estim_delta_t

      if (  (estim_layers_before != state_order[i, "n_layers_before"])
         || (estim_layers_after != state_order[i, "n_layers_after"]) )
        {
        if (grepl( ",", state_order[i, "n_layers"], fixed=T) )
          updated_val <- paste0 (estim_layers_before, ",", estim_layers_after)
        else
          updated_val <- estim_layers_before + 1 + estim_layers_after

        miic_warning( "temporal window assesment",
          "the number of layers for the variable ",
           state_order[i, "var_names"]," has been reduced from ",
           state_order[i, "n_layers"], " to ", updated_val,
          " as the trajectories are too small.")
        state_order[i, c("n_layers", "n_layers_before", "n_layers_after")] <-
          c(updated_val, estim_layers_before, estim_layers_after)
        }
      }

    if (verbose)
      {
      if (length( unique(state_order$delta_t) ) <= 1)
        miic_msg( "Estimated time steps between layers: ", unique(state_order$delta_t), "." )
      else
        {
        miic_msg( "Estimated time steps between layers:")
        for ( i in 1:nrow(state_order) )
          {
          miic_msg( "- ", state_order[i, "var_names"], "\t",
            state_order[i, "delta_t"], " (for ",
            state_order[i, "n_layers"], " layers)")
          }
        }
      }
    }
  else if ( ! any (is.na (state_order$delta_t) ) )
    {
    # Cases where delta_t known and n_layers not defined
    #
    miic_msg ("Estimating number of layers to cover the trajectories ...")
    for ( i in 1:nrow(state_order) )
      {
      estim_delta_t <- state_order[i, "delta_t"]
      #
      # Even if delta_t is supplied by user, it can not exceed the available
      # time steps
      #
      max_delta_t_before <- NA_integer_
      if (nb_timesteps_before > 0)
        max_delta_t_before <- nb_timesteps_before
      max_delta_t_after <- NA_integer_
      if (nb_timesteps_after > 0)
        max_delta_t_after <- nb_timesteps_after
      max_delta_t <- min (max_delta_t_before, max_delta_t_after, na.rm=TRUE)

      if (estim_delta_t > max_delta_t)
        estim_delta_t <- max_delta_t

      if (estim_delta_t != state_order[i, "delta_t"])
        {
        miic_warning( "temporal window assesment",
          "the number of time steps between layers for the variable ",
           state_order[i, "var_names"]," has been reduced from ",
           state_order[i, "delta_t"], " to ", estim_delta_t,
          " as the trajectories are too small.")
        state_order[i, "delta_t"] <- estim_delta_t
        }
      #
      # Delta t is compatible with the trajectories
      # (=> we can have at least one layer before or/and after)
      #
      # Determine max number of layers
      #
      estim_layers_before <- nb_timesteps_before %/% estim_delta_t
      estim_layers_after <- nb_timesteps_after %/% estim_delta_t

      max_layers <- params$max_nodes %/% n_vars
      while (estim_layers_before + 1 + estim_layers_after > max_layers)
        {
        if (estim_layers_before >= estim_layers_after)
          estim_layers_before <- estim_layers_before - 1
        else
          estim_layers_after <- estim_layers_after - 1
        }

      if (  (params$window_position  == "around")
         && (estim_layers_before != estim_layers_after) )
        estim_layers <- paste0 (estim_layers_before, ",", estim_layers_after)
      else
        estim_layers <- estim_layers_before + 1 + estim_layers_after
      state_order[i, c("n_layers", "n_layers_before", "n_layers_after")] <-
        c(estim_layers, estim_layers_before, estim_layers_after)
      }

    if (verbose)
      {
      if (length( unique(state_order$n_layers) ) <= 1)
        miic_msg( "Estimated number of layers: ", unique(state_order$n_layers), "." )
      else
        {
        miic_msg( "Estimated number of layers:")
        for ( i in 1:nrow(state_order) )
          {
          miic_msg( "- ", state_order[i, "var_names"], "\t",
            state_order[i, "n_layers"], " (with ",
            state_order[i, "delta_t"], " time steps between layers)")
          }
        }
      }
    }
  else
    {
    # Last case where both n_layers and delta_t are not supplied by the user
    #
    miic_msg ("Estimating number of layers and time between layers to cover the trajectories ...")
    #
    # The number of layers is first computed from the max nodes in the network
    #
    estim_layers <- params$max_nodes %/% n_vars
    if (estim_layers < 1)
      {
      if (params$window_position == "around")
        estim_layers <- 3
      else
        estim_layers <- 2
      miic_warning( "temporal window assesment",
        "impossible to respect the maximum nodes of ", params$max_nodes,
        ". The number of layers has been set to ", estim_layers,
        " (", estim_layers * n_vars, " nodes).")
      }
    #
    # Reduce this number of layers if too much for the time steps available
    #
    nb_timestep_tot <- nb_timesteps_before + 1 + nb_timesteps_after
    if (estim_layers > nb_timestep_tot)
      estim_layers <- nb_timestep_tot
    estim_delta_t <- max (1, (nb_timestep_tot - 1) %/% (estim_layers - 1) )
    #
    # We have base n_layers/delta_t, look if we have to tune them
    # in regard of the position of interest:
    # If the pos_interest is not start or the end or the trajectories
    # and is not a multiple of delta t, it can be needed to reduce the
    # n_layers or delta_t:
    #
    # i.e.: n_layers = 3 on an aligned trajectory of 7 time steps
    # base delta_t : (7 - 1) // (3 - 1) => 3 time steps
    # if pos_interest = 1 : t1 t2 t3 t4 t5 t6 t7
    #                       L  x  x  L  x  x  L => OK, fit nicely
    # if pos_interest = 4 : t1 t2 t3 t4 t5 t6 t7
    #                       L  x  x  L  x  x  L => OK, no need to reduce delta_t
    # if pos_interest = 5 : t1 t2 t3 t4 t5 t6 t7
    #                       x  L  x  x  L  x  x => KO, only 2 layers
    #                       L  x  L  x  L  x  L => KO, 4 layers if delta_t reduced
    #                                           => Chose not reduced delta_t
    #                                              (minimal number of nodes)
    # if pos_interest = 6 : t1 t2 t3 t4 t5 t6 t7
    #                       x  x  L  x  x  L  x => KO, only 2 layers
    #                       x  L  x  L  x  L  x => OK if delta_t reduced
    # if pos_interest = 7 : t1 t2 t3 t4 t5 t6 t7
    #                       L  x  x  L  x  x  L => OK, fit nicely
    #
    estim_layers_before <- nb_timesteps_before %/% estim_delta_t
    estim_layers_after <- nb_timesteps_after %/% estim_delta_t
    if (estim_layers_before + 1 + estim_layers_after != estim_layers)
      {
      other_delta_t <- estim_delta_t - 1
      other_layers_before <- nb_timesteps_before %/% other_delta_t
      other_layers_after <- nb_timesteps_after %/% other_delta_t
      if (other_layers_before + 1 + other_layers_after == estim_layers)
        {
        estim_delta_t <- other_delta_t
        estim_layers_before <- other_layers_before
        estim_layers_after <- other_layers_after
        }
      }
    if (  (estim_layers_before > 0)
       && (estim_layers_after > 0)
       && (estim_layers_before != estim_layers_after) )
      estim_layers <- paste0 (estim_layers_before, ",", estim_layers_after)
    else
      estim_layers <- estim_layers_before + 1 + estim_layers_after
    #
    # Estimation completed
    #
    state_order$n_layers <- estim_layers
    state_order$delta_t <- estim_delta_t
    state_order$n_layers_before <- estim_layers_before
    state_order$n_layers_after <- estim_layers_after
    if (verbose)
      miic_msg( "Estimated parameters to cover the trajectories: ", estim_layers,
                " layers with ", estim_delta_t, " time steps between layers.")
    }

  return (state_order)
  }

#-------------------------------------------------------------------------------
# tmiic_non_stat_tune_inputs
#-------------------------------------------------------------------------------
# In temporal non stationary mode, we adjust parameter, input data and state
# order after the temporal window check or estimation
# - list_inputs: a list of input data frames, each data frame is a trajectory
# - pos_interest : position of the time step of interest in aligned trajectories
#-------------------------------------------------------------------------------
tmiic_non_stat_tune_inputs <- function (list_traj, state_order, params, pos_interest)
  {
  # Tune data alignment (is not necessary to perform the data lagging
  # but return the true minimal input data with these parameters)
  #
  nb_ts_before <- max (state_order$n_layers_before * state_order$delta_t)
  nb_ts_after <- max (state_order$n_layers_after * state_order$delta_t)
  for (i in 1:length(list_traj) )
    list_traj[[i]] <- list_traj[[i]][
      (pos_interest - nb_ts_before):(pos_interest + nb_ts_after), ,drop=FALSE]
  pos_interest <- as.integer (nb_ts_before + 1)
  #
  # Tune state order
  #
  if (  (params$window_position  != "start")
     && (all (state_order$n_layers_before == 0)) )
    {
      miic_warning( "tuning inputs",
        "the window position has been set to 'start'",
        " as there is no time step used before the variable of interest.")
      params$window_position  <- "start"
      state_order$n_layers <- state_order$n_layers_after + 1
      }
  if (  (params$window_position  != "end")
     && (all (state_order$n_layers_after == 0)) )
    {
    miic_warning( "tuning inputs",
      "the window position has been set to 'end'",
      " as there is no time step used after the variable of interest.")
    params$window_position  <- "end"
    state_order$n_layers <- state_order$n_layers_before + 1
    }

  list_ret <- list( "list_traj" = list_traj, "state_order"= state_order,
                    "params" = params, "pos_interest" = pos_interest )
  return (list_ret)
  }

#-------------------------------------------------------------------------------
# tmiic_non_stat_lag_state_order
#-------------------------------------------------------------------------------
# Modify the state order into a lagged version in non stationary mode:
# the lagged variables are completed and/or repeated with t-X, t, t+X to
# match the n_layers before/after and delta_t parameters
#
# - state_order: a dataframe, the state order returned by
#   tmiic_non_stat_tune_inputs, with n_layers before/after and delta_t defined
# Returns: a dataframe: the lagged state_order
#-------------------------------------------------------------------------------
tmiic_non_stat_lag_state_order <- function (state_order)
  {
  n_vars <- nrow (state_order)
  state_order$lag <- -9999
  state_order$var_idx_data <- -9999
  #
  # Put lag0 and not lagged variable first
  #
  state_lagged <- state_order
  for (var_idx in 1:n_vars)
    {
    state_lagged [var_idx, "var_names"] <- paste0 (state_order [var_idx, "var_names"], "_lag0")
    state_lagged [var_idx, "lag"] <- 0
    state_lagged [var_idx, "var_idx_data"] <- var_idx
    }
  #
  # Duplicate rows for lagged variables
  #
  state_lagged_nrows <- nrow (state_lagged)
  for (var_idx in 1:n_vars)
    {
    n_layers_before <- state_lagged[var_idx, "n_layers_before"]
    while (n_layers_before > 0)
      {
      state_lagged_nrows <- state_lagged_nrows + 1
      state_lagged [state_lagged_nrows,] <- state_order [var_idx,]
      lag <- n_layers_before * state_order[var_idx, "delta_t"]
      state_lagged [state_lagged_nrows, "var_names"] <- paste0 (
        state_order [var_idx, "var_names"], "_lag-", lag)
      state_lagged [state_lagged_nrows, "lag"] <- -lag
      state_lagged [state_lagged_nrows, "var_idx_data"] <- var_idx
      n_layers_before <- n_layers_before - 1
      }
    n_layers_after <- state_lagged[var_idx, "n_layers_after"]
    while (n_layers_after > 0)
      {
      state_lagged_nrows <- state_lagged_nrows + 1
      state_lagged [state_lagged_nrows,] <- state_order [var_idx,]
      lag <- n_layers_after * state_order[var_idx, "delta_t"]
      state_lagged [state_lagged_nrows, "var_names"] <- paste0 (
        state_order [var_idx, "var_names"], "_lag+", lag)
      state_lagged [state_lagged_nrows, "lag"] <- lag
      state_lagged [state_lagged_nrows, "var_idx_data"] <- var_idx
      n_layers_after <- n_layers_after - 1
      }
    }
  #
  # Apply a reordering (not truly important but easier to read for an human,
  # this order will be used also when lagging input data)
  #
  order_rows <- order(state_lagged$lag, state_lagged$var_names)
  state_lagged <- state_lagged[order_rows,]
  rownames(state_lagged) <- NULL
  return (state_lagged)
  }

#-------------------------------------------------------------------------------
# tmiic_non_stat_lag_other_df
#-------------------------------------------------------------------------------
# Modify the complementary df int a lagged version: the 4 columns dataframes are
# transformed into a 2 columns one, in which variables are transformed into
# their lagged representation. i.e:
# - var1, -2, var2, 1 becomes var1_lag-2 - var2_lag+1
# inputs:
# - df: the dataframe to transform in its lagged version
# Returns: a dataframe: the lagged dataframe
#-------------------------------------------------------------------------------
tmiic_non_stat_lag_other_df <- function (df)
  {
  if ( is.null(df) )
    return( df )
  if ( nrow(df) > 0 )
    {
    for (i in 1:nrow (df))
      {
      if (df [i, 2] <= 0)
        df[i, 1] <- paste0 (df [i, 1], "_lag", df [i, 2])
      else
        df[i, 1] <- paste0 (df [i, 1], "_lag+", df [i, 2])
      if (df [i, 4] <= 0)
        df[i, 3] <- paste0 (df [i, 3], "_lag", df [i, 4])
      else
        df[i, 3] <- paste0 (df [i, 3], "_lag+", df [i, 4])
      }
    }
  df <- df[,c(1,3)]
  return (df)
  }

#-------------------------------------------------------------------------------
# tmiic_non_stat_lag_input_data
#-------------------------------------------------------------------------------
# Reorganizes the inputs in a format usable by miic: input data are lagged
# using the history to create lagged variables
# The function slices the input data according to the information supplied in
# the state_order n_layers and delta_t.
#
# The number of variables is increased and renamed on n_layers
# layers by delta_t. steps.
# i.e. with n_layers_before=2, n_layers_after=1 and delta_t.=3 :
# var1, var2 => var1_lag-6, var2_lag-6, var1_lag-3, var2_lag-3,
#               var1_lag0, var2_lag0, var1_lag+3, var2_lag+3
#
# Note that the lagging can be different for each input variable
# if different values of n_layers or delta_t are supplied
#
# inputs:
# - list_traj: the list of timeseries
# - state_order: a dataframe, the lagged state order returned by
#   tmiic_non_stat_lag_state_order
# - keep_max_data: boolean flag, optional, FALSE by default
#   When FALSE, the rows containing NA introduced by the lagging process
#   are deleted, otherwise when TRUE, the rows are kept
#-------------------------------------------------------------------------------
tmiic_non_stat_lag_input_data <- function (list_traj, lagged_order, pos_interest)
  {
  n_vars <- nrow (lagged_order)
  na_count <- 0
  list_ret <- list()
  for ( ts_idx in 1:length(list_traj) )
    {
    df <- list_traj[[ts_idx]]
    #
    # Lag the df into a df with one row and more columns
    #
    list_tmp <- list()
    var_idx <- 1
    for ( var_idx in 1:n_vars )
      {
      source_idx <- lagged_order[var_idx, "var_idx_data"]
      lag <- lagged_order[var_idx, "lag"]
      list_tmp[[var_idx]] <- df[(pos_interest+lag),(source_idx+1),F]
      }
    df <- as.data.frame (do.call (cbind, list_tmp) )
    if (nrow (df) != 1)
      stop("aze")
    colnames(df) <- lagged_order$var_names
    #
    # Check if the row has only NAs
    #
    rows_only_na <- ( rowSums (is.na (df)) == ncol (df) )
    df <- df [!rows_only_na, , F]
    na_count <- na_count + sum (rows_only_na)
    list_ret[[ts_idx]] <- df
    }
  if (na_count > 0)
    miic_warning ("data lagging", "the lagged data contains ", sum(na_count),
             " row(s) with only NAs. These row(s) have been removed.")
  return (list_ret)
  }

#-------------------------------------------------------------------------------
# tmiic_non_stat_to_ml
#-------------------------------------------------------------------------------
# Convert non stationary parameters into layers.
# With the introduction of the multi-layered version of MIIC, temporal non
# stationary mode can be considered as a specific setup of multi-layers network
#-------------------------------------------------------------------------------
tmiic_non_stat_to_ml <- function (list_in)
  {
  # Associate each node to a layer corresponding to its lag
  #
  node_layers = rep ( NA_integer_, nrow (list_in$state_order) )
  #
  # Define layers so that contributors of a layers can only be past or present
  # (as we can not use the future to condition an edge)
  # and set  pre-orientation according to time: edges will be pre-oriented
  # from a past layer toward the more recent layer
  #
  df_ml <- data.frame ("contributors" = character(),
                       "preoriented" = character(),
                       stringsAsFactors = F)

  uniq_lags = sort (unique (list_in$state_order$lag) )
  for (layer_idx in 1:length(uniq_lags))
    {
    one_lag <- uniq_lags[[layer_idx]]
    node_layers[list_in$state_order$lag == one_lag ] <- layer_idx
    contribs_str <- paste (1:layer_idx, collapse=";")
    if (layer_idx <= 1)
      preoriented_str <- as.character (FALSE)
    else
      {
      preoriented_bool <- rep (TRUE, layer_idx)
      preoriented_bool[[layer_idx]] <- FALSE
      preoriented_str <- paste (preoriented_bool, collapse=";")
      }
    df_ml[layer_idx,] <- list (contribs_str, preoriented_str)
    }
  # df_ml[3,] <- list ("3", "F")
  list_in$state_order$layers = node_layers
  list_in$layers = df_ml
  return (list_in)
  }
