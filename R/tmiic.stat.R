#*******************************************************************************
# Filename   : tmiic.stat.R                   Creation date: 31 jan 2024
#
# Description: Utility functions for temporal MIIC in stationary mode
#
# Author     : Franck SIMON
#
# Changes history:
# - 31 jan 2024 : initial version
#*******************************************************************************

#===============================================================================
# FUNCTIONS
#===============================================================================
# tmiic_stat_ajust_window_for_nb_samples
#-----------------------------------------------------------------------------
# Reduce the window size (n_layers or delta_t) if the foreseen n_layers and
# delta_t would lead to too few samples after data lagging
# params:
# - list_traj: a list of data frame, each item representing a trajectory.
# - n_layers: a list, the n_layers in the state_order column
# - delta_t: a list, the delta_t in the state_order column
# - reduced_param: a string, can be "n_layers" or "delta_t". Indicates with
#   parameter will be reduced if the number of samples is too small
# - verbose: boolean, if TRUE, display a message if the window size is reduced
# returns:
# - a list: the n_layers or delta_t, depending of the reduced_param value.
#   The value are possibly decreased to reduce the window size
#-----------------------------------------------------------------------------
tmiic_stat_ajust_window_for_nb_samples <- function (list_traj, n_layers, delta_t,
                                                    reduced_param, verbose)
  {
  tau_per_var <- (n_layers - 1) * delta_t
  tau_max <- max (tau_per_var)
  ts_lengths <- unlist ( lapply (list_traj, nrow) )
  tot_ts <- sum (ts_lengths)
  nb_samples <- sum ( unlist (lapply (ts_lengths, FUN=function (x) {
                                            max (0, x - tau_max) } ) ) )
  target <- min (1000, tot_ts / 10)

  if (nb_samples < target)
    {
    # Look for the best value to reach the target recursivly
    # At each iteration, we keep the half part where the best value is
    # until we can not divide in half further
    #
    recurs_eval <- function (target, tau_low, tau_high, ts_lengths)
      {
      if (tau_high - tau_low <= 1)
        return (tau_low)
      tau <- round ( (tau_low + tau_high) / 2, 0)
      nb_samples <- sum ( unlist (lapply (ts_lengths, FUN=function (x) {
                                            max (0, x - tau) } ) ) )
      if (nb_samples >= target)
        tau_ret <- recurs_eval (target, tau, tau_high, ts_lengths)
      else
        tau_ret <- recurs_eval (target, tau_low, tau, ts_lengths)
      return (tau_ret)
      }
    tau_red <- recurs_eval (target, 1, tau_max, ts_lengths)
    #
    # Max time steps back in time found, (try to) reduce n_layers or delta_t
    #
    if (reduced_param == "n_layers")
      n_layers[tau_per_var > tau_red] <- max ( 2,
        floor (tau_red / delta_t[tau_per_var > tau_red]) + 1)
    else # reduce delta_t
      delta_t[tau_per_var > tau_red] <- max ( 1,
        floor (tau_red / (n_layers[tau_per_var > tau_red] - 1)) )
    #
    # Check the effect of reduction and feed back to user
    #
    if (reduced_param == "n_layers")
      fixed_param <- "delta_t"
    else
      fixed_param <- "n_layers"
    tau_max_red <- max ( (n_layers - 1) * delta_t )
    nb_samples_red <- sum ( unlist (lapply (ts_lengths, FUN=function (x) {
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
# tmiic_stat_estimate_dynamic
#-------------------------------------------------------------------------------
# Estimate tau (the number of total time steps back to cover the dynamic,
# the number of layers and delta t parameters from the data
# - list_traj: list of data frame, each item representing a trajectory.
#   Each data frame must contain the time step information in the 1st column
#   and the variables in the other columns.
# - state_order: the state_order data frame. This state_order is expected
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
tmiic_stat_estimate_dynamic <- function (list_traj, state_order, max_nodes=50,
                                         verbose_level=1)
  {
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
  if (verbose_level >= 2)
    miic_msg ("Estimating the temporal dynamic...")
  n_ts <- length (list_traj)
  n_vars_tot <- ncol (list_traj[[1]]) - 1
  n_vars_ctx <- sum (state_order$is_contextual)
  n_vars_lag <- n_vars_tot - n_vars_ctx
  #
  # Remove time step, contextual and discrete variables
  #
  if ( ! any ( (state_order$is_contextual == 0) & (state_order$var_type == 1) ) )
    miic_error ("dynamic estimation", "no variable to estimate the temporal dynamic",
                " (all variables are discrete or contextual).",
                " Consider specifying the n_layers and delta_t parameters.")
  for (ts_idx in 1:n_ts)
    {
    if (nrow (list_traj[[ts_idx]]) == 1)
      miic_warning ("dynamic estimation", "trajectory ", ts_idx,
                    " with only 1 time step is ignored for dynamic estimation.")
    list_traj[[ts_idx]] <- list_traj[[ts_idx]][, c(F, ( (state_order$is_contextual == 0)
                                                  & (state_order$var_type == 1) ) ), F]
    }
  #
  # Compute mean alpha per variable
  #
  n_vars <- ncol (list_traj[[1]])
  var_names <- colnames(list_traj[[1]])
  length_to_test <- min (unlist (lapply (list_traj, FUN=function (x) {
    ifelse ( nrow(x) <= 1, NA, nrow(x) ) } ) ), na.rm=T)
  alphas_per_var <- rep (NA, n_vars)
  taus_per_var <- rep (NA, n_vars)
  var_idx <- 1
  for (var_idx in 1:n_vars)
    {
    alphas_per_ts <- rep (NA, n_ts)
    ts_idx <- 1
    for (ts_idx in 1:n_ts)
      {
      if (nrow (list_traj[[ts_idx]]) == 1)
        next
      acf_res <- acf (list_traj[[ts_idx]][,(var_idx)], na.action=na.pass,
                      lag.max=length_to_test-1, plot=F)
      if ( all (is.na(acf_res$acf) ) )
        next
      acf_vanish <- which (acf_res$acf[,1,1] < 0.05)
      if ( length (acf_vanish) == 0 )
        acf_vanish <- length_to_test
      lag_vanish <- acf_res$lag[min (acf_vanish), 1, 1]
      lag_4_alpha <- max ( 1, round (lag_vanish / 2) )
      alphas_per_ts[[ts_idx]] <- acf_res$acf[lag_4_alpha+1,1,1] ^ (1/lag_4_alpha)
      }
    alphas_per_var[[var_idx]] <- mean (alphas_per_ts, na.rm=T)
    taus_per_var[[var_idx]] <- round ( (1+alphas_per_var[[var_idx]])
                                    / (1-alphas_per_var[[var_idx]]) )
    }
  if (verbose_level >= 2)
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

  tau_min  <- max ( 1, min (taus_per_var, na.rm=T) )
  tau_mean <- max ( 1, round (mean (taus_per_var, na.rm=T), 0) )
  tau_max  <- min ( length_to_test, max  (taus_per_var, na.rm=T) )
  tau_max_kept <- min (length_to_test, tau_max, tau_mean * 2)
  tau <- tau_max_kept
  if (verbose_level >= 1)
    miic_msg ("Automatic estimation of parameters:\n",
      "- Relaxation times goes from ", tau_min, " to ", tau_max,
      " with a mean of ", tau_mean, ", tau max considered = ", tau_max_kept)
  #
  # We know tau : the average maximum time steps back in time to use for the
  # temporal discovery. Now estimate the number of layers 'n_layers'
  # and/or number of time steps between two layers 'delta_t'
  #
  if (  all (!is.na (state_order$n_layers)) ) # n_layers known => NAs in delta_t
    {
    state_order$delta_t[is.na(state_order$delta_t)] <- max ( 1,
        ceiling (tau / (state_order$n_layers[is.na(state_order$delta_t)] - 1)) )

    state_order$delta_t <- tmiic_stat_ajust_window_for_nb_samples (list_traj,
      state_order$n_layers, state_order$delta_t, reduced_param="delta_t",
      verbose=(verbose_level >= 1) )

    uniq_n_layers <- unique (state_order$n_layers[ (state_order$is_contextual == 0)
                                                & (state_order$var_type == 1)] )
    uniq_delta_t <- unique (state_order$delta_t[ (state_order$is_contextual == 0)
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
    n_layers_max <- max (2, floor ( (max_nodes - n_vars_ctx) / n_vars_lag ) )
    #
    # The final number of layers will (tau / delta_t) + 1 unless if greater
    # than the max number of layers
    #
    state_order$n_layers[is.na(state_order$n_layers)] <- min (n_layers_max,
      ceiling (tau / state_order$delta_t[is.na(state_order$n_layers)]) + 1)

    state_order$n_layers <- tmiic_stat_ajust_window_for_nb_samples (list_traj,
      state_order$n_layers, state_order$delta_t, reduced_param="n_layers",
      verbose=(verbose_level >= 1) )

    uniq_n_layers <- unique (state_order$n_layers[ (state_order$is_contextual == 0)
                                                 & (state_order$var_type == 1)] )
    uniq_delta_t <- unique (state_order$delta_t[ (state_order$is_contextual == 0)
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
    delta_t <- 1
    if ( (tau + 1) * n_vars_lag + n_vars_ctx <= max_nodes)
      {
      # If when using delta_t = 1, the n_layers (= tau + 1) does not lead to
      # a graph with a total number of nodes > max => OK, nothing more to do
      #
      n_layers <- tau + 1
      }
    else
      {
      # We need reduce the number of layers to respect the maximum nodes number
      # and increase de delta t to still cover all the dynamic tau.
      # => Compute the max number of layers and deduce the delta t
      #
      n_layers <- max (2, floor ( (max_nodes - n_vars_ctx) / n_vars_lag ) )
      if (n_layers > 2)
        {
        delta_t <- max (1, ceiling ( tau / (n_layers-1)  ) )
        tau <- (n_layers - 1) * delta_t
        }
      else
        delta_t <- tau
      }

    state_order$n_layers[is.na(state_order$n_layers)] <- n_layers
    state_order$delta_t[is.na(state_order$delta_t)] <- delta_t

    state_order$delta_t <- tmiic_stat_ajust_window_for_nb_samples (list_traj,
      state_order$n_layers, state_order$delta_t, reduced_param="delta_t",
      verbose=(verbose_level >= 1) )
    delta_t <- unique (state_order$delta_t[ (state_order$var_type == 1)
                                          & (state_order$is_contextual == 0) ])

    if (verbose_level >= 1)
      miic_msg ("- For a final graph with a target of ", max_nodes,
        " nodes having ", n_vars_lag, " lagged variables",
        ifelse (n_vars_ctx > 0, paste0 ("\n  and ", n_vars_ctx, " contextual variables"), ""),
        ":\n  ", n_layers,  " layers spaced by ", delta_t, " time steps",
        ", dynamic covered goes over t, t-", delta_t,
        ifelse (n_layers > 3, ", ...", ""),
        ifelse (n_layers > 2, paste0 (", t-", tau), "") )
    }

  return (state_order)
  }

#-------------------------------------------------------------------------------
# tmiic_stat_lag_state_order
#-------------------------------------------------------------------------------
# Modify the state order into a lagged version: the lagged variables are
# completed and/or repeated with lagX to match the lagged temporal graph.
# inputs:
# - state_order: a dataframe, the state order returned by
#   tmiic_check_state_order_part2
# Returns: a dataframe: the lagged state_order
#-------------------------------------------------------------------------------
tmiic_stat_lag_state_order <- function (state_order)
  {
  n_vars <- nrow (state_order)
  state_order$lag <- -1
  state_order$var_idx_data <- -1
  #
  # Put lag0 and not lagged variable first
  #
  state_lagged <- state_order
  for (var_idx in 1:n_vars)
    {
    if (state_lagged [var_idx, "is_contextual"] == 0)
      {
      state_lagged [var_idx, "var_names"] = paste0 (state_order [var_idx, "var_names"], "_lag0")
      state_lagged [var_idx, "lag"] = 0
      state_lagged [var_idx, "var_idx_data"] = var_idx
      }
    else
      {
      state_lagged [var_idx, "lag"] = 0
      state_lagged [var_idx, "var_idx_data"] = var_idx
      }
    }
  #
  # Duplicate rows for lagged variables
  #
  state_lagged_nrows <- nrow (state_lagged)
  n_layers_back_max <- max ( (state_lagged$n_layers - 1) )
  n_layers_back_idx = 1
  for (n_layers_back_idx in 1:n_layers_back_max)
    {
    var_idx = 1
    for (var_idx in 1:n_vars)
      {
      n_layers_back_of_var <- state_lagged[var_idx, "n_layers"]  - 1
      if (n_layers_back_idx <= n_layers_back_of_var)
        {
        state_lagged_nrows <- state_lagged_nrows + 1
        state_lagged [state_lagged_nrows,] <- state_order [var_idx,]
        lag <- n_layers_back_idx * state_order[var_idx, "delta_t"]
        state_lagged [state_lagged_nrows, "var_names"] <- paste0 (
          state_order [var_idx, "var_names"], "_lag", lag)
        state_lagged [state_lagged_nrows, "lag"] <- lag
        state_lagged [state_lagged_nrows, "var_idx_data"] <- var_idx
        }
      }
    }
  return (state_lagged)
  }

#-------------------------------------------------------------------------------
# tmiic_stat_lag_other_df
#-------------------------------------------------------------------------------
# Modify the complementary df int a lagged version: the 3 column dataframes are
# transformed into a 2 columns one, in which variables are transformed into
# their lagged representation. i.e:
# - normal_var1 - normal_var2 - 1 becomes normal_var1_lag1 - normal_var2_lag0
# - ctx_var1 - normal_var2 - NA becomes ctx_var1 - normal_var2_lag0
# inputs:
# - state_order: a dataframe, the state order returned by
#   tmiic_check_state_order_part2
# - df: the dataframe to transform in its lagged version
# Returns: a dataframe: the lagged dataframe
#-------------------------------------------------------------------------------
tmiic_stat_lag_other_df <- function (state_order, df)
  {
  if (is.null (df))
    return (df)

  if (nrow (df) > 0)
    {
    for (i in 1:nrow (df))
      {
      orig_node_idx <- which (state_order$var_names == df[i, 1])
      if (state_order[orig_node_idx, "is_contextual"] == 0)
        df[i, 1] = paste0 (df [i, 1], "_lag", df [i, 3])
      df[i, 2] = paste0 (df [i, 2], "_lag0")
      }
    }
  df <- df[,c(1,2)]
  return (df)
  }

#-----------------------------------------------------------------------------
# tmiic_stat_lag_input_data
#-------------------------------------------------------------------------------
# Reorganizes the inputs in a format usable by miic: input data are lagged
# using the history to create lagged variables
# The function slices the input data according to the information supplied in
# the state_order n_layers and delta_t.
#
# The number of variables is increased and renamed on n_layers
# layers by delta_t. steps.
# i.e. with n_layers=3 and delta_t=3 : var1, var2 =>
# var1_lag0, var2_lag0, var1_lag3, var2_lag3, var1_lag6, var2_lag6.
#
# Every time step (until number of time steps - (n_layers  - 1) * delta_t.)
# is converted into a sample in the lagged data.
#
# Example with n_layers=3 and delta_t.=3:
#
# Timestep Var & value    Var & value  => Sample  Var & value   Var & value
#   t-6    Var1_val(t-6) Var2_val(t-6) =>   i    Var1_lag6_val Var2_lag6_val
#   t-3    Var1_val(t-3) Var2_val(t-3) =>   i    Var1_lag3_val Var2_lag3_val
#    t      Var1_val(t)   Var2_val(t)  =>   i    Var1_lag0_val Var2_lag0_val
#
#   t-7    Var1_val(t-7) Var2_val(t-7) =>   i'   Var1_lag6_val Var2_lag6_val
#   t-4    Var1_val(t-4) Var2_val(t-4) =>   i'   Var1_lag3_val Var2_lag3_val
#   t-1    Var1_val(t-1) Var2_val(t-1) =>   i'   Var1_lag0_val Var2_lag0_val
#
#   t-8    Var1_val(t-8) Var2_val(t-8) =>   i"   Var1_lag6_val Var2_lag6_val
#   t-5    Var1_val(t-5) Var2_val(t-5) =>   i"   Var1_lag3_val Var2_lag3_val
#   t-2    Var1_val(t-2) Var2_val(t-2) =>   i"   Var1_lag0_val Var2_lag0_val
#
#   ...    ............. ............. => ...... ............. ............
#
# until number of time steps - (n_layers - 1) * delta_t is reached.
# The same process is applied to all input time series.
#
# Note that the lagging can be different for each input variable
# if different values of n_layers or delta_t are supplied and some
# variables can be not lagged at all like contextual ones.
#
# inputs:
# - list_traj: the list of time series
# - state_order: a dataframe, the lagged state order returned by
#   tmiic_stat_lag_state_order
# - keep_max_data: boolean flag, optional, FALSE by default
#   When FALSE, the rows containing NA introduced by the lagging process
#   are deleted, otherwise when TRUE, the rows are kept
#
# TODO: remove lagged variables with columns full of NAs ?
# If the lagging results in a lagged variable full of NAs, this variable
# will not be connected in the reconstructed graph. So it is useless to send it
# to the C++ reconstruct part.
# Caution however: if some lagged variables are removed, it will likely
# jeopardize the lagged graph plotting
#-------------------------------------------------------------------------------
tmiic_stat_lag_input_data <- function (list_traj, state_order, keep_max_data=FALSE)
  {
  tau_max = max(state_order$lag)
  na_count = 0
  list_ret = list()
  for ( ts_idx in 1:length(list_traj) )
    {
    df = list_traj[[ts_idx]]
    #
    # Check if the df has enough rows = timsteps to be lagged
    #
    if (nrow (df) <= tau_max)
      {
      if (!keep_max_data)
        {
        miic_warning ("data lagging", "the trajectory ", ts_idx, " has only ",
          nrow (df), " time steps and will be ignored.")
        # NB : the number of columns is the same as the initial df !
        list_ret[[ts_idx]] = df[FALSE,]
        next
        }
      miic_warning ("data lagging", "the trajectory ", ts_idx, " has only ",
        nrow (df), " time steps and can not be lagged over ", tau_max,
        " time steps back.")
      }
    #
    # Lag the df
    #
    list_tmp = list()
    for ( var_idx in 1:nrow (state_order) )
      {
      if (state_order[var_idx, "lag"] == 0)
        list_tmp[[var_idx]] = df[,(var_idx+1)]
      else
        {
        max_row = nrow(df) - state_order[var_idx, "lag"]
        if (max_row <= 0)
          list_tmp[[var_idx]] = rep (NA, nrow(df) )
        else
          list_tmp[[var_idx]] = c ( rep (NA, state_order[var_idx, "lag"]),
                                    df [1:max_row,
                                        state_order[var_idx, "var_idx_data"]+1] )
        }
      }
    names(list_tmp) = state_order$var_names
    #
    # do.call (cbind, ...) was used for speed but returns a matrix
    # => problem when different types are used.
    # data.frame() preservers data type and seems as fast
    #
    # df <- as.data.frame (do.call (cbind, list_tmp) )
    df <- data.frame (list_tmp)
    # for (var_idx in 1:ncol(df))
    #   print (paste0 ("var : ", state_order[var_idx, "var_names"],
    #                 ", class :", class (df[,var_idx])))

    if (!keep_max_data)
      df = df [(tau_max+1):nrow(df),]
    #
    # Check rows with only NAs
    #
    rows_only_na <- ( rowSums (is.na (df)) == ncol (df) )
    df <- df [!rows_only_na, ]
    na_count = na_count + sum (rows_only_na)

    list_ret[[ts_idx]] = df
    }
  if (na_count > 0)
    miic_warning ("data lagging", "the lagged data contains ", sum(na_count),
             " row(s) with only NAs. These row(s) have been removed.")
  n_tot_samples <- sum (unlist (lapply (list_ret, FUN=nrow) ) )
  if (n_tot_samples <= 0)
    miic_error ("data lagging", "the data lagging produced no sample.",
                " Consider to review the data and/or the temporal window used.")
  if (n_tot_samples <= 1)
    miic_error ("data lagging", "the data lagging produced only 1 sample.",
                " Consider to review the data and/or the temporal window used.")
  return (list_ret)
  }

#===============================================================================
# FUNCTIONS (exported)
#===============================================================================
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
#' The expected data frame layout is variables as columns and
#' time series/time steps as rows.
#' The time step information must be supplied in the first column and,
#' for each time series, be consecutive (increment of 1) and in ascending order.
#' Multiple trajectories can be provided, the function will consider that a
#' new trajectory starts each time a smaller time step than the one of the
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
  input_data <- check_input_data (input_data, "TS")
  state_order <- check_state_order (input_data, state_order, "TS")
  state_order$n_layers <- NULL
  state_order$delta_t <- NULL
  state_order <- tmiic_check_state_order_part1 (state_order)
  list_ret <- tmiic_check_parameters (state_order = state_order,
                                     params = list(),
                                     n_layers = NULL,
                                     delta_t = NULL,
                                     movavg = movavg,
                                     keep_max_data = F,
                                     max_nodes = max_nodes)
  state_order <- tmiic_check_state_order_part2 (list_ret$state_order)

  list_traj <- tmiic_extract_trajectories (input_data)
  list_traj <- tmiic_movavg (list_traj, state_order$movavg, verbose_level=verbose_level)

  state_order <- tmiic_stat_estimate_dynamic (list_traj, state_order, max_nodes=max_nodes,
                                              verbose_level=verbose_level)
  n_layers <- unique (state_order$n_layers[ (state_order$var_type == 1)
                                         & (state_order$is_contextual == 0) ])
  delta_t <- unique (state_order$delta_t[ (state_order$var_type == 1)
                                         & (state_order$is_contextual == 0) ])
  return ( list ("n_layers"=n_layers, "delta_t"=delta_t))
  }
