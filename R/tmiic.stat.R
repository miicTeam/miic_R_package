#*******************************************************************************
# Filename   : tmiic.stat.R                         Creation date: 31 jan 2024
#
# Description: Utility functions for temporal MIIC in stationary mode
#
# Author     : Franck SIMON
#*******************************************************************************

#===============================================================================
# FUNCTIONS
#===============================================================================
# tmiic_ajust_window_for_nb_samples
#-------------------------------------------------------------------------------
# Reduce the window size (n_layers or delta_t) if the foreseen n_layers and
# delta_t would lead to too few samples after data lagging
# params:
# - list_traj: a list of data frames, each item representing a trajectory.
# - n_layers: a list, the n_layers in the state_order column
# - delta_t: a list, the delta_t in the state_order column
# - reduced_param: a string, can be "n_layers" or "delta_t". Indicates which
#   parameter will be reduced if the number of samples is too small
# - verbose: boolean, if TRUE, display a message if the window size is reduced
# returns:
# - a list: the n_layers or delta_t, depending of the reduced_param value.
#   The value are possibly decreased to reduce the window size
#-------------------------------------------------------------------------------
tmiic_ajust_window_for_nb_samples <- function (list_traj, n_layers, delta_t,
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
    # Look for the best value to reach the target recursively
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
    # Max time steps back in time found, try to reduce n_layers or delta_t
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
    else if (verbose >= 1)
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
# Estimate tau (the number of total time steps back to cover the dynamic,
# the number of layers and delta t parameters from the data
# - list_traj: list of data frames, each item representing a trajectory.
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
# - verbose: integer in the range [0,2], 1 by default. The level of verbosity:
#   0 = no display, 1 = summary display, 2 = maximum display.
#-------------------------------------------------------------------------------
tmiic_estimate_dynamic <- function (list_traj, state_order, max_nodes=50,
                                         verbose=1)
  {
  # If n_layers and delta_t all defined, nothing to do (need only to return
  # the state_order item as this case can only happen when called by tMiicStat)
  #
  if (  ( ! any (is.na (state_order$n_layers) ) )
     && ( ! any (is.na (state_order$delta_t ) ) ) )
    return (list ("state_order" = state_order) )
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
  if (verbose >= 1)
    miic_msg ("Estimating the temporal dynamic ...")
  n_ts <- length (list_traj)
  n_vars_tot <- ncol (list_traj[[1]]) - 1
  n_vars_ctx <- sum (state_order$is_contextual)
  n_vars_lag <- n_vars_tot - n_vars_ctx
  #
  # Remove time step, contextual and discrete variables
  #
  if ( ! any ( (state_order$is_contextual == 0) & (state_order$var_type == 1) ) )
    miic_error ("dynamic estimation",
                "no variable to estimate the temporal dynamic",
                " (all variables are discrete or contextual).",
                " Consider specifying the n_layers and delta_t parameters.")
  for (ts_idx in 1:n_ts)
    {
    if (nrow (list_traj[[ts_idx]]) == 1)
      miic_warning ("dynamic estimation", "trajectory ", ts_idx,
                    " with only 1 time step is ignored for dynamic estimation.")
    list_traj[[ts_idx]] <- list_traj[[ts_idx]][,
      c(F, ( (state_order$is_contextual == 0) & (state_order$var_type == 1) ) ),
      drop=F]
    }
  #
  # Compute mean alpha per variable
  #
  n_vars <- ncol (list_traj[[1]])
  var_names <- colnames(list_traj[[1]])
  length_to_test <- min (unlist (lapply (list_traj, FUN=function (x) {
    ifelse ( nrow(x) <= 1, NA, nrow(x) ) } ) ), na.rm=T)
  mat_lags_vanish <- matrix (NA_integer_, nrow=n_ts, ncol=n_vars)
  colnames(mat_lags_vanish) <- var_names
  mat_lags_4_alpha <- matrix (NA_integer_, nrow=n_ts, ncol=n_vars)
  colnames(mat_lags_4_alpha) <- var_names
  mat_alphas <- matrix (NA_real_, nrow=n_ts, ncol=n_vars)
  colnames(mat_alphas) <- var_names
  alphas_per_var <- rep (NA_real_, n_vars)
  names(alphas_per_var) <- var_names
  taus_per_var <- rep (NA_integer_, n_vars)
  names(taus_per_var) <- var_names
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
      mat_lags_vanish[ts_idx, var_idx] <- acf_res$lag[min (acf_vanish), 1, 1]
      mat_lags_4_alpha[ts_idx, var_idx] <- max ( 1, round (mat_lags_vanish[ts_idx, var_idx] / 2) )
      mat_alphas[ts_idx, var_idx] <- acf_res$acf[
        mat_lags_4_alpha[ts_idx, var_idx]+1,1,1] ^ (1/mat_lags_4_alpha[ts_idx, var_idx])
      }
    alphas_per_var[[var_idx]] <- mean (mat_alphas[,var_idx], na.rm=T)
    taus_per_var[[var_idx]] <- round ( (1+alphas_per_var[[var_idx]])
                                     / (1-alphas_per_var[[var_idx]]) )
    }
  if (verbose >= 2)
    {
    miic_msg ("Tau per variable:")
    for (i in 1:length (var_names))
      miic_msg ("- ", var_names[[i]], ": ", taus_per_var[[i]])
    }
  #
  # Compute alphas range and deduce taus range
  #
  tau_min  <- max ( 1, min (taus_per_var, na.rm=T) )
  tau_mean <- max ( 1, round (mean (taus_per_var, na.rm=T), 0) )
  tau_max  <- min ( length_to_test, max  (taus_per_var, na.rm=T) )
  tau_max_kept <- min (length_to_test, tau_max, tau_mean * 2)
  tau <- tau_max_kept
  if (verbose >= 1)
    miic_msg ("Automatic estimation of parameters:\n",
      "- Relaxation times goes from ", tau_min, " to ", tau_max,
      " with a mean of ", tau_mean, ", tau max considered = ", tau_max_kept)
  #
  # We know tau : the average maximum time steps back in time to use for the
  # temporal discovery. Now estimate the number of layers 'n_layers'
  # and/or number of time steps between two layers 'delta_t'
  #
  n_layers = NULL
  delta_t = NULL
  if (  all (!is.na (state_order$n_layers)) ) # n_layers known => NAs in delta_t
    {
    state_order$delta_t[is.na(state_order$delta_t)] <- max ( 1,
        ceiling (tau / (state_order$n_layers[is.na(state_order$delta_t)] - 1)) )

    state_order$delta_t <- tmiic_ajust_window_for_nb_samples (list_traj,
      state_order$n_layers, state_order$delta_t, reduced_param="delta_t",
      verbose=verbose)

    uniq_n_layers <- unique (state_order$n_layers[ (state_order$is_contextual == 0)
                                                & (state_order$var_type == 1)] )
    uniq_delta_t <- unique (state_order$delta_t[ (state_order$is_contextual == 0)
                                              & (state_order$var_type == 1)] )
    if (verbose >= 1)
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

    state_order$n_layers <- tmiic_ajust_window_for_nb_samples (list_traj,
      state_order$n_layers, state_order$delta_t, reduced_param="n_layers",
      verbose=verbose)

    uniq_n_layers <- unique (state_order$n_layers[ (state_order$is_contextual == 0)
                                                 & (state_order$var_type == 1)] )
    uniq_delta_t <- unique (state_order$delta_t[ (state_order$is_contextual == 0)
                                               & (state_order$var_type == 1)] )
    if (verbose >= 1)
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

    state_order$delta_t <- tmiic_ajust_window_for_nb_samples (list_traj,
      state_order$n_layers, state_order$delta_t, reduced_param="delta_t",
      verbose=verbose)
    delta_t <- unique (state_order$delta_t[ (state_order$var_type == 1)
                                          & (state_order$is_contextual == 0) ])

    if (verbose >= 1)
      miic_msg ("- For a final graph with a target of ", max_nodes,
        " nodes having ", n_vars_lag, " lagged variables",
        ifelse (n_vars_ctx > 0, paste0 ("\n  and ", n_vars_ctx, " contextual variables"), ""),
        ":\n  ", n_layers,  " layers spaced by ", delta_t, " time steps",
        ", dynamic covered goes over t, t-", delta_t,
        ifelse (n_layers > 3, ", ...", ""),
        ifelse (n_layers > 2, paste0 (", t-", tau), "") )
    }

  ret = list ("state_order" = state_order,
              "n_layers" = n_layers,
              "delta_t" = delta_t,
              "lags_vanish" = mat_lags_vanish,
              "lags_alphas" = mat_lags_4_alpha,
              "alphas" =  mat_alphas,
              "alphas_mean" = alphas_per_var,
              "taus" = taus_per_var)
  return (ret)
  }

#-------------------------------------------------------------------------------
# tmiic_lag_state_order
#-------------------------------------------------------------------------------
# Modify the state order into a lagged version: the lagged variables are
# completed and/or repeated with lagX to match the lagged temporal graph.
# Params:
# - state_order: a data frame, the state order returned by
#   tmiic_check_state_order
# Returns: a data frame: the lagged state_order
#-------------------------------------------------------------------------------
tmiic_lag_state_order <- function (state_order)
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
      state_lagged [var_idx, "var_names"] <- paste0 (state_order [var_idx, "var_names"], "_lag0")
      state_lagged [var_idx, "lag"] <- 0
      state_lagged [var_idx, "var_idx_data"] <- var_idx
      }
    else
      {
      state_lagged [var_idx, "lag"] <- 0
      state_lagged [var_idx, "var_idx_data"] <- var_idx
      }
    }
  #
  # Duplicate rows for lagged variables
  #
  state_lagged_nrows <- nrow (state_lagged)
  n_layers_back_max <- max ( (state_lagged$n_layers - 1) )
  n_layers_back_idx <- 1
  for (n_layers_back_idx in 1:n_layers_back_max)
    {
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
# tmiic_lag_bb_te
#-------------------------------------------------------------------------------
# Modify the complementary df int a lagged version: the 3 column data frames are
# transformed into a 2 columns one, in which variables are transformed into
# their lagged representation. e.g.
# - normal_var1 - normal_var2 - 1 becomes normal_var1_lag1 - normal_var2_lag0
# - ctx_var1 - normal_var2 - NA becomes ctx_var1 - normal_var2_lag0
# Params:
# - state_order: a data frame, the state order returned by
#   tmiic_check_state_order
# - df: the data frame to transform in its lagged version
# Returns: a data frame: the lagged data frame
#-------------------------------------------------------------------------------
tmiic_lag_bb_te <- function (state_order, df)
  {
  if ( is.null (df) )
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

#-------------------------------------------------------------------------------
# tmiic_lag_input_data
#-------------------------------------------------------------------------------
# Reorganizes the inputs in a format usable by miic: input data are lagged
# using the history to create lagged variables
# The function slices the input data according to the information supplied in
# the state_order n_layers and delta_t.
#
# The number of variables is increased and renamed on n_layers
# layers by delta_t steps.
# e.g. with n_layers=3 and delta_t=3 : var1, var2 =>
# var1_lag0, var2_lag0, var1_lag3, var2_lag3, var1_lag6, var2_lag6.
#
# Every time step (until number of time steps - (n_layers  - 1) * delta_t)
# is converted into a sample in the lagged data.
#
# Example with n_layers=3 and delta_t=3:
#
# Time step Var & value    Var & value  => Sample  Var & value   Var & value
#   t-6     Var1_val(t-6) Var2_val(t-6) =>   i    Var1_lag6_val Var2_lag6_val
#   t-3     Var1_val(t-3) Var2_val(t-3) =>   i    Var1_lag3_val Var2_lag3_val
#    t       Var1_val(t)   Var2_val(t)  =>   i    Var1_lag0_val Var2_lag0_val
#
#   t-7     Var1_val(t-7) Var2_val(t-7) =>   i'   Var1_lag6_val Var2_lag6_val
#   t-4     Var1_val(t-4) Var2_val(t-4) =>   i'   Var1_lag3_val Var2_lag3_val
#   t-1     Var1_val(t-1) Var2_val(t-1) =>   i'   Var1_lag0_val Var2_lag0_val
#
#   t-8     Var1_val(t-8) Var2_val(t-8) =>   i"   Var1_lag6_val Var2_lag6_val
#   t-5     Var1_val(t-5) Var2_val(t-5) =>   i"   Var1_lag3_val Var2_lag3_val
#   t-2     Var1_val(t-2) Var2_val(t-2) =>   i"   Var1_lag0_val Var2_lag0_val
#
#   ...     ............. ............. => ...... ............. ............
#
# until number of time steps - (n_layers - 1) * delta_t is reached.
# The same process is applied to all input time series.
#
# Note that the lagging can be different for each input variable
# if different values of n_layers or delta_t are supplied and some
# variables can be not lagged at all like contextual ones.
#
# Params:
# - list_traj: the list of time series
# - state_order: a data frame, the lagged state order returned by
#   tmiic_lag_state_order
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
tmiic_lag_input_data <- function (list_traj, state_order, keep_max_data=FALSE)
  {
  tau_max <- max(state_order$lag)
  na_count <- 0
  list_ret <- list()
  for ( ts_idx in 1:length(list_traj) )
    {
    df <- list_traj[[ts_idx]]
    #
    # Check if the df has enough rows = timsteps to be lagged
    #
    if (nrow (df) <= tau_max)
      {
      if (!keep_max_data)
        {
        miic_warning ("data lagging", "the trajectory ", ts_idx, " has only ",
          nrow (df), " time steps and will be ignored.")
        # NB/TODO? : the number of columns is the same as the initial df !
        list_ret[[ts_idx]] <- df[FALSE,]
        next
        }
      #
      # When keep_max_data is T but not enough time steps to lag completely,
      # lag till the maximum available in this trajectory
      #
      miic_warning ("data lagging", "the trajectory ", ts_idx, " has only ",
        nrow (df), " time steps and can not be lagged over ", tau_max,
        " time steps back.")
      }
    #
    # Lag the df
    #
    list_tmp <- list()
    for ( var_idx in 1:nrow (state_order) )
      {
      if (state_order[var_idx, "lag"] == 0)
        list_tmp[[var_idx]] <- df[,(var_idx+1)]
      else
        {
        max_row <- nrow(df) - state_order[var_idx, "lag"]
        if (max_row <= 0)
          list_tmp[[var_idx]] <- rep (NA, nrow(df) )
        else
          list_tmp[[var_idx]] <- c ( rep (NA, state_order[var_idx, "lag"]),
                                     df [1:max_row,
                                         state_order[var_idx, "var_idx_data"]+1] )
        }
      }
    names(list_tmp) <- state_order$var_names
    #
    # do.call (cbind, ...) was used for speed but returns a matrix
    # => problem when different types are used.
    # data.frame() preserves data type and seems as fast
    #
    # df <- as.data.frame (do.call (cbind, list_tmp) )
    df <- data.frame (list_tmp)
    # for (var_idx in 1:ncol(df))
    #   print (paste0 ("var : ", state_order[var_idx, "var_names"],
    #                 ", class :", class (df[,var_idx])))
    if (!keep_max_data)
      df <- df [(tau_max+1):nrow(df), , drop=F]
    #
    # Check rows with only NAs
    #
    rows_only_na <- ( rowSums (is.na (df)) == ncol (df) )
    df <- df [!rows_only_na, , drop=F]
    na_count <- na_count + sum (rows_only_na)

    list_ret[[ts_idx]] <- df
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
#' Estimation of the temporal stationary causal discovery parameters.
#'
#' @description This function estimates the number of layers and number of
#' time steps between each layer that are needed to cover the dynamic of a
#' stationary temporal dataset when reconstructing a temporal causal graph.
#' Using autocorrelation decay, the function computes the average relaxation
#' time of the variables and, based on a maximum number of nodes, deduces the
#' number of layers and number of time steps between each layer to be used.
#'
#' @param input_data [a data frame]
#' A data frame containing the observational data.\cr
#' The expected data frame layout is variables as columns and
#' time series/time steps as rows.
#' The time step information must be supplied in the first column and,
#' for each time series, be consecutive and in ascending order (increment of 1).
#' Multiple trajectories can be provided, the function will consider that a
#' new trajectory starts each time a smaller time step than the one of the
#' previous row is encountered.
#'
#' @param state_order [a data frame] An optional data frame providing extra
#' information about variables. It must have d rows where d is the number of
#' input variables, excluding the time step one.\cr
#' For optional columns, if they are not provided or contain missing
#' values, default values suitable for \emph{input_data} will be used.
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
#' "mov_avg" (optional) contains an integer value that specifies the size of
#' the moving average window to be applied to the variable.
#' Note that if "mov_avg" column is present in the \emph{state_order},
#' its values will overwrite the function parameter.
#'
#' @param mov_avg [an integer] Optional, NULL by default.\cr
#' When an integer>= 2 is supplied, a moving average operation is applied
#' to all the non discrete and not contextual variables.
#' If no \emph{state_order} is provided, the discrete/continuous variables
#' are deduced from the input data.
#' If you want to apply a moving average only on specific columns, consider
#' to use a \emph{mov_avg} column in the \emph{state_order} parameter.
#'
#' @param max_nodes [a positive integer] The maximum number of nodes in the
#' final time-unfolded causal graph. The more nodes allowed in the temporal
#' causal discovery, the more precise will be the discovery but at the cost
#' of longer execution time. The default is set to 50 for fast causal
#' discovery. On recent computers, values up to 200 or 300 nodes are usually
#' possible (depending on the number of trajectories and time steps in the
#' input data).
#'
#' @param verbose [an integer value in the range [0,2], 1 by default]
#' The level of verbosity: 0 = no display, 1 = summary display, 2 = full
#' display.
#'
#' @return A named list with :
#' \itemize{
#'  \item{\emph{n_layers}: the number of layers. }
#'  \item{\emph{delta_t}: the number of time steps between the layers. }
#' }
#'
#' These extra items are also available for more details on how \emph{n_layers}
#' and \emph{delta_t} have been estimated:
#'
#' \itemize{
#'  \item{\emph{lags_vanish}: the matrix (trajectories * variables) of lags
#'  where autocorrelaton vanishes. }
#'  \item{\emph{lags_alphas}: the matrix (trajectories * variables) of lags
#'  used to estimate the alphas. }
#'  \item{\emph{alphas}: the matrix (trajectories * variables) of alphas. }
#'  \item{\emph{alphas_mean}: the mean of alphas for each variable. }
#'  \item{\emph{taus}: the relaxation time for each variable. }
#' }
#'
#' @export
#-------------------------------------------------------------------------------
estimateTemporalDynamic <- function (input_data, state_order=NULL,
                                     mov_avg=NULL, max_nodes=50, verbose=1)
  {
  # We check/prepare the data almost as prepare_inputs and tmiic_prepare_inputs
  # would do until we get the list of trajectories from the input_data
  #
  list_ret = list()
  list_ret$input_data <- check_input_data (input_data, "TS")
  list_ret$params <- check_parameters (input_data = list_ret$input_data,
                                       n_threads = 1,
                                       cplx = "nml",
                                       orientation = TRUE,
                                       ort_proba_ratio = 1,
                                       ort_consensus_ratio = NULL,
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
                                       mode = "TS",
                                       verbose = verbose)

  list_ret$state_order <- check_state_order (
    list_ret$input_data, state_order, list_ret$params$mode)
  state_order$n_layers <- NULL
  state_order$delta_t <- NULL
  #
  # Extra checks needing several inputs
  #
  list_ret = check_cross_inputs (list_ret)
  #
  # Still as prepare_inputs to avoid any issue with the functions called
  # =>  move the non lagged inputs into a nested list 'non_lagged',
  #
  list_ret <- list ("params" = list_ret$params, "non_lagged" = list_ret)
  list_ret$non_lagged$params = NULL
  #
  # Check the temporal params with layer and delta_t not defined
  #
  list_ret <- tmiic_check_parameters (
    list_in = list_ret,
    n_layers = NULL,
    delta_t = NULL,
    mov_avg = mov_avg,
    keep_max_data = F,
    max_nodes = max_nodes)
  #
  # Still as prepare_inputs to avoid any issue with the functions called
  # Init / check (/ harmonize if possible) n_layers, delta_t and mov_avg
  #
  list_ret <- tmiic_check_state_order (list_in = list_ret)
  #
  # Prepare trajectories
  #
  list_traj <- tmiic_extract_trajectories (list_ret$non_lagged$input_data)
  list_traj <- tmiic_mov_avg (list_traj, list_ret$non_lagged$state_order$mov_avg,
    keep_max_data=list_ret$params$keep_max_data, verbose=list_ret$params$verbose)
  #
  # We have the trajectories, we can estimate the dynamic
  #
  ret <- tmiic_estimate_dynamic (list_traj,
    list_ret$non_lagged$state_order, max_nodes=list_ret$params$max_nodes,
    verbose=list_ret$params$verbose)

  return ( ret [c ("n_layers", "delta_t", "lags_vanish", "lags_alphas",
                   "alphas", "alphas_mean", "taus") ] )

  }

#-------------------------------------------------------------------------------
# tMiicStat
#-------------------------------------------------------------------------------
#' tMiicStat, temporal version of miic to learn temporal causal networks
#' including latent variables from stationary time series.
#'
#' @description tMiicStat extends the miic method (Multivariate
#' Information-based Inductive Causation) to stationary time series.
#' It combines constraint-based, information-theoretic approaches and
#' temporality to disentangle direct from indirect effects amongst correlated
#' contemporaneous or lagged variables, including cause-effect relationships
#' and the effect of unobserved latent causes.
#'
#' @details tMiicStat reorganizes the dataset using the \emph{n_layers} and
#' \emph{delta_t} parameters (that are estimated automatically if not
#' supplied) to transform the time steps into lagged samples.
#' As starting point, a lagged graph is created with only edges having at
#' least one node laying on the last time step.
#' Then, miic standard algorithm is applied to remove dispensable edges.
#' The remaining edges are then duplicated to ensure time invariance
#' (stationary dynamic) and oriented using the temporality and the
#' signature of causality in observational data. The use of temporal mode
#' is presented in Simon 2024.
#'
#' tMiicStat is mostly compatible with the classsical miic method and, as of,
#' it relies on information theoretic principles which replace (conditional)
#' independence tests as described in Affeldt 2015, Cabeli 2020,
#' Cabeli 2021 and Ribeiro-Dantas 2024. It deals with both categorical and
#' continuous variables by performing optimal context-dependent discretization.
#' As such, the input data frame may contain both numerical columns which will
#' be treated as continuous, or character / factor columns which will be treated
#' as categorical. For further details on the optimal discretization method and
#' the conditional independence test, see the function discretizeMutual.
#' The user may also choose to run tMiicStat with scheme presented in Li 2019
#' and Ribeiro-Dantas 2024 to improve the end result's interpretability
#' by ensuring consistent separating sets.
#'
#' As tMiicStat provides the same optional features as miic, most of the miic
#' parameters can be used exactly in the same way and are not described
#' extensively here. Detailed information about the common parameters is
#' available in the \code{\link{miic}} documentation.
#'
#' @seealso \code{\link{miic}} for the non temporal method
#' and \code{\link{discretizeMutual}} for optimal discretization and
#' (conditional) independence test.
#'
#' @references
#' \itemize{
#' \item Simon \emph{et al.}, eLife 2024, \href{https://www.biorxiv.org/content/10.1101/2024.02.06.579177v1.abstract}{CausalXtract: a flexible pipeline to extract causal effects from live-cell time-lapse imaging data}
#' \item Ribeiro-Dantas \emph{et al.}, iScience 2024, \href{https://arxiv.org/pdf/2303.06423}{Learning interpretable causal networks from very large datasets, application to 400,000 medical records of breast cancer patients}
#' \item Cabeli \emph{et al.}, NeurIPS 2021, \href{https://why21.causalai.net/papers/WHY21_24.pdf}{Reliable causal discovery based on mutual information supremum principle for finite dataset}
#' \item Cabeli \emph{et al.}, PLoS Comput. Biol. 2020, \href{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1007866}{Learning clinical networks from medical records based on information estimates in mixed-type data}
#' \item Li \emph{et al.}, NeurIPS 2019, \href{http://papers.nips.cc/paper/9573-constraint-based-causal-structure-learning-with-consistent-separating-sets.pdf}{Constraint-based causal structure learning with consistent separating sets}
#' \item Verny \emph{et al.}, PLoS Comput. Biol. 2017, \href{https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005662}{Learning causal networks with latent variables from multivariate information in genomic data}
#' \item Affeldt \emph{et al.}, UAI 2015, \href{https://auai.org/uai2015/proceedings/papers/293.pdf}{Robust Reconstruction of Causal Graphical Models based on Conditional 2-point and 3-point Information}
#' }
#'
#' @param input_data [a data frame, required]
#'
#' A n*(1+d) data frame (n samples, 1 column for the time steps, d variables)
#' that contains the observational data.
#'
#' The expected data frame layout is variables as columns
#' and time series/time steps as rows.
#' The time step information must be supplied in the first column and,
#' for each time series, be consecutive and in ascending order (increment of 1).
#' Multiple trajectories can be provided, tMiicStat will consider that a new
#' trajectory starts each time a smaller time step than the one of the
#' previous row is encountered.
#'
#' Each extra column after the first corresponds to one variable.
#' The column names correspond to the names of the observed variables.
#' By default, after the lagging step, numeric columns with at least 5
#' distinct values will be treated as continuous whilst numeric columns
#' with less than 5 distinct values, factors and characters will be
#' considered as categorical.
#'
#' @param state_order [a data frame, optional, NULL by default]
#'
#' A data frame providing extra information for variables. It is expected to
#' have d rows where d is the number of input variables and possible columns
#' are described below. If some variables are missing and, for optional columns,
#' if they are not provided or contain missing values, default values suitable
#' for \emph{input_data} will be used.
#'
#' \emph{"var_names"} (required) contains the name of each variable as specified
#' by colnames(input_data)[2:(d+1)]. The time steps column is not considered
#' as a variable and should not be mentioned in the variables list.
#'
#' \emph{"var_type"} (optional) contains a binary value that specifies if each
#' variable is to be considered as discrete (0) or continuous (1).
#'
#' \emph{"levels_increasing_order"} (optional) contains a single character
#' string with all of the unique levels of the ordinal variable in
#' increasing order, delimited by comma ','. It will be used during
#' the post-processing to compute the sign of an edge using Spearman's rank
#' correlation. If a variable is continuous or is categorical but not ordinal,
#' this column should be NA.
#'
#' \emph{"is_contextual"} (optional) contains a binary value that specifies
#' if a variable is to be considered as a contextual variable (1) or not (0).
#' Contextual variables cannot be the child node of any other variable (cannot
#' have edge with arrowhead pointing to them).
#'
#' \emph{"is_consequence"} (ignored) the consequence prior is not compatible
#' with tMiicStat, such column will be ignored with a warning.
#'
#' Several other columns are possible in temporal mode:
#'
#' \emph{"n_layers"} (optional) contains an integer value that specifies the
#' number of layers to be considered for the variable.
#' Note that if a \emph{"n_layers"} column is present in the \emph{state_order},
#' its values will overwrite the function parameter.
#'
#' \emph{"delta_t"} (optional) contains an integer value that specifies the
#' number of time steps between each layer for the variable.
#' Note that if a \emph{"delta_t"} column is present in the \emph{state_order},
#' its values will overwrite the function parameter.
#'
#' \emph{"mov_avg"} (optional) contains an integer value that specifies the size
#' of the moving average window to be applied to the variable.
#' Note that if \emph{"mov_avg"} column is present in the \emph{state_order},
#' its values will overwrite the function parameter.
#'
#' @param true_edges [a data frame, optional, NULL by default]
#'
#' A data frame containing the edges of the true graph for computing
#' performance after the run.\cr
#' The expected layout is a three columns data frame,
#' with the first two columns being variable names and the third the lag.
#' Variables names must exist in the \emph{input_data} data frame and the lag
#' must be valid in the time unfolded graph. e.g. a row var1, var2, 3 is valid
#' with \emph{n_layers} = 4 + \emph{delta_t} = 1 or
#' \emph{n_layers} = 2 + \emph{delta_t} = 3
#' but not for \emph{n_layers} = 2 + \emph{delta_t} = 2 as there is no matching
#' edge in the time unfolded graph.\cr
#' Please note that the order is important: "var1 var2 3" is interpreted as
#' var1_lag3 -> var2_lag0.
#' Please note also that, for contextual variables that are not lagged,
#' the expected value in the third column for the time lag is NA.
#'
#' @param black_box [a data frame, optional, NULL by default]
#'
#' A data frame containing pairs of variables that will be considered
#' as independent during the network reconstruction. In practice, these edges
#' will not be included in the skeleton initialization and cannot be part of
#' the final result.\cr
#' The expected layout is a three columns data frame,
#' with the first two columns being variable names and the third the lag.
#' Variables names must exist in the \emph{input_data} data frame and the lag
#' must be valid in the time unfolded graph. e.g. a row var1, var2, 3 is valid
#' with \emph{n_layers} = 4 + \emph{delta_t} = 1 or
#' \emph{n_layers} = 2 + \emph{delta_t} = 3
#' but not for \emph{n_layers} = 2 + \emph{delta_t} = 2 as there is no matching
#' edge in the time unfolded graph.\cr
#' Please note that the order is important: var1, var2, 3 is interpreted as
#' var1_lag3 - var2_lag0.
#' Please note also that, for contextual variables that are not lagged,
#' the expected value in the third column for the time lag is NA.
#'
#' @param n_threads [a positive integer, optional, 1 by default, see \code{\link{miic}}]
#' @param cplx [a string, optional, "nml" by default, possible values:
#' "nml", "bic", see \code{\link{miic}}]
#' @param orientation [a boolean value, optional, TRUE by default, see \code{\link{miic}}]
#' @param ort_proba_ratio [a floating point between 0 and 1, optional,
#' 1 by default, see \code{\link{miic}}]
#' @param ort_consensus_ratio [a floating point between 0 and 1, optional,
#' NULL by default, see \code{\link{miic}}]
#' @param propagation [a boolean value, optional, FALSE by default, see \code{\link{miic}}]
#' @param latent [a string, optional, "orientation" by default, possible
#' values: "orientation", "no", "yes", see \code{\link{miic}}]
#' @param n_shuffles [a positive integer, optional, 0 by default, see \code{\link{miic}}]
#' @param conf_threshold [a positive floating point, optional, 0 by default, see \code{\link{miic}}]
#' @param test_mar [a boolean value, optional, TRUE by default, see \code{\link{miic}}]
#' @param max_iteration [a positive integer, optional, 100 by default, see \code{\link{miic}}]
#' @param negative_info [a boolean value, optional, FALSE by default, see \code{\link{miic}}]
#' @param verbose [an integer value, optional, 1 by default, see \code{\link{miic}}]
#'
#' @param n_eff [a positive integer, optional, -1 by default]
#'
#' The \emph{n_eff}  parameter has a specific usage in tMiicStat.
#' In non temporal miic, \emph{n_eff} is the number of effective
#' samples and can be provided when dealing with correlated samples.
#' In temporal datasets, samples are expected to be correlated as the past of
#' each variable is likely correlated with its value of the next time step.
#' However, there is no correction to apply has the temporal autocorrelation
#' is taken into account during the lagged network reconstruction.
#' So, the typically case of use of the \emph{n_eff} parameter in temporal mode
#' is not the auto-correlation but, when the \emph{delta_t} value is > 1,
#' as the number of effective samples is divided by \emph{delta_t}
#' after the lagging process.
#' When set to its default (-1), tMiicStat will automatically adjust this value
#' to the total number of time steps / \emph{delta_t}.
#'
#' @param n_layers [an integer, optional, NULL by default, must be >= 2
#' if supplied]
#'
#' \emph{n_layers} defines the number of layers
#' that will be considered for the variables in the time unfolded graph.
#' The layers will be distant of \emph{delta_t} time steps.
#' If not supplied, the number of layers is estimated from the dynamic of the
#' dataset and the maximum number of nodes \emph{max_nodes} allowed in the
#' final lagged graph.
#'
#' @param delta_t [an integer, optional, NULL by default, must be >= 1
#' if supplied]
#'
#' \emph{delta_t} defines the number of time steps between each layer.
#' e.g. on 1000 time steps with \emph{n_layers} = 3 and \emph{delta_t} = 7,
#' the time steps kept for the samples conversion will be 1, 8, 15
#' for the first sample, the next sample will use 2, 9, 16 and so on.
#' If not supplied, the number of time steps between layers is estimated
#' from the dynamic of the dataset and the number of layers.
#'
#' @param mov_avg [an integer, optional, NULL by default, must be >= 2
#' if supplied]
#'
#' When supplied, a moving average operation is applied to all integer
#' and numeric variables that are not contextual variables.
#' If the moving average should only be applied on a part of the numeric
#' variables, the moving average can be specified in the \emph{state_order}
#' per variable.
#'
#' @param keep_max_data [a boolean value, optional, FALSE by default]
#'
#' If TRUE, rows where some NAs have been introduced during the moving averages
#' and lagging will be kept whilst they will be dropped if FALSE.
#'
#' @param max_nodes [an integer, optional, 50 by default]
#'
#' Used only if the \emph{n_layers} or \emph{delta_t}
#' parameters are not supplied. \emph{max_nodes} is used as the maximum number
#' of nodes in the final time-unfolded graph to compute \emph{n_layers}
#' and/or \emph{delta_t}.
#' The default is 50 to produce quick runs and can be increased up to 200
#' or 300 on recent computers to produce more precise results.
#'
#' @return As tMiicStat is the extension of miic to stationary time series,
#' the object returned is a \emph{miic-like} object enriched with
#' extra information.
#'
#' These following items describe the time unfolded network inferred by
#' tMiicStat. As they are identical to the ones returned by \emph{miic},
#' please refer to \code{\link{miic}} for their description:
#'
#' \itemize{
#'  \item{\emph{summary:} a data frame with information about the
#'  relationship between relevant pair of variables. }
#'
#'  \item{\emph{edges:} a data frame with the raw edges output coming from
#'  the C++ core function. }
#'
#'  \item{\emph{triples:} this data frame lists the orientation
#'  probabilities of the two edges of all unshielded triples of the
#'  reconstructed network. }
#'
#'  \item {\emph{adj_matrix:} the adjacency matrix is a square matrix used to
#'  represent the inferred graph. }
#'
#'  \item {\emph{proba_adj_matrix:} the probability adjacency matrix is
#'  a square matrix used to represent the orientation probabilities associated
#'  to the edges tips. }
#'
#'  \item {\emph{adj_matrices:} present only when consistency is activated.
#'  The list of the adjacency matrices, one for each graph
#'  which is part of the resulting cycle of graphs. }
#'
#'  \item {\emph{proba_adj_matrices:} present only when consistency is
#'  activated. The list of the probability adjacency matrices,
#'  one for each graph which is part of the resulting cycle of graphs. }
#'
#'  \item {\emph{proba_adj_average:} present only when consistency is activated.
#'  The average probability adjacency matrix is a square matrix used to
#'  represent the orientation probabilities associated to the edges tips
#'  of the consensus graph. }
#'
#'  \item {\emph{is_consistent:} present only when consistency is activated.
#'  TRUE if the returned graph is consistent, FALSE otherwise. }
#'
#'  \item {\emph{time:} execution time of the different steps and total run-time
#'  of the causal graph reconstruction. }
#'
#'  \item {\emph{interrupted:} TRUE if causal graph reconstruction has been
#'  interrupted, FALSE otherwise. }
#'
#'  \item {\emph{scores:} present only when true edges have been supplied.
#'  Contains the scores of the returned graph in regard of the ground truth. }
#'
#'  \item {\emph{params:} the list of parameters used for the network
#'  reconstruction. The parameters not supplied are initialized to their default
#'  values. Otherwise, the parameters are checked and corrected if necessary. }
#' }
#'
#' These following items are similar to the ones returned by miic,
#' but with a specific organization for tMiicStat:
#'
#' \itemize{
#'  \item {\emph{input_data:} the input data before the lagging process,
#'  checked and corrected if necessary.
#'  Please note that these input data correspond to the ones supplied
#'  as parameter of the tMiicStat function but is not the data used internally
#'  for the network reconstruction as it is not lagged.
#'  The lagged input data are available in the \emph{tmiic} sub list. }
#'
#'  \item {\emph{state_order:} the state order used before the lagging process.
#'  If no state order is supplied, it is generated by using default values.
#'  Otherwise, it is the state order checked and corrected if necessary.
#'  Please note that this state order layout corresponds to the one supplied
#'  (or can be supply) as parameter to the tMiicStat function but is not
#'  the state order used internally for the network reconstruction
#'  as it is not lagged.
#'  The lagged stater order is available in the \emph{tmiic} sub list. }
#'
#'  \item {\emph{black_box:} present only if a black box has been supplied,
#'  the black box, before the lagging process, checked and corrected
#'  if necessary.
#'  Please note that this black box layout corresponds to the one supplied
#'  as parameter to the tMiicStat function but is not the black box used
#'  internally for the network reconstruction as it is not lagged.
#'  The lagged black box is available in the \emph{tmiic} sub list. }
#'
#'  \item {\emph{true_edges:} present only if the true edges have been supplied,
#'  the true edges, before the lagging process, checked and corrected
#'  if necessary.
#'  Please note that this true edges layout corresponds to the one supplied
#'  as parameter to the tMiicStat function but is not the true edges used
#'  internally for the network reconstruction as it is not lagged.
#'  The lagged true edges are available in the \emph{tmiic} sub list. }
#' }
#'
#' In addition, tMiicStat provides extra information dedicated to the network
#' inferrence in temporal mode:
#'
#' \itemize{
#'   \item {\emph{tmiic:} named list containing:
#'   \itemize{
#'     \item {\emph{input_data:} the lagged data
#'     used internally to perform the time unfolded graph reconstruction. }
#'     \item {\emph{state_order:} the lagged state order
#'     used internally to perform the time unfolded graph reconstruction. }
#'     \item {\emph{black_box:} if a black box has been supplied, the lagged
#'     version used internally to perform the time unfolded graph
#'     reconstruction. }
#'     \item {\emph{true_edges:} if true edges have been supplied, the lagged
#'     version used internally to perform the time unfolded graph
#'     reconstruction. }
#'     \item {\emph{stationary:} the inferred network with the list of edges
#'     completed by stationarity. }
#'     }
#'   }
#' }
#'
#' @export
#' @useDynLib miic
#' @import Rcpp
#'
#' @examples
#' library(miic)
#'
#' # Example on Covid cases (time series toy demo)
#'
#' data(covidCases)
#' # execute MIIC (reconstruct graph in temporal mode)
#' tmiic_obj <- tMiicStat (input_data = covidCases,
#' n_layers = 3, delta_t = 1, mov_avg = 14)
#'
#' # to plot the default graph (compact)
#' if(require(igraph)) {
#'  plot(tmiic_obj)
#' }
#'
#' # to plot the raw temporal network (lagged)
#' if(require(igraph)) {
#'   plot(tmiic_obj, display="raw")
#' }
#'
#' # to plot the full temporal network  (lagged and completed by stationarity)
#' if(require(igraph)) {
#'   plot(tmiic_obj, display="lagged")
#' }
#-------------------------------------------------------------------------------
tMiicStat <- function (input_data,
                       state_order = NULL,
                       true_edges = NULL,
                       black_box = NULL,
                       n_threads = 1,
                       cplx = "nml",
                       orientation = TRUE,
                       ort_proba_ratio = 1,
                       ort_consensus_ratio = NULL,
                       propagation = FALSE,
                       latent = "orientation",
                       n_eff = -1,
                       n_shuffles = 0,
                       conf_threshold = 0,
#                       sample_weights = NULL,
                       test_mar = TRUE,
                       consistent = "no",
                       max_iteration = 100,
                       consensus_threshold = 0.8,
                       negative_info = FALSE,
                       n_layers = NULL,
                       delta_t = NULL,
                       mov_avg = NULL,
                       keep_max_data = FALSE,
                       max_nodes = 50,
                       verbose = 1)
  {
  return (miic_private (input_data = input_data,
                        state_order = state_order,
                        true_edges = true_edges,
                        black_box = black_box,
                        n_threads = n_threads,
                        cplx = cplx,
                        orientation = orientation,
                        ort_proba_ratio = ort_proba_ratio,
                        ort_consensus_ratio = ort_consensus_ratio,
                        propagation = propagation,
                        latent = latent,
                        n_eff = n_eff,
                        n_shuffles = n_shuffles,
                        conf_threshold = conf_threshold,
                        sample_weights = NULL,
                        test_mar = test_mar,
                        consistent = consistent,
                        max_iteration = max_iteration,
                        consensus_threshold = consensus_threshold,
                        negative_info = negative_info,
                        mode = "TS",
                        n_layers = n_layers,
                        delta_t = delta_t,
                        mov_avg = mov_avg,
                        keep_max_data = keep_max_data,
                        max_nodes = max_nodes,
                        verbose = verbose) )
  }

