#*******************************************************************************
# Filename   : tmiic.wrapper.R                    Creation date: 24 march 2020
#
# Description: Data transformation of time series for miic
#
# Author     : Franck SIMON
#
# Changes history:
# - 24 mar 2020 : initial version
# - 04 jun 2020 : add tmiic_flatten_network
# - 15 jun 2020 : add delta_tau and moving average
# - 27 jul 2020 : rewrite of tmiic.lag_inputs to allow variable
#                 number of time steps between time series
# - 06 oct 2020 : add tmiic_repeat_edges_over_history to duplicate the edges
#                 over history
# - 07 jan2021 : replace tmiic.transform_data_for_miic by tmiic.lag_inputs:
#                 allow different tau, delta_tau or movavg per variable
#                 add support for contextual variables (not lagged)
#                 transfer input data lagging into C++ function
#                 lag state_order, true_edges and black_box
# - 10 fev 2021 : add support in flattening and repeat of edges over history
#                 of different tau, delta_tau per variable
# - 09 aou 2023 : lag the inputs in the R part,
#                 change tau, delta_tau parameters into n_layers and delta_t
#*******************************************************************************

#-------------------------------------------------------------------------------
# tmiic_lag_state_order
#-------------------------------------------------------------------------------
# Modify the state order into a lagged version: the lagged variables are
# completed and/or repeated with lagX to match the lagged temporal graph.
# inputs:
# - state_order: a dataframe, the state order returned by
#   tmiic_check_state_order_part2
# Returns: a dataframe: the lagged state_order
#-------------------------------------------------------------------------------
tmiic_lag_state_order <- function (state_order)
  {
  n_vars <- nrow (state_order)
  state_order$lag <- -1
  state_order$initial_var <- -1
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
      state_lagged [var_idx, "initial_var"] = var_idx
      }
    else
      {
      state_lagged [var_idx, "lag"] = 0
      state_lagged [var_idx, "initial_var"] = var_idx
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
        state_lagged [state_lagged_nrows, "initial_var"] <- var_idx
        }
      }
    }
  return (state_lagged)
  }

#-------------------------------------------------------------------------------
# tmiic_lag_other_df
#-------------------------------------------------------------------------------
# Modify the complementary df int a lagged version: the 3 column dataframes are
# transformed into a 2 columns one, in which variables are transformed into
# their lagged representation. i.e:
# - lagged_var1 - lagged_var2 - 1 becomes lagged_var1_lag1 - lagged_var2_lag0
# - ctx_var1 - lagged_var2 - NA becomes ctx_var1 - lagged_var2_lag0
# inputs:
# - df: the dataframe to transform in its lagged version
# - state_order: a dataframe, the state order returned by
#   tmiic_check_state_order_part2
# Returns: a dataframe: the lagged dataframe
#-------------------------------------------------------------------------------
tmiic_lag_other_df <- function (state_order, df)
  {
  if ( (is.null (df)) || (nrow (df) <= 0) )
    return (df)

  for (i in 1:nrow (df))
    {
    orig_node_idx <- which (state_order$var_names == df[i, 1])
    if (state_order[orig_node_idx, "is_contextual"] == 0)
      df[i, 1] = paste0 (df [i, 1], "_lag", df [i, 3])
    df[i, 2] = paste0 (df [i, 2], "_lag0")
    }
  df <- df[,-3]
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
# layers by delta_t. steps.
# i.e. with n_layers=3 and delta_t.=3 : var1, var2 =>
# var1_lag0, var2_lag0, var1_lag3, var2_lag3, var1_lag6, var2_lag6.
#
# Every time step (until number of time steps - (n_layers  - 1) * delta_t.)
# is converted into a sample in the lagged data.
#
# Exemple with n_layers=3 and delta_t.=3:
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
# inputs:
# - list_ts: the list of time series
# - state_order: a dataframe, the lagged state order returned by
#   tmiic_lag_state_order
# - keep_max_data: boolean flag, optional, FALSE by default
#   When FALSE, the rows containing NA introduced by the lagging process
#   are deleted, otherwise when TRUE, the rows are kept
#-------------------------------------------------------------------------------
tmiic_lag_input_data <- function (list_ts, state_order, keep_max_data=FALSE)
  {
  tau_max = max(state_order$lag)
  na_count = 0
  list_ret = list()
  for ( ts_idx in 1:length(list_ts) )
    {
    df = list_ts[[ts_idx]]
    #
    # Check if the df has enough rows = timsteps to be lagged
    #
    if (nrow (df) <= tau_max)
      {
      if (!keep_max_data)
        {
        miic_warning ("data lagging", "the trajectory ", ts_idx, " has only ",
          nrow (df), " time steps and will be ignored.")
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
                                        state_order[var_idx, "initial_var"]+1] )
        }
      }
    names(list_tmp) = state_order$var_names
    # df <- as.data.frame (do.call (cbind, list_tmp) )
    df <- data.frame (list_tmp)
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
  return (list_ret)
  }

#-----------------------------------------------------------------------------
# tmiic_precompute_lags_layers_and_shifts
#-----------------------------------------------------------------------------
# Utility function to precompute lags, layers and shifts of nodes in the
# lagged network
#
# params: tmiic_res [a tmiic object] The object returned by miic's
# execution in temporal mode.
#
# returns: a dataframe with lagged nodes as row name and 3 columns:
#  - lags: the lag of each lagged node
#  - corresp_nodes: the corresponding non lagged node
#  - shifts: the shift to apply to find the next lagged node
#-----------------------------------------------------------------------------
tmiic_precompute_lags_layers_and_shifts <- function (tmiic_res)
  {
  list_nodes_not_lagged = tmiic_res$state_order$var_names
  n_nodes_not_lagged = length (list_nodes_not_lagged)
  list_n_layers_back <- tmiic_res$state_order$n_layers - 1
  list_delta_t <- tmiic_res$state_order$delta_t
  #
  # Identify lag and layer of each node
  #
  list_lags <- rep(0, n_nodes_not_lagged)
  list_nodes_lagged <- c()
  list_corresp_nodes <- c()
  i = 1
  for (node_idx in 1:n_nodes_not_lagged)
    {
    node_name <- list_nodes_not_lagged[[node_idx]]
    list_corresp_nodes[[i]] <- node_name

    if (list_n_layers_back[[node_idx]] >= 1)
      node_name <- paste0 (node_name, "_lag0")
    list_nodes_lagged [[i]] <- node_name
    i <- i + 1
    }

  n_layers_back_max <- max (list_n_layers_back)
  for (n_layers_back_max in 1:n_layers_back_max)
    {
    for (node_idx in 1:n_nodes_not_lagged)
      {
      n_layers_back_of_var <- list_n_layers_back[[node_idx]];
      if (n_layers_back_max <= n_layers_back_of_var)
        {
        node_name <- list_nodes_not_lagged[[node_idx]]
        list_corresp_nodes[[i]] <- node_name

        lag <- n_layers_back_max * list_delta_t[[node_idx]];
        node_name <- paste0 (node_name, "_lag", lag)
        list_nodes_lagged [[i]] <- node_name
        list_lags[[i]] <- lag
        i <- i + 1
        }
      }
    }
  #
  # Precompute the index shifts from a node to its first lagged counterpart
  #
  n_nodes_shifts = n_nodes_not_lagged
  end_reached = rep (FALSE, n_nodes_not_lagged)

  list_shifts <- c()
  for (n_layers_back_idx in 1:(n_layers_back_max+1) )
    {
    for (node_idx in 1:n_nodes_not_lagged)
      {
      n_layers_back_of_var <- list_n_layers_back[[node_idx]];
      if (n_layers_back_idx <= n_layers_back_of_var)
        list_shifts <- append (list_shifts, n_nodes_shifts)
      else if (!end_reached[[node_idx]])
        {
        end_reached[[node_idx]] = TRUE;
        list_shifts <- append (list_shifts, 0)
        n_nodes_shifts <- n_nodes_shifts - 1
        }
      }
    }

  df_ret <- data.frame (lags=as.integer(unlist(list_lags)),
                        corresp_nodes=unlist(list_corresp_nodes),
                        shifts=as.integer(unlist(list_shifts)),
                        stringsAsFactors=FALSE)
  rownames (df_ret) <- list_nodes_lagged
  return (df_ret)
  }

#-------------------------------------------------------------------------------
#  tmiic_combine_lag
#-------------------------------------------------------------------------------
# Utility function to combine lags  when flattening the network.
#
# param: df, a non empty dataframe with the edges to combine
#-------------------------------------------------------------------------------
tmiic_combine_lag <- function (df)
  {
  # Reverse inverted edges (orient == -2) and duplicate lags of bidrectional
  # temporal edges (lag != 0)
  # NB: for non lag 0 edges, such cases are more than likely errors
  # and will generate negative lags however showing them will allow the user
  # to identify possible issues
  #
  for (idx in 1:nrow(df) )
    {
    if (df[idx,"infOrt"] == -2)
      df[idx, c("x","y","lag")] <- c (df[idx,"y"], df[idx,"x"],
                                      -as.integer (df[idx,"lag"]) )

    if ( (df[idx,"infOrt"] == 6) & (as.integer (df[idx,"lag"]) != 0) )
      df[nrow(df)+1, c("x","y","lag")] <- c (df[idx,"y"], df[idx,"x"],
                                               -as.integer (df[idx,"lag"]) )
    }
  #
  # Combine lags from node1 < node2, lag0 and node1 >= node2
  #
  list_lag1 <- unique (df[ ( (df$x < df$y) & (df$lag != 0) ),]$lag)
  list_lag2 <- unique (df[ (df$lag == 0),]$lag)
  list_lag3 <- unique (df[ ( (df$x >= df$y) & (df$lag != 0) ),]$lag)

  list_lag1 <- paste (unlist (list_lag1[order (list_lag1, decreasing=TRUE)]), collapse=",")
  list_lag2 <- as.character (list_lag2)
  list_lag3 <- paste (unlist (list_lag3[order (list_lag3)]), collapse=",")

  list_lags <- c(list_lag1[[1]], list_lag2, list_lag3[[1]])
  list_lags <- lapply (list_lags, function(z) { z[!is.na(z) & z != ""]})
  if ( (list_lag1[[1]] != "") & (list_lag3[[1]] != "") )
    list_lags <- paste (unlist (list_lags), collapse="/")
  else
    list_lags <- paste (unlist (list_lags), collapse=",")

  return (list_lags)
  }

#-------------------------------------------------------------------------------
# tmiic_combine_orient
#-------------------------------------------------------------------------------
# Utility function to combine edges orientations when flattening the network.
#
# params:
# - df: a dataframe with the edges to combine
# - col_name: string, the orientation column
#-------------------------------------------------------------------------------
tmiic_combine_orient   <- function (df, col_name)
  {
  df <- df[!is.na (df[[col_name]]),]
  if (nrow (df) <= 0)
    return (NA)
  #
  # We set orientations as if node X <= node Y
  # NB: we do not care of x, y, lag and proba columns as they are not used
  # later in the function
  #
  for (idx in 1:nrow(df) )
    if ( (df[idx,"x"] > df[idx,"y"]) & (!is.na (df[idx, col_name])) )
        if (abs(df[idx, col_name]) == 2)
          df[idx, col_name] <- -(df[idx, col_name])

  col_min <- min (df[, col_name])
  col_max <- max (df[, col_name])

  if ( (col_max == 6) | ((col_min == -2) & (col_max == 2)) )
    return (6)
  else if (col_max == 2)
    return (2)
  else if (col_min == -2)
    return (-2)
  else
    return (1)
  }

#-----------------------------------------------------------------------------
# tmiic_combine_probas
#-----------------------------------------------------------------------------
# Utility function to combine edges probabilities when flattening the network.
# Depending on the combined edges orientation, chose the appropriate max, min
# or mean probabilities to compute the combined edge probabilities
#
# params:
# - df: the dataframe with the edges to combine
# - comb_orient: integer, the orientation of the combined edge
#-----------------------------------------------------------------------------
tmiic_combine_probas <- function (df, comb_orient)
  {
  valid_probas <- grepl (';', df$proba, fixed=TRUE)
  df <- df[valid_probas,]
  if (nrow (df) <= 0)
    return (NA)
  #
  # We set probas like if we have node X <= node Y
  #
  for ( idx in 1:nrow(df) )
    if (df[idx,"x"] > df[idx,"y"])
      {
      proba_split <- strsplit (df[idx, "proba"], ';' )[[1]]
      df[idx, "proba"] <- paste (proba_split[[2]], proba_split[[1]], sep=";")
      }
  #
  # Split proba column so we can do maths on it
  #
  probas_split <- strsplit (df$proba, ';' )
  df_probas <- do.call(rbind, probas_split)
  df_probas <- data.frame (x=as.numeric ( as.character (df_probas[,1]) ),
                           y=as.numeric ( as.character (df_probas[,2]) ) )
  #
  # Depending on the pre-computed combined orientation, keep max/min/avg
  #
  if (comb_orient == 6)
    return (paste (max(df_probas[,1]), max(df_probas[,2]), sep=";") )
  if (comb_orient == 2)
    return (paste (min(df_probas[,1]), max(df_probas[,2]), sep=";") )
  if (comb_orient == -2)
    return (paste (max(df_probas[,1]), min(df_probas[,2]), sep=";") )
  return (paste (mean(df_probas[,1]), mean(df_probas[,2]), sep=";"))
  }

#-----------------------------------------------------------------------------
# tmiic_flatten_network
#-----------------------------------------------------------------------------
# Flattten the lagged network returned by tmiic for plotting
#
# In temporal mode, the network returned by miic contains lagged nodes
# (X_lag0, X_lag1, ...). This function flatten the  network depending
# of the flatten_mode parameter.
# Note that only the summary data frame is flattened and the adjacency matrix
# is reduced to non lagged nodes and filled with NA during the process
#
# params:
# - tmiic_res: a tmiic object, returned by tmiic
#
# - flatten_mode: string, optional, default value "compact".
#   Possible values are "compact", "combine", "unique", "drop":
#   * "compact": the default. Nodes and edges are converted into a flattened
#     version preserving all the initial information.
#     i.e.: X_lag1->Y_lag0, X_lag0<-Y_lag2 become respectively X->Y lag=1,
#     X<-Y lag=2.
#   * "combine": one edge will be kept per couple of nodes.
#     The info_shifted will be the highest of the summarized edges whilst
#     the lag and orientation of the summarized edge will be an aggregation.
#     i.e.: X_lag2->Y_lag0, X_lag0<-Y_lag1 will become X<->Y lag=1,2 with
#     the info_shifted of X_lag2->Y_lag0 if info_shifted of
#     X_lag2->Y_lag0 > X_lag0<-Y_lag1.
#   * "unique": only the edges having the highest info_shifted for a couple
#     of nodes are kept in the flattened network. If several edges between
#     the sames nodes have the same info_shifted, then the edge kept is
#     the one with the minimum lag.
#     i.e.: X_lag1->Y_lag0, X_lag0<-Y_lag2 with info_shifted of
#     X_lag1->Y_lag0 > X_lag0<-Y_lag2 become X->Y lag=1.
#   * "drop"}, only the edges having the
#     highest info_shifted for a couple of nodes are kept in the flattened
#     network. If several edges between the sames nodes have the same
#     info_shifted, then the edge kept is the one with the minimum lag.\cr
#     i.e. :  X_lag1->Y_lag0, X_lag0<-Y_lag2 with info_shifted of
#     X_lag1->Y_lag0 > X_lag0<-Y_lag2 become X->Y. The lag information is
#     "lost" after flattening
#
#   Note that for all modes other than "drop", lag is a new column added
#   in the dataframe.
#
# - keep_edges_on_same_node: boolean, optional, TRUE by default.
#   When TRUE, the edges like X_lag0-X_lag1 are kept during flattening
#   (it becomes an X-X edge). When FALSE, only edges having different nodes
#   are kept in the flatten network.
#
# returns: a tmiic object. The returned tmiic object is the one received
# as input where the summary dataframe has been flattened and the adjacency
# matrix reduced to the non lagged nodes
#-----------------------------------------------------------------------------
tmiic_flatten_network <- function (tmiic_res, flatten_mode="compact",
                                   keep_edges_on_same_node=TRUE)
  {
  # Reduce size of adj_matrix to non lagged nodes
  # (we don't care about content as it is not used for plotting)
  #
  list_nodes <- tmiic_res$state_order$var_names
  tmiic_res$adj_matrix <- matrix(NA, nrow=0, ncol=length (list_nodes))
  colnames(tmiic_res$adj_matrix) <- list_nodes
  #
  # Keep only edges found by miic
  #
  df_edges <- tmiic_res$all.edges.summary[tmiic_res$all.edges.summary$type %in% c('P', 'TP', 'FP'), ]
  if (nrow(df_edges) <= 0)
    {
    if (flatten_mode != "drop")
      df_edges$lag = numeric(0)
    tmiic_res$all.edges.summary <- df_edges
    return (tmiic_res)
    }
  #
  # Precompute lag and layer of each node
  #
  df_precomputed <- tmiic_precompute_lags_layers_and_shifts (tmiic_res)
  #
  # First step, perform flatten_mode="compact":
  # from summary, remove lag info from nodes names and put it into a lag column
  #
  df_edges$lag <- -1
  for ( edge_idx in 1:nrow(df_edges) )
    {
    one_edge <- df_edges[edge_idx,]
    lag_x  <- df_precomputed [[one_edge$x, "lags"]]
    node_x <- df_precomputed [[one_edge$x, "corresp_nodes"]]
    lag_y  <- df_precomputed [[one_edge$y, "lags"]]
    node_y <- df_precomputed [[one_edge$y, "corresp_nodes"]]
    #
    # Ensure the lag is > 0 (=> we put nodes as x=oldest to y=newest)
    #
    lag = lag_x - lag_y
    if (lag >= 0)
      df_edges [edge_idx, c("x","y","lag")] <- c(node_x, node_y, lag)
    else
      {
      df_edges [edge_idx, c("x","y","lag")] <- c(node_y, node_x, -lag)
      if (abs (one_edge$infOrt) == 2)
         df_edges [edge_idx,"infOrt"] <- -one_edge$infOrt
      if ( !is.na (one_edge$trueOrt ) )
        if (abs (one_edge$trueOrt ) == 2)
          df_edges [edge_idx,"trueOrt"] <- -one_edge$trueOrt
      if ( !is.na (one_edge$proba ) )
        {
        df_edges [edge_idx, "proba"] = paste0 (rev (
          strsplit (df_edges [edge_idx, "proba"], ";")[[1]]), collapse=";")
        }
      }
    }
  df_edges <- transform ( df_edges, lag = as.integer (lag) )
  #
  # Exclude self loops if requested
  #
  if (!keep_edges_on_same_node)
    df_edges <- df_edges[df_edges$x != df_edges$y, ]
  if (nrow(df_edges) <= 0)
    {
    if (flatten_mode == "drop")
      df_edges$lag <- NULL
    tmiic_res$all.edges.summary <- df_edges
    return (tmiic_res)
    }
  #
  # "compact" mode is done
  #
  if (flatten_mode != "compact")
    {
    # if mode != "compact", we want only one edge per couple of nodes:
    # the edges kept per couple of nodes will be the one having the max
    # info_shifted and if several edges have the same info_shifted,
    # the one with the minimum lag.
    #
    # Identify the couples of same X-Y or Y-X whatever the lag or orientation
    #
    df_xy <- df_edges[,c("x", "y")]
    list_rows_to_swap <- (df_xy$x > df_xy$y)
    df_xy [list_rows_to_swap, c("x","y")] <- df_xy [list_rows_to_swap, c("y","x")]
    df_xy <- unique (df_xy)
    #
    # Keep one edge per couple of nodes
    #
    df_group <- df_edges[FALSE,]
    for ( xy_idx in 1:nrow(df_xy) )
      {
      ref_x <- df_xy[xy_idx,"x"]
      ref_y <- df_xy[xy_idx,"y"]
      cond_same_edges = ( ( (df_edges[["x"]] == ref_x) & (df_edges[["y"]] == ref_y) )
                        | ( (df_edges[["x"]] == ref_y) & (df_edges[["y"]] == ref_x) ) )
      df_same <- df_edges[cond_same_edges,]

      if (nrow (df_same) > 1)
        {
        if (flatten_mode == "combine")
          {
          # Combine lag, orient and proba
          #
          df_same$new_lag <-  tmiic_combine_lag (df_same)
          comb_infOrt <- tmiic_combine_orient (df_same, "infOrt")
          df_same$proba <- tmiic_combine_probas (df_same, comb_infOrt)
          df_same$trueOrt <- tmiic_combine_orient (df_same, "trueOrt")
          df_same$infOrt <- comb_infOrt
          #
          # Orientations and probas have been computed for x <= y,
          # so force x <= y on all rows
          #
          if (ref_x <= ref_y)
            {
            df_same[,"x"] <- ref_x
            df_same[,"y"] <- ref_y
            }
          else
            {
            df_same[, "x"] <- ref_y
            df_same[, "y"] <- ref_x
            }
          }

        max_info <- max (df_same[["info_shifted"]])
        df_same <- df_same[ (df_same[["info_shifted"]] == max_info),]
        }
      if (nrow(df_same) > 1)
        {
        min_lag <- min (df_same[["lag"]])
        df_same <- df_same[ (df_same[["lag"]] == min_lag),]
        }
      if ("new_lag" %in% colnames(df_same) )
        {
        df_same$lag <- df_same$new_lag
        df_same$new_lag <- NULL
        }
      df_group <- rbind (df_group, df_same)
      }
    df_edges <- df_group
    }
  #
  # Remove lag info when not wanted
  #
  if (flatten_mode == "drop")
    #
    # We do not want to keep info about lag at all
    #
    df_edges$lag <- NULL
  else
    {
    # For contextual variable, we clean the lag info
    #
    is_contextual <- tmiic_res$state_order$is_contextual
    if (!is.null(is_contextual))
      {
      list_nodes_not_lagged = tmiic_res$state_order$var_names
      for ( edge_idx in 1:nrow(df_edges) )
        {
        one_edge <- df_edges[edge_idx,]
        x_idx <- which (list_nodes_not_lagged == one_edge$x)
        y_idx <- which (list_nodes_not_lagged == one_edge$y)
        if (is_contextual[[x_idx]] | is_contextual[[y_idx]])
          df_edges[edge_idx, "lag"] <- ""
        }
      }
    }
  #
  # returns the tmiic structure where network summary has been flattened
  #
  tmiic_res$all.edges.summary <- df_edges
  return (tmiic_res)
  }

#-----------------------------------------------------------------------------
# tmiic_repeat_edges_over_history
#-----------------------------------------------------------------------------
# Duplicates edges found by miic over the history assuming stationarity
#
# In temporal mode, the network returned by miic contains only edges
# with at least one contemporaneous node (lag0). This function duplicates
# the edges over the history.
# i.e: assuming that we used nlayers=4 and delta_t=1, the edge X_lag0-X_lag1
# will be copied as X_lag1-X_lag2 and X_lag2-X_lag3.
#
# param: tmiic_res, the object returned by tmiic
#
# returns: a dataframe with edges completed by stationarity
#-----------------------------------------------------------------------------
tmiic_repeat_edges_over_history <- function (tmiic_res)
  {
  # Consider only edges found by miic  type = "P", "TP", "FP"
  #
  df_edges <- tmiic_res$all.edges.summary[tmiic_res$all.edges.summary$type %in% c('P', 'TP', 'FP'), ]
  if (nrow(df_edges) <= 0)
    return (df_edges)
  #
  # Precompute lag, layer and shift of each node
  #
  df_precomp <- tmiic_precompute_lags_layers_and_shifts (tmiic_res)
  list_n_layers_back <- tmiic_res$state_order$n_layers - 1
  list_nodes_not_lagged <- tmiic_res$state_order$var_names
  #
  # Duplicate the edges over all layers of history
  #
  n_edges <- nrow(df_edges)
  for (edge_idx in 1:n_edges)
    {
    node_x <- df_edges[edge_idx,"x"]
    node_y <- df_edges[edge_idx,"y"]
    node_x_pos = which (rownames (df_precomp) == node_x)
    node_y_pos = which (rownames (df_precomp) == node_y)
    #
    # If one of the variable is not lagged, the lag is not constant
    #
    sav_lag = df_precomp [node_x_pos, "lags"] - df_precomp [node_y_pos, "lags"]
    node_x_base <- df_precomp [node_x_pos, "corresp_nodes"]
    node_y_base <- df_precomp [node_y_pos, "corresp_nodes"]
    n_layers_back_x <- list_n_layers_back [[which (list_nodes_not_lagged == node_x_base)]]
    n_layers_back_y <- list_n_layers_back [[which (list_nodes_not_lagged == node_y_base)]]
    same_lag_needed = TRUE
    if (n_layers_back_x <= 0)
      same_lag_needed = FALSE
    if (n_layers_back_y <= 0)
      same_lag_needed = FALSE
    #
    # Duplication of the edge
    #
    while (TRUE)
      {
      # We shift the nodes positions using pre-computed nodes shifts
      #
      node_x_shift = df_precomp [node_x_pos, "shifts"]
      node_y_shift = df_precomp [node_y_pos, "shifts"]
      if ( (node_x_shift <= 0) & (node_y_shift <= 0) )
        break
      node_x_pos = node_x_pos + node_x_shift
      node_y_pos = node_y_pos + node_y_shift
      #
      # Ensure if both variable are lagged than we keep the same lag when duplicating
      #
      same_lag_impossible = FALSE
      if (same_lag_needed)
        {
        new_lag = df_precomp [node_x_pos, "lags"] - df_precomp [node_y_pos, "lags"]
        while (sav_lag != new_lag)
          {
          if (sav_lag < new_lag)
            {
            node_y_shift = df_precomp [node_y_pos, "shifts"]
            if (node_y_shift <= 0)
              {
              same_lag_impossible = TRUE
              break
              }
            node_y_pos = node_y_pos + node_y_shift;
            }
          else # sav_lag > new_lag
            {
            node_x_shift = df_precomp [node_x_pos, "shifts"]
            if (node_x_shift <= 0)
              {
              same_lag_impossible = TRUE
              break
              }
            node_x_pos = node_x_pos + node_x_shift;
            }
          new_lag = df_precomp [node_x_pos, "lags"] - df_precomp [node_y_pos, "lags"]
          }
        }
      if (same_lag_impossible)
        break
      #
      # Add the duplicated edge
      #
      df_edges [ nrow(df_edges)+1, ] <- df_edges [ edge_idx, ]
      df_edges [nrow(df_edges),"x"] <- rownames(df_precomp)[[node_x_pos]]
      df_edges [nrow(df_edges),"y"] <- rownames(df_precomp)[[node_y_pos]]
      }
    }
  return (df_edges)
  }

