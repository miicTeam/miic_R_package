#*******************************************************************************
# Filename   : tmiic.wrapper.R                     Creation date: 24 march 2020
#
# Description: Data transformation of time series for miic
#
# Author     : Franck SIMON
#*******************************************************************************

#===============================================================================
# FUNCTIONS
#===============================================================================
# tmiic_extract_trajectories
#-------------------------------------------------------------------------------
# Extract the trajectories from a data frame and return them in a list
# - input_data: a data frame with the time steps in the 1st column
#   A new trajectory is identified when time step < previous time step
# - check: optional, default=T. Emit warnings when:
#   * there is a gap between 2 consecutive time steps
#   * the time step value is not incremented  between 2 consecutive rows
#   * the 1st time step of a trajectory is not 1
# Returns:
# - a list: the list of trajectories
#   Note that the time step information in each trajectory is renumbered from 1
#   to number of time steps of the trajectory (so no gap, no unchanged time step)
#-------------------------------------------------------------------------------
tmiic_extract_trajectories <- function (input_data, check=T)
  {
  timesteps <- input_data[, 1]
  if ( any ( is.na (timesteps) ) )
    miic_error ("trajectories check", "the time steps column (column 1) contains NA(s).")
  if ( ! all (is.numeric (timesteps)) )
    miic_error ("trajectories check", "the time steps column (column 1) is not integer.")
  if ( ! all (round (timesteps, 0) == timesteps) )
    miic_error ("trajectories check", "the time steps column (column 1) is not integer.")
  timesteps <- as.integer(timesteps)
  timesteps_next <- c (timesteps[2:length(timesteps)], 0)
  breaks <- which (timesteps_next < timesteps)

  list_traj <- list()
  row_prev <- 1
  for ( i in 1:length (breaks) )
    {
    row_new <- breaks[[i]]
    list_traj[[i]] <- input_data[row_prev:row_new,,F]
    row_prev <- row_new + 1
    }

  if (check)
    {
    no_inc <- which (timesteps_next == timesteps)
    if (length (no_inc) > 0)
      miic_warning ("check trajectories", "time steps unchanged at ",
                    length (no_inc), " position(s).")
    gaps <- which (timesteps_next > timesteps + 1)
    if (length (gaps) > 0)
      miic_warning ("check trajectories", "gap in time steps at ",
                    length (gaps), " position(s).")
    wrong_starts <- which ( unlist (lapply (list_traj,
                      FUN=function (x) { return (x[1,1] != 1) } ) ) )
    if (length (wrong_starts) > 0)
      miic_warning ("check trajectories", length (wrong_starts),
        " trajectorie(s) don't start with 1 as first time step.")
    max_nb_ts <- max (unlist (lapply (list_traj, FUN=nrow) ) )
    if (max_nb_ts == 1)
      miic_error ("trajectories check",
                  "all trajectories have only 1 time step.")
    }
  #
  # Renum time steps
  #
  for ( i in 1:length (list_traj) )
    list_traj[[i]][,1] <- 1:nrow (list_traj[[i]])
  return (list_traj)
  }

#-------------------------------------------------------------------------------
# tmiic_group_trajectories
#-------------------------------------------------------------------------------
# Merge a list of trajectories into a data frame
# - list_traj: the list of trajectories
# Returns:
# - a data frame: data frame with all the trajectories
#-------------------------------------------------------------------------------
tmiic_group_trajectories <- function (list_traj)
  {
  # Pre-allocate the data frame with the same structure as trajectories
  # and the same number of rows as all the trajectories
  # VOIR data <- data[1:n_row_tot,]
  #
  df <- list_traj[[1]][FALSE,]
  n_row_tot <- sum (unlist (lapply(list_traj, nrow)))
  df <- df[seq_len(n_row_tot), , drop=F]
  rownames(df) <- NULL

  row_idx <- 1
  for (i in 1:length(list_traj) )
    if (nrow(list_traj[[i]]) > 0)
      {
      df[row_idx:(row_idx-1+nrow(list_traj[[i]])),] <- list_traj[[i]]
      row_idx <- row_idx + nrow(list_traj[[i]])
      }
  return (df)
  }

#-------------------------------------------------------------------------------
# tmiic_precompute_lags_layers_and_shifts
#-------------------------------------------------------------------------------
# Utility function to precompute lags, layers and shifts of nodes in the
# lagged network
#
# params: tmiic_obj [a tmiic object] The object returned by miic's
# execution in temporal mode.
#
# returns: a dataframe with lagged nodes as row name and 3 columns:
#  - lags: the lag of each lagged node
#  - corresp_nodes: the corresponding non lagged node
#  - shifts: the shift to apply to find the next lagged node
#-----------------------------------------------------------------------------
tmiic_precompute_lags_layers_and_shifts <- function (tmiic_obj)
  {
  list_nodes_not_lagged = tmiic_obj$state_order$var_names
  is_contextual = tmiic_obj$state_order$is_contextual
  n_nodes_not_lagged = length (list_nodes_not_lagged)
  list_n_layers_back <- tmiic_obj$state_order$n_layers - 1
  list_delta_t <- tmiic_obj$state_order$delta_t
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
      n_layers_back_of_var <- list_n_layers_back[[node_idx]]
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
    if (df[idx,"ort_inferred"] == -2)
      df[idx, c("x","y","lag")] <- c (df[idx,"y"], df[idx,"x"],
                                      -as.integer (df[idx,"lag"]) )

    if ( (df[idx,"ort_inferred"] == 6) & (as.integer (df[idx,"lag"]) != 0) )
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
# Utility function to combine edge probabilities when flattening the network.
# Depending on the combined edges orientation, chooses the appropriate max, min
# or mean probabilities to compute the combined edge probabilities
#
# params:
# - df: the data frame with the edges to combine
# - comb_orient: integer, the orientation of the combined edge
#-----------------------------------------------------------------------------
tmiic_combine_probas <- function (df, comb_orient)
  {
  df <- df[ (!is.na (df[, "p_y2x"]) )
          & (!is.na (df[, "p_x2y"]) ), , drop=F]
  if (nrow (df) <= 0)
    return ( c(NA_real_, NA_real_) )
  #
  # We set probas like if we have node X <= node Y
  #
  for ( idx in 1:nrow(df) )
    if (df[idx,"x"] > df[idx,"y"])
      {
      temp <- df[idx, "p_y2x"]
      df[idx, "p_y2x"] <- df[idx, "p_x2y"]
      df[idx, "p_x2y"] <- temp
      }
  #
  # Depending on the pre-computed combined orientation, keep max/min/avg
  #
  if (comb_orient == 6)
    return (c (max(df$p_y2x), max(df$p_x2y) ) )
  if (comb_orient == 2)
    return (c (min(df$p_y2x), max(df$p_x2y) ) )
  if (comb_orient == -2)
    return (c (max(df$p_y2x), min(df$p_x2y) ) )
  return (c (mean(df$p_y2x), mean(df$p_x2y) ) )
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
# - tmiic_obj: a tmiic object, returned by tmiic
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
tmiic_flatten_network <- function (tmiic_obj, flatten_mode="compact",
                                   keep_edges_on_same_node=TRUE)
  {
  # Reduce size of adj_matrix to non lagged nodes
  # (we don't care about content as it is not used for plotting)
  #
  list_nodes <- tmiic_obj$state_order$var_names
  tmiic_obj$adj_matrix <- matrix(NA, nrow=0, ncol=length (list_nodes))
  colnames(tmiic_obj$adj_matrix) <- list_nodes
  #
  # Keep only edges found by miic
  #
  df_edges <- tmiic_obj$summary[tmiic_obj$summary$type %in% c('P', 'TP', 'FP'), ]
  if (nrow(df_edges) <= 0)
    {
    if (flatten_mode != "drop")
      df_edges$lag = numeric(0)
    tmiic_obj$summary <- df_edges
    return (tmiic_obj)
    }
  #
  # Precompute lag and layer of each node
  #
  df_precomputed <- tmiic_precompute_lags_layers_and_shifts (tmiic_obj)
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
      if (abs (one_edge$ort_inferred) == 2)
         df_edges [edge_idx,"ort_inferred"] <- -one_edge$ort_inferred
      if ( !is.na (one_edge$ort_ground_truth ) )
        if (abs (one_edge$ort_ground_truth ) == 2)
          df_edges [edge_idx,"ort_ground_truth"] <- -one_edge$ort_ground_truth
      if (  (!is.na(one_edge$p_y2x))
         && (!is.na(one_edge$p_x2y)) )
        {
        temp <- one_edge$p_y2x
        df_edges[edge_idx, "p_y2x"] <- one_edge$p_x2y
        df_edges[edge_idx, "p_x2y"] <- temp
        }
      }
    }
  df_edges <- transform ( df_edges, lag = as.integer (lag) )
  #
  # Exclude self loops if requested
  #
  if (!keep_edges_on_same_node)
    df_edges <- df_edges[df_edges$x != df_edges$y, , drop=F]
  if (nrow(df_edges) <= 0)
    {
    if (flatten_mode == "drop")
      df_edges$lag <- NULL
    tmiic_obj$summary <- df_edges
    return (tmiic_obj)
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
    df_group <- df_edges[FALSE, , drop=F]
    for ( xy_idx in 1:nrow(df_xy) )
      {
      ref_x <- df_xy[xy_idx,"x"]
      ref_y <- df_xy[xy_idx,"y"]
      cond_same_edges = ( ( (df_edges[["x"]] == ref_x) & (df_edges[["y"]] == ref_y) )
                        | ( (df_edges[["x"]] == ref_y) & (df_edges[["y"]] == ref_x) ) )
      df_same <- df_edges[cond_same_edges, , drop=F]

      if (nrow (df_same) > 1)
        {
        if (flatten_mode == "combine")
          {
          # Combine lag, orient and proba
          #
          df_same$new_lag <-  tmiic_combine_lag (df_same)
          comb_ort_inferred <- tmiic_combine_orient (df_same, "ort_inferred")
          tmp_ret <- tmiic_combine_probas (df_same, comb_ort_inferred)
          df_same$p_y2x <- tmp_ret[[1]]
          df_same$p_x2y <- tmp_ret[[2]]
          df_same$ort_ground_truth <- tmiic_combine_orient (df_same, "ort_ground_truth")
          df_same$ort_inferred <- comb_ort_inferred
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
        df_same <- df_same[ (df_same[["info_shifted"]] == max_info), , drop=F]
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
    is_contextual <- tmiic_obj$state_order$is_contextual
    if (!is.null(is_contextual))
      {
      list_nodes_not_lagged = tmiic_obj$state_order$var_names
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
  tmiic_obj$summary <- df_edges
  return (tmiic_obj)
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
# param: tmiic_obj, the object returned by tmiic
#
# returns: a dataframe with edges completed by stationarity
#-----------------------------------------------------------------------------
#
tmiic_repeat_edges_over_history <- function (tmiic_obj)
  {
  # Consider only edges found by miic  type = "P", "TP", "FP"
  #
  df_edges <- tmiic_obj$summary[tmiic_obj$summary$type %in% c('P', 'TP', 'FP'), , drop=F]
  if (nrow(df_edges) <= 0)
    return (df_edges)
  #
  # Precompute lag, layer and shift of each node
  #
  df_precomp <- tmiic_precompute_lags_layers_and_shifts (tmiic_obj)
  list_n_layers_back <- tmiic_obj$state_order$n_layers - 1
  list_nodes_not_lagged <- tmiic_obj$state_order$var_names
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

