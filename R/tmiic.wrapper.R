#*****************************************************************************
# Filename   : tmiic.wrapper.R                   Creation date: 24 march 2020
#
# Description: Data transformation of time series for miic
#
# Author     : Franck SIMON (fsimon.informaticien@wanadoo.fr)
#
# Changes history:
# - 24 march 2020 : initial version
# - 04 june 2020 : add tmiic.flatten_network
# - 15 june 2020 : add delta_tau and moving average
# - 27 july 2020 : rewrite of tmiic.transform_data_for_miic to allow variable
#                  number of timesteps between timeseries
#*****************************************************************************

#-----------------------------------------------------------------------------
# tmiic.transform_data_for_miic
#-----------------------------------------------------------------------------
# tmiic.transform_data_for_miic
#
# @description
# Reorganizes the data using the history to create lagged nodes and  
# samples in a format usable by miic
#
# @details 
# The function slices the input data according to the \emph{tau}
# and \emph{delta_tau} parameters. Data are expected to be received in a 
# dataframe with variables as columns and timeseries/timesteps as rows. 
# The timestep information must be suplied in the first column and, 
# for each timeseries, be in ascending order.
# 
# The number of variables is increased and renamed on \emph{tau} 
# / \emph{delta_tau} layers.\cr 
# i.e. with \emph{tau}=6 and \emph{delta_tau}=3 : node1, node2 => 
# node1_lag0, node2_lag0, node1_lag3, node2_lag3, node1_lag6, node2_lag6.
# 
# Every timestep (until number of timesteps - \emph{tau}) is converted
# into a sample in the lagged data. Exemple with tau=6 and delta_tau=3:
# 
# \tabular{ccccccc}{
# Timestep \tab  Node & value  \tab  Node & value  \tab => \tab Sample \tab  Node & value  \tab  Node & value \cr
#   t-6    \tab node1_val(t-6) \tab node2_val(t-6) \tab => \tab   i    \tab node1_lag6_val \tab node2_lag6_val\cr
#   t-3    \tab node1_val(t-3) \tab node2_val(t-3) \tab => \tab   i    \tab node1_lag3_val \tab node2_lag3_val\cr
#    t     \tab  node1_val(t)  \tab  node2_val(t)  \tab => \tab   i    \tab node1_lag0_val \tab node2_lag0_val\cr
#   \cr    \tab                \tab                \tab    \tab        \tab                \tab               \cr
#   t-7    \tab node1_val(t-7) \tab node2_val(t-7) \tab => \tab   i'   \tab node1_lag6_val \tab node2_lag6_val\cr
#   t-4    \tab node1_val(t-4) \tab node2_val(t-4) \tab => \tab   i'   \tab node1_lag3_val \tab node2_lag3_val\cr
#   t-1    \tab node1_val(t-1) \tab node2_val(t-1) \tab => \tab   i'   \tab node1_lag0_val \tab node2_lag0_val\cr
#   \cr    \tab                \tab                \tab    \tab        \tab                \tab               \cr
#   t-8    \tab node1_val(t-8) \tab node2_val(t-8) \tab => \tab   i"   \tab node1_lag6_val \tab node2_lag6_val\cr
#   t-5    \tab node1_val(t-5) \tab node2_val(t-5) \tab => \tab   i"   \tab node1_lag3_val \tab node2_lag3_val\cr
#   t-2    \tab node1_val(t-2) \tab node2_val(t-2) \tab => \tab   i"   \tab node1_lag0_val \tab node2_lag0_val\cr
#   \cr    \tab                \tab                \tab    \tab        \tab                \tab               \cr
#   ...    \tab .............. \tab .............. \tab => \tab ...... \tab .............. \tab ............. \cr
# }
# until number of timesteps - \emph{tau} * \emph{delta_tau} is reached. 
# The same process is applied to all input timeseries.
#
# @param input_data [a dataframe] 
# A dataframe with the time series with variables as columns and
# timeseries/timesteps as rows. The timestep information must be suplied 
# in the first column and, for each timeseries, in an ascending order.
#
# @param tau [an int > 0] A strictly positive int defining the max lag.
# Note that if \emph{delta_tau} is also supplied, \emph{tau} must be a
# multiple of \emph{delta_tau}.
# 
# @param state_order [a data frame] Optional, NULL by default. 
# A data frame giving information about how to order the various states of 
# categorical variables. This data frame will be lagged as the input data 
# on \emph{tau} / \emph{delta_tau} layers.
# 
# @param movavg [an integer] Optional, -1 by default.\cr
# When \emph{movavg} is supplied (integer > 1), a moving average 
# operation is applied to each timeseries.\cr
# 
# @param delta_tau [an integer] Optional, 1 by default.\cr
# When \emph{delta_tau} is supplied (integer > 1), the samples will be 
# constructed using 1 timestep every \emph{delta_tau} timesteps starting 
# from the last.\cr 
# i.e.: on 1000 timesteps with  \emph{tau} = 14 and \emph{delta_tau} = 7, 
# the timesteps kept for the samples conversion will be 1000, 993, 986 
# for the first sample, the next sample will use 999, 992, 985 and so on.\cr
# 
# @return a list with two elements:
# \itemize{
#  \item \emph{input_data:} the samples generated from the timeseries
#  \item \emph{state_order:} the lagged state_order 
#  }
#                                    
#-----------------------------------------------------------------------------
tmiic.transform_data_for_miic <- function (input_data, tau, state_order=NULL, 
                                           movavg=-1, delta_tau=1) {
  n_rows <- nrow(input_data)
  n_nodes <- ncol(input_data) - 1
  list_nodes <- colnames(input_data)[-1]
  #
  # Lag the nodes
  #
  list_nodes_lagged <- list()
  for (tau_idx in seq(0,tau,by=delta_tau) )
    for (node_idx in 1:n_nodes)
      list_nodes_lagged <- append (list_nodes_lagged, paste (list_nodes[[node_idx]], "_lag", tau_idx, sep="" ) )
  n_nodes_lagged <- length(list_nodes_lagged)
  #
  # Lag the state_order if supplied
  #
  if (!is.null(state_order)) {
    state_lagged <- state_order [FALSE,]
    state_col1 <- colnames(state_order)[[1]]
    for (tau_idx in seq(0,tau,by=delta_tau) ) {
      for (old_state_idx in 1:nrow(state_order)) {
        lagged_state_idx <- nrow(state_lagged) + 1
        state_lagged [lagged_state_idx,] <- state_order [old_state_idx,]
        
        node_name <- state_order [old_state_idx, state_col1]
        node_name_lagged <- paste (node_name, "_lag", tau_idx, sep="" )
        state_lagged [lagged_state_idx, state_col1] <- node_name_lagged
      }
    }
    state_order <- state_lagged
  }
  #
  # Iterate over input data to create a lagged dataset. 
  # The principle is to  identify the different timeseries using the fact
  # that the timestep information in the first column is ascending 
  # for each time series. Once a complete time series has been identified, 
  # we convert it into a lagged form
  #
  df_lagged <- data.frame (matrix(ncol=n_nodes_lagged, nrow=0), 
                           stringsAsFactors=FALSE)
  colnames (df_lagged) <- list_nodes_lagged
  previous_row_idx <- 1
  previous_timestep <- -Inf
  for (row_idx in 1:n_rows) {
    timestep = input_data[row_idx,1]
    if (timestep < previous_timestep) {
      one_timeseries <- input_data[previous_row_idx:(row_idx-1), list_nodes]
      df_ret <- tmiic.lag_one_timeseries (one_timeseries, tau, list_nodes_lagged,
                                          movavg, delta_tau)
      df_lagged <- rbind (df_lagged, df_ret)
      previous_row_idx <- row_idx
    }
    previous_timestep <- timestep
  }
  one_timeseries <- input_data[previous_row_idx:row_idx, list_nodes]
  df_ret <- tmiic.lag_one_timeseries (one_timeseries, tau, list_nodes_lagged,
                                      movavg, delta_tau)
  df_lagged <- rbind (df_lagged, df_ret)
  #
  # returns the dataframe as miic expects and state_order (if supplied) lagged
  #
  return (list (input_data=df_lagged, state_order=state_order) )
  }
  
#-----------------------------------------------------------------------------
# tmiic.lag_one_timeseries
#-----------------------------------------------------------------------------
# tmiic.lag_one_timeseries
#
# @description
# Reorganizes the data of one timeseries using the history to create lagged
# samples in a format usable by miic
#
# @details 
# The function slices the input data according to  \emph{tau} and 
# \emph{delta_tau} parameters. Data are expected to be received in a dataframe 
# with variables as columns and timesteps as rows. 
# 
# Every timestep (until number of timesteps - \emph{tau}) is converted 
# into a sample in the lagged graph. Exemple with tau=6 and delta_tau=3:
# 
# \tabular{ccccccc}{
# Timestep \tab  Node & value  \tab  Node & value  \tab => \tab Sample \tab  Node & value  \tab  Node & value \cr
#   t-6    \tab node1_val(t-6) \tab node2_val(t-6) \tab => \tab   i    \tab node1_lag6_val \tab node2_lag6_val\cr
#   t-3    \tab node1_val(t-3) \tab node2_val(t-3) \tab => \tab   i    \tab node1_lag3_val \tab node2_lag3_val\cr
#    t     \tab  node1_val(t)  \tab  node2_val(t)  \tab => \tab   i    \tab node1_lag0_val \tab node2_lag0_val\cr
#   \cr    \tab                \tab                \tab    \tab        \tab                \tab               \cr
#   t-7    \tab node1_val(t-7) \tab node2_val(t-7) \tab => \tab   i'   \tab node1_lag6_val \tab node2_lag6_val\cr
#   t-4    \tab node1_val(t-4) \tab node2_val(t-4) \tab => \tab   i'   \tab node1_lag3_val \tab node2_lag3_val\cr
#   t-1    \tab node1_val(t-1) \tab node2_val(t-1) \tab => \tab   i'   \tab node1_lag0_val \tab node2_lag0_val\cr
#   \cr    \tab                \tab                \tab    \tab        \tab                \tab               \cr
#   t-8    \tab node1_val(t-8) \tab node2_val(t-8) \tab => \tab   i"   \tab node1_lag6_val \tab node2_lag6_val\cr
#   t-5    \tab node1_val(t-5) \tab node2_val(t-5) \tab => \tab   i"   \tab node1_lag3_val \tab node2_lag3_val\cr
#   t-2    \tab node1_val(t-2) \tab node2_val(t-2) \tab => \tab   i"   \tab node1_lag0_val \tab node2_lag0_val\cr
#   \cr    \tab                \tab                \tab    \tab        \tab                \tab               \cr
#   ...    \tab .............. \tab .............. \tab => \tab ...... \tab .............. \tab ............. \cr
# }
# until number of timesteps - \emph{tau} * \emph{delta_tau} is reached. 
#
# @param df_timeseries [a dataframe] 
# A dataframe containing the timeseries with variables as columns and
# timesteps as rows. 
#
# @param tau [an int > 0] A strictly positive int defining the max lag.
# Note that if \emph{delta_tau} is also supplied, \emph{tau} must be a
# multiple of \emph{delta_tau}.
# 
# @param list_nodes_lagged [a list] 
# The list of variables lagged over \emph{tau} / \emph{delta_tau}
# 
# @param movavg [an integer] Optional, -1 by default.\cr
# When \emph{movavg} is supplied (integer > 1), a moving average 
# operation is applied to the time series.
# 
# @param delta_tau [an integer] Optional, 1 by default.\cr
# When \emph{delta_tau} is supplied (integer > 1), the samples will be 
# construted using 1 timestep every \emph{delta_tau} timesteps starting 
# from the last.\cr 
# i.e.: on 1000 timesteps with  \emph{tau} = 14 and \emph{delta_tau} = 7, 
# the timesteps kept for the samples conversion will be 1000, 993, 986 
# for the first sample, the next sample will use 999, 992, 985 and so on.\cr
# 
# 
# @return a dataframe with the samples generated form the timeseries
#-----------------------------------------------------------------------------
tmiic.lag_one_timeseries <- function (df_timeseries, tau, list_nodes_lagged, 
                                      movavg, delta_tau) {
  n_timesteps <- nrow(df_timeseries)
  n_nodes <- ncol(df_timeseries) 
  list_nodes <- colnames(df_timeseries)
  n_nodes_lagged <- length(list_nodes_lagged)
  #
  # Init df to be returned
  #
  df_lagged <- data.frame (matrix(ncol=n_nodes_lagged, nrow=0), 
                           stringsAsFactors=FALSE)
  colnames (df_lagged) <- list_nodes_lagged
  #
  # Apply moving average if requested
  #
  if (movavg > 1) {
    df_movavg <- data.frame (matrix(ncol=n_nodes, nrow=n_timesteps), 
                             stringsAsFactors=FALSE)
    colnames (df_movavg) <- list_nodes
    for (node_idx in 1:n_nodes)
      df_movavg[, node_idx] <- data.table::frollmean (df_timeseries[,node_idx], movavg)
    #
    # Do not keep the first timesteps as NA due to the moving average
    #
    df_timeseries <- df_movavg[movavg:n_timesteps,]
    n_timesteps <- nrow(df_timeseries)
  }
  #
  # Loop over timesteps to create lagged samples
  #
  for (row_idx in seq(from=n_timesteps, to=tau+1, by=-1) ) {
    new_line <- nrow(df_lagged)+1
    for (node_idx in 1:n_nodes)
      for (tau_idx in seq(0,tau,by=delta_tau) ) {
        lagged_node_idx <- (tau_idx %/% delta_tau) * n_nodes + node_idx
        df_lagged [new_line,lagged_node_idx] <- df_timeseries[row_idx-tau_idx,node_idx]
      }
  }
  return (df_lagged)
}

#-----------------------------------------------------------------------------
# tmiic.flatten_network
#-----------------------------------------------------------------------------
#' tmiic.flatten_network
#'
#' @description
#' Flattten the lagged network returned by tmiic
#'
#' @details 
#' In temporal mode, the network returned by miic contains lagged nodes 
#' (X_lag0, X_lag1, ...). This function flatten the  network depending 
#' of the \emph{flatten_mode} parameter.\cr
#' Note that only the adj_matrix and summary data frames are flatened.
#' 
#' @param miic_result [a miic object] The object returned by miic's 
#' execution.
#' 
#' @param flatten_mode [a string]. Optional, default value "normal".
#' Possible values are \emph{"normal"}, \emph{"combine"}, \emph{"unique"}, 
#' \emph{"drop"}:
#' \itemize{
#' \item When \emph{flatten_mode} = \emph{"normal"}, the default. Nodes 
#'   and edges are converted into a flattened version preserving all 
#'   the initial information.\cr
#'   i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 become respectively X->Y lag=1, 
#'   X<-Y lag=2. 
#' \item When \emph{flatten_mode} = \emph{"combine"}, one edge will be kept
#'   per couple of nodes. The log_confidence will be the highest one 
#'   of the summarized edges whilst the lag and orientation of the
#'   summarized edge will be an agregation.\cr 
#'   i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 will become X<->Y lag=1,2 with
#'   the log_confidence of X_lag1->Y_lag0 if log_confidence of 
#'   X_lag1->Y_lag0 > X_lag2<-Y_lag0.
#' \item When \emph{flatten_mode} = \emph{"unique"}, only the edges having the
#'   highest log_confidence for a couple of nodes are kept in the flattened 
#'   network. If several edges between the sames nodes have the same
#'   log_confidence, then the edge kept is the one with the minimum lag.\cr 
#'   i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 with log_confidence of 
#'   X_lag1->Y_lag0 > X_lag2<-Y_lag0 become X->Y lag=1.
#' \item When \emph{flatten_mode} = \emph{"drop"}, only the edges having the 
#'   highest log_confidence for a couple of nodes are kept in the flattened 
#'   network. If several edges between the sames nodes have the same
#'   log_confidence, then the edge kept is the one with the minimum lag.\cr
#'   i.e. :  X_lag1->Y_lag0, X_lag2<-Y_lag0 with log_confidence of 
#'   X_lag1->Y_lag0 > X_lag2<-Y_lag0 become X->Y. The lag information is 
#'   "lost" after flattening
#' }
#' Note that for all modes other than \emph{"drop"}, lag is a new column 
#' added in the dataframe. 
#' 
#' @param keep_edges_on_same_node [a boolean] Optional, TRUE by default.
#' When TRUE, the edges like X_lag0-X_lag1 are kept during flattening
#' (it becomes an X-X edge). When FALSE, only edges having different nodes 
#' are kept in the flatten network.
#' 
#' @return [a list] The returned list is the one received as input where
#' the adjacency and summary dataframes has been flattened.\cr
#' Note that the adjacency matrix does not contain any more the orientation
#' as this information make less sense after flattening. Orientations are 
#' replaced by lag(s) in the adjacency matrix.
#'     
#' @export
#-----------------------------------------------------------------------------
tmiic.flatten_network <- function (miic_result, flatten_mode="normal", 
                                   keep_edges_on_same_node=TRUE) {
  #
  # Construct the list of nodes not lagged 
  #
  list_nodes_lagged = colnames (miic_result$adj_matrix)
  list_nodes_not_lagged = list ()
  for (node_idx in 1:length (list_nodes_lagged) ) {
    #
    # As long we find "_lag0" at the end of the node name
    #
    one_node = list_nodes_lagged[[node_idx]]
    if ( endsWith(one_node, "_lag0") ) {
      #
      # Remove "_lag0" information from the nodes names and add it to list
      #
      one_node <- gsub("_lag0", "", one_node)
      list_nodes_not_lagged[[node_idx]] <- one_node
    }
    else
      break
  }
  n_nodes = length (list_nodes_not_lagged)
  #
  # First step, perform flatten_mode="normal" : 
  # - on adjacency matrix : construct a n_nodes x n_nodes with lags in each cell
  # - on summary: remove lag info from nodes names and put it into a lag column
  #
  df_adj <- array ( data="", dim=c(n_nodes ,n_nodes),
                    dimnames=list(list_nodes_not_lagged, list_nodes_not_lagged) )
  df_edges <- miic_result$all.edges.summary
  df_edges$lag <- -1
  for (edge_idx in 1:nrow(df_edges) ) {
    #
    # Get nodes of the edge without "_lagX"
    #
    one_edge <- df_edges[edge_idx,]
    #
    # Edges with type != "P", "TP", "FP" are not true edges
    #
    if ( ! (one_edge$type %in% c("P", "TP", "FP") ) ) {
      df_edges[edge_idx,]$x <- "DROP"
      df_edges[edge_idx,]$y <- "DROP"
      next
    }
    #
    # The edge is a real one, extract nodes names and orient
    #
    node_x <- one_edge$x
    node_y <- one_edge$y
    orient <- one_edge$infOrt
    #
    # Get the lag for each node of the edge
    #
    pos_lag_x <- stringr::str_locate(node_x, "_lag")
    tau_idx_x <- 0
    if ( !is.na(pos_lag_x[1]) ) {
      tau_idx_x <- stringr::str_replace(node_x, ".*_lag", "")
      tau_idx_x <- strtoi (tau_idx_x)
    }
    pos_lag_y <- stringr::str_locate(node_y, "_lag")
    tau_idx_y <- 0
    if ( !is.na(pos_lag_y[1]) ) {
      tau_idx_y <- stringr::str_replace(node_y, ".*_lag", "")
      tau_idx_y <- strtoi (tau_idx_y)
    }
    #
    # Ensure to order from oldest to newest
    #
    lag = tau_idx_x - tau_idx_y
    if (lag < 0) {
      lag <- -lag
      
      temp_var <- tau_idx_x
      tau_idx_x <- tau_idx_y
      tau_idx_y <- temp_var
      
      temp_var <- node_x
      node_x <- node_y
      node_y <- temp_var
      
      if ( (orient == 2) | (orient == -2) )
        orient = -orient
    }
    #
    # Remove "_lag" information from the nodes names
    #
    node_x <- gsub("_lag.*","",node_x)
    node_y <- gsub("_lag.*","",node_y)
    #
    # If we oblige the nodes to be different, discard edges between the same node
    #
    if ( (!keep_edges_on_same_node) & (node_x == node_y) ) {
      df_edges[edge_idx,]$x <- "DROP"
      df_edges[edge_idx,]$y <- "DROP"
      next
    }
    #
    # Edge kept, Update it
    #
    df_edges[edge_idx,]$x <- node_x
    df_edges[edge_idx,]$y <- node_y
    df_edges[edge_idx,]$infOrt <- orient
    df_edges[edge_idx,]$lag <- lag
    #
    # Update adjacency matrix
    #
    if (df_adj[node_x,node_y] == "")
      df_adj[node_x,node_y] = paste (lag, sep="")
    else
      df_adj[node_x,node_y] = paste (df_adj[node_x,node_y], ",", lag, sep="")
    df_adj[node_y,node_x] = df_adj[node_x,node_y]
  }
  #
  # Drop rows flagged 
  #
  df_edges <- df_edges[!(df_edges$x=="DROP"),]
  #
  # Order the list of lag in the adjacency matrix
  #
  for (i in 1:ncol (df_adj) ) {
    for (j in i:ncol (df_adj) ) {
      lags_str <- df_adj[i,j]
      if (nchar (lags_str) > 0)
        if ( grepl (',', lags_str, fixed = TRUE) ) {
          split_str <- strsplit (lags_str, ",")[[1]]
          split_num <- as.integer (split_str)
          split_sorted <- sort (split_num)
          df_adj[i,j] <- paste (split_sorted, collapse=',', sep="")
          df_adj[j,i] <- df_adj[i,j]
        }
    }
  }
  #
  # Step 1 "normal" node done
  #
  # If we want only one edge per couple of nodes 
  #
  if (flatten_mode != "normal") {
    # At first, we ensure X string < Y string and we create an extra
    # column containing x-y so we have an id to perform a duplicate check
    #
    for (edge_idx in 1:nrow(df_edges))  {
      one_edge <- df_edges[edge_idx,]
      if (one_edge$x < one_edge$y) {
        df_edges[edge_idx,"x"] <- one_edge$y
        df_edges[edge_idx,"y"] <- one_edge$x
        # inverse edge oriention if oriented and not bidirectional
        # not sure -4,4 orientations still exists, it was present 
        # in the old code of plot
        if ( (abs(one_edge$infOrt) == 2) | (abs(one_edge$infOrt) == 4) )
          df_edges[edge_idx, "infOrt"] <- -one_edge$infOrt
        if ( !is.na(one_edge$trueOrt) ) 
          if ( (abs(one_edge$trueOrt) == 2) | (abs(one_edge$trueOrt) == 4) )
            df_edges[edge_idx, "trueOrt"] <- -one_edge$trueOrt
      }
    }
    #
    # Keep the rows having max log_confidence when grouped by on x, y
    # If some rows between same nodes have the same log_confidence 
    # We addd a second step keeping the edges with the minimum lag 
    #
    df_rownames <- as.data.frame (rownames(df_edges), stringsAsFactors = FALSE)
    names (df_rownames) <- c("row_names")    
    df_edges_rownames <- merge (x = df_rownames, y = df_edges, by.x = "row_names", by.y = 0)
    df_group <- stats::aggregate (data.frame(log_confidence = df_edges_rownames$log_confidence, 
                                             stringsAsFactors = FALSE), 
                                  by = list(x = df_edges_rownames$x, y = df_edges_rownames$y), 
                                  max)
    df_group <- merge (x=df_group, y=df_edges_rownames, by=c("x", "y", "log_confidence") )
    df_group <- stats::aggregate (data.frame(lag = df_group$lag, stringsAsFactors = FALSE), 
        by = list(x = df_group$x, y = df_group$y, log_confidence = df_group$log_confidence), 
        min)
    df_group <- merge (x=df_group, y=df_edges_rownames, by=c("x", "y", "lag", "log_confidence") )
    df_group <- stats::aggregate (data.frame(row_names = df_group$row_names, stringsAsFactors = FALSE), 
        by = list(x = df_group$x, y = df_group$y, log_confidence = df_group$log_confidence, lag = df_group$lag), 
        min)
    df_group <- merge (x=df_group, y=df_edges_rownames, by=c("x", "y", "lag", "log_confidence", "row_names") )
    #
    # If the flatten_mode is "unique" or "drop", nothing more to do for now
    #
    # If the flatten_mode is "combine", combine orientations
    # 
    if (flatten_mode == "combine") {
      #
      # Update the orientations (in case the grouped edges have different 
      # orientations) and put in the lag column the list of lags
      #
      for (edge_idx in 1:nrow(df_group) ) {
        one_edge <- df_group[edge_idx,]
        node_x <- one_edge$x
        node_y <- one_edge$y
        row_name <- one_edge$row_names
        #
        # Update the list of lags (value has been computed in lag matrix)
        #
        df_group[edge_idx,]$lag <- df_adj[node_x, node_y]
        #
        # Select the other edges between same nodes for orientation update
        #
        cond_for_orient <- ( (df_edges_rownames[["x"]] == node_x) 
                           & (df_edges_rownames[["y"]] == node_y) 
                           & (df_edges_rownames[["row_names"]] != row_name) )
        #
        # If edge was unique (no different lag between the nodes), nothing to do
        #
        if (sum(cond_for_orient) == 0) 
          next
        #
        # Edge was having two or more lags, update proba and orientation
        #
        head_proba <- 0
        tail_proba <- 0
        if ( (abs(one_edge$infOrt) >= 2) & (!is.na(one_edge$proba)) ) {
          edge_probas <- strsplit (one_edge$proba, ';')[[1]]
          if (one_edge$infOrt > 0) {
            head_proba <- edge_probas[[1]]
            tail_proba <- edge_probas[[2]]
          }
          else {
            head_proba <- edge_probas[[2]]
            tail_proba <- edge_probas[[1]]
          }
        }
        
        df_other <- df_edges_rownames[cond_for_orient,]
        for (other_idx in 1:nrow(df_other) ) {
          other_edge <- df_other[other_idx,]
          #
          # Update infOrt and trueOrt
          #
          for (col_to_update in c("infOrt", "trueOrt")) {
            other_edge_orient <- other_edge[[col_to_update]]
            group_edge_orient <-  df_group[edge_idx,][[col_to_update]]
            #
            # If other edge is na has no orient, nothing to do
            #
            if ( (is.na(other_edge_orient)) | (other_edge_orient == 1) ) {
              next            
            }
            #
            # If grouped edge is na or has no orient, we can always update
            #
            if (  (is.na(group_edge_orient)) | (group_edge_orient == 1) ) {
              df_group[edge_idx,][[col_to_update]] <- other_edge_orient
              next            
            }
            #
            # The 2 edges have an orientation (unidirectional or birectional)
            #
            # If the other edge is bidirectional or opposite to the grouped one 
            # -> update to bidirectional
            #
            if ( (other_edge_orient == 6) | (other_edge_orient == -group_edge_orient) ) {
              df_group[edge_idx,][[col_to_update]] <- 6
              next            
            }
          }
          #
          # Keep max proba of tail/queue
          #
          other_edge_head_proba <- 0
          other_edge_tail_proba <- 0
          if ( (abs(other_edge$infOrt) >= 2) & (!is.na(other_edge$proba)) ) {
            other_edge_probas <- strsplit (other_edge$proba, ';')[[1]]
            if (other_edge$infOrt > 0) {
              other_edge_head_proba <- other_edge_probas[[1]]
              other_edge_tail_proba <- other_edge_probas[[2]]
            }
            else {
              other_edge_head_proba <- other_edge_probas[[2]]
              other_edge_tail_proba <- other_edge_probas[[1]]
            }
          }
          if (other_edge_head_proba > head_proba)
            head_proba <- other_edge_head_proba
          if (other_edge_tail_proba > tail_proba)
            tail_proba <- other_edge_tail_proba
        }
        #
        # Update the grouped edge with the new probas
        #
        if (df_group[edge_idx,]$infOrt >= 0)
          df_group[edge_idx,]$proba <- paste(head_proba, tail_proba, sep=";")
        else
          df_group[edge_idx,]$proba <- paste(tail_proba, head_proba, sep=";")
      }
    }
    df_edges <- within (df_group, rm(row_names))
  }
  #
  # If we do not want to keep info about lag at all
  #
  if (flatten_mode == "drop")
    df_edges$lag <- NULL
  #
  # returns the list of dataframes where network has been flattened
  #
  miic_result$adj_matrix <- df_adj
  miic_result$all.edges.summary <- df_edges
  miic_result$tmiic_specific[["is_lagged"]] = FALSE
  return (miic_result)
}

