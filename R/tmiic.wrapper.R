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
# - 27 july 2020 : rewrite of tmiic.lag_inputs to allow variable
#                  number of timesteps between timeseries
# - 06 oct  2020 : add tmiic.repeat_edges_over_history to duplicate the edges 
#                  over history
# - 07 jan  2021 : replace tmiic.transform_data_for_miic by tmiic.lag_inputs:
#                  allow different tau, delta_tau or movavg per variable
#                  add support for contextual variables (not lagged)
#                  transfer input data lagging into C++ function 
#                  lag state_order, true_edges and black_box
# - 10 fev  2021 : add support in flattening and repeat of edges over history 
#                  of different tau, delta_tau per variable
#*****************************************************************************

#-----------------------------------------------------------------------------
# tmiic.lag_inputs
#-----------------------------------------------------------------------------
# tmiic.lag_inputs
#
# @description
# Reorganizes the inputs in a format usable by miic: input data are lagged
# using the history to create lagged variables and the optional extra inputs 
# state_order, true_edges and black_box are modified to match the updated 
# list of variables
#
# @details 
# The function slices the input data according to the \emph{tau}
# and \emph{delta_tau} parameters. Data are expected to be received in a 
# dataframe with variables as columns and timeseries/timesteps as rows. 
# The timestep information must be supplied in the first column and, 
# for each timeseries, be consecutive (incremented by 1) and in ascending 
# order.
# 
# The number of variables is increased and renamed on \emph{tau} 
# layers by \emph{delta_tau} steps.\cr 
# i.e. with \emph{tau}=2 and \emph{delta_tau}=3 : var1, var2 => 
# var1_lag0, var2_lag0, var1_lag3, var2_lag3, var1_lag6, var2_lag6.
# 
# Every timestep (until number of timesteps - \emph{tau} * \emph{delta_tau}) 
# is converted into a sample in the lagged data. 
# Exemple with tau=2 and delta_tau=3:
# \tabular{ccccccc}{
# Timestep \tab  Var & value  \tab  Var & value  \tab => \tab Sample \tab  Var & value  \tab  Var & value \cr
#   t-6    \tab Var1_val(t-6) \tab Var2_val(t-6) \tab => \tab   i    \tab Var1_lag6_val \tab Var2_lag6_val\cr
#   t-3    \tab Var1_val(t-3) \tab Var2_val(t-3) \tab => \tab   i    \tab Var1_lag3_val \tab Var2_lag3_val\cr
#    t     \tab  Var1_val(t)  \tab  Var2_val(t)  \tab => \tab   i    \tab Var1_lag0_val \tab Var2_lag0_val\cr
#   \cr    \tab               \tab               \tab    \tab        \tab               \tab              \cr
#   t-7    \tab Var1_val(t-7) \tab Var2_val(t-7) \tab => \tab   i'   \tab Var1_lag6_val \tab Var2_lag6_val\cr
#   t-4    \tab Var1_val(t-4) \tab Var2_val(t-4) \tab => \tab   i'   \tab Var1_lag3_val \tab Var2_lag3_val\cr
#   t-1    \tab Var1_val(t-1) \tab Var2_val(t-1) \tab => \tab   i'   \tab Var1_lag0_val \tab Var2_lag0_val\cr
#   \cr    \tab               \tab               \tab    \tab        \tab               \tab              \cr
#   t-8    \tab Var1_val(t-8) \tab Var2_val(t-8) \tab => \tab   i"   \tab Var1_lag6_val \tab Var2_lag6_val\cr
#   t-5    \tab Var1_val(t-5) \tab Var2_val(t-5) \tab => \tab   i"   \tab Var1_lag3_val \tab Var2_lag3_val\cr
#   t-2    \tab Var1_val(t-2) \tab Var2_val(t-2) \tab => \tab   i"   \tab Var1_lag0_val \tab Var2_lag0_val\cr
#   \cr    \tab               \tab               \tab    \tab        \tab               \tab              \cr
#   ...    \tab ............. \tab ............. \tab => \tab ...... \tab ............. \tab ............ \cr
# }
# until number of timesteps - \emph{tau} * \emph{delta_tau} is reached. 
# The same process is applied to all input timeseries.\cr
# \cr
# Note that the lagging can be different for each input variable
# if different values of tau or delta_tau are supplied and some 
# variables can be not lagged at all like contextual ones.
#
# if supplied, optional inputs are modified to mach modified input
# data:
# \itemize{
# \item \emph{state_order} has extra rows added for lagged variables:
#  Var_names var_type levels_increasing_order
#  var1         0	              0
#  var2         1	              NA
#  will become with tau=2 and delta_tau=50:
#  Var_names   var_type levels_increasing_order
#  var1_lag0      0	              0
#  var2_lag0      1	              NA
#  var1_lag50     0	              0
#  var2_lag50     1	              NA
#  var1_lag100    0	              0
#  var2_lag100    1	              NA
# \item \emph{true_edges} is modified to match the name of lagged variables:
#  var1 var2 0
#  var3 var4 2
#  will become:
#  var1_lag0 var2_lag0
#  var3_lag2 var4_lag0
# \item \emph{black_box} is modified to match the name of lagged variables:
#  var1 var2 0
#  var3 var4 2
#  will become:
#  var1_lag0 var2_lag0
#  var3_lag2 var4_lag0
#  }
#
# @param input_data [a dataframe] 
# A dataframe of the time series with variables as columns and
# timeseries/timesteps as rows. The timestep information must be supplied in 
# the first column and, for each timeseries, be consecutive (increment of 1)
# and in an ascending order.
#
# @param tau [an int > 0 or a list of int] Optional, NULL by default.\cr
# Tau(s) define(s) the number of layers that will be considered for the variables. 
# The layers will be distant of \emph{delta_tau} timesteps.\cr
# \itemize{
# \item When NULL is supplied, the state_order parameter must contains a tau
#       column to specify the number of layers of each variable. 
# \item When an integer is supplied, it must be strictly positive and all the 
#       variables will be lagged using this value as the number of layers.
# \item When a list is supplied, the size of the list must match the number of 
#       variables (so excluding the timesteps column) and each variable will 
#       be lagged by the corresponding number of layers. 
#       Variables associated with a \emph{tau} values < 1  will not be lagged.
#  }
# 
# @param delta_tau [an integer or a list of int] Optional, 1 by default.\cr
# Delta_tau(s) define(s) the number of timesteps between each layer.\cr
# i.e.: on 1000 timesteps with  \emph{tau} = 2 and \emph{delta_tau} = 7, 
# the timesteps kept for the samples conversion will be 1000, 993, 986 
# for the first sample, the next sample will use 999, 992, 985 and so on.\cr
# \itemize{
# \item When an integer is supplied (integer > 0). All the variables will be 
#       lagged using this value as the number of timesteps between layers.
# \item When a list is supplied, the size of the list must match the number of 
#       variables (so excluding the timesteps column) and each variable will 
#       be lagged with corresponding \emph{delta_tau} timesteps between layers. 
#       If a variable is not lagged ( \emph{tau} < 1), the \emph{delta_tau}
#       is ignored.
# 
# @param movavg [an integer or a list of int] Optional, -1 by default.\cr
# \itemize{
# \item When an integer is supplied (integer > 1), a moving average 
#       operation is applied to each integer and numeric variable.
# \item When a list is supplied, the size of the list must match the number of 
#       variables. A moving average will be applied to each variable having
#       a movavg value greater than 1
#  }
# 
# @param state_order [a data frame] Optional, NULL by default. 
# A data frame giving extra information about how to process the data.
# If the dataframe contains columns named 'tau', 'delta_tau', 'movavg' or
# 'is_contextual', it will modifies the lagging process:
# If some 'tau', 'delta_tau' or 'movavg' columns are present, the column 
# values from the state_order will overwrite the function parameters.
# If 'is_contextual' column is present, the information will be used
# to prevent contextual variables from being lagged.
# The content of the state_order will be lagged by the 
# \emph{tmiic.lag_inputs} function to match the lagged input data
# (see details section for an exemple).
# 
# @param true_edges [a data frame] Optional, NULL by default. If supplied, 
# the dataframe must have 3 columns: origin node, destination node and lag.
# The content of the true_edges dataframe will be lagged to match the lagged 
# input data (see details section for an exemple).
# 
# @param black_box [a data frame] Optional, NULL by default. If supplied, 
# the dataframe must have 3 columns : origin node, destination node and lag.
# the content of the black_boxs dataframe will be lagged to match the lagged
# input data (see details section for an exemple).
# 
# @param keep_max_data [a boolean] Optional, TRUE by default. 
# Lagging, and if used, moving average will produce NAs during the process.
# When keep_max_data is TRUE, rows contained cells with such introduced NAs
# will be kept, whilst these rows will be dropped when keep_max_data is FALSE.
# 
# @return a list with four elements:
# \itemize{
#  \item \emph{input_data:} the samples generated from the timeseries
#  \item \emph{state_order:} the lagged state_order 
#  \item \emph{true_edges:} the lagged true_edges 
#  \item \emph{black_box:} the lagged black_box 
#  }
#                                    
#' @useDynLib miic
#-----------------------------------------------------------------------------
tmiic.lag_inputs <- function (input_data, tau=NULL, delta_tau=1, movavg=-1, 
                              state_order=NULL, true_edges=NULL, 
                              black_box=NULL, keep_max_data=TRUE) {
  list_vars <- colnames (input_data)[-1]
  n_vars <- length (list_vars)
  #
  # Associate each variable with its tau, delta_tau and movavg
  #
  is_contextual <- NULL
  if ( !is.null (state_order) ) {
    factor_cols <- sapply (state_order, is.factor)
    state_order[factor_cols] <- lapply (state_order[factor_cols], as.character)
    
    if ( "tau" %in% colnames(state_order) ) 
      tau <- sapply (as.numeric (state_order$tau), 
                     function (x) ifelse (is.na(x), -1, x) )
    if ( "delta_tau" %in% colnames(state_order) ) 
      delta_tau <- sapply (as.numeric (state_order$delta_tau), 
                           function (x) ifelse (is.na(x), 1, x) )
    if ( "movavg" %in% colnames(state_order) ) 
      movavg <- sapply (as.numeric (state_order$movavg), 
                        function (x) ifelse (is.na(x), -1, x) )
    if ( "is_contextual" %in% colnames(state_order) )
      is_contextual <- sapply (as.numeric (state_order$is_contextual),
                               function (x) ifelse (is.na(x), 0, x) )
  }
  if (is.null(tau))
    stop ("Error: Tau not supplied as parameter, nor in the state_order")
  
  if (length (tau) == 1) 
    tau <- rep (tau, n_vars)
  if (length (delta_tau) == 1)
    delta_tau <- rep (delta_tau, n_vars)
  if (length (movavg) == 1) {
    movavg <- rep (movavg, n_vars)
    # if a global movavg is supplied, apply it only to int and numeric vars
    list_classes <- lapply (input_data[,-1], class)
    list_not_int_or_num <- which (  (list_classes != "integer") 
                                  & (list_classes != "numeric") )
    movavg [list_not_int_or_num] <- -1
  }
  if (length (tau) != n_vars)
    stop (paste0 ("Error: Tau list has ", length(tau), 
                  " items while data has ", n_vars, " variables") )
  if (length (delta_tau) != n_vars)
    stop (paste0 ("Error: Delta_tau list has ", length(delta_tau), 
                  " items while data has ", n_vars, " variables") )
  if (length (movavg) != n_vars)
    stop (paste0 ("Error: movavg list has ", length(movavg), 
                  " items while data has ", n_vars, " variables") )
  
  if (! is.null (is_contextual) )
    # Contextual variables are not lagged
    tau [which (is_contextual == TRUE)] <- -1
  tau [which (tau < 1)] <- -1
  delta_tau [which (tau < 1)] <- 0
  delta_tau [which (delta_tau <= -1)] <- 0
  #
  # Lag the state_order if supplied
  #
  if (!is.null (state_order)) {
    if (nrow(state_order) > 0) {
      #
      # Associate state_order row lines to var_names
      #
      row_of_var <- list()
      for (var_idx in 1:n_vars) {
        var_name <- list_vars[[var_idx]]
        row_in_state_order <- which (state_order$var_names == var_name)
        if (! row_in_state_order)
          stop (paste0 ("Error: Variable ", var_name, " not found in the state order file" ) )
        row_of_var[[var_name]] <- row_in_state_order
      }
      #
      # Put lag0 and not lagged variable first
      # 
      state_lagged <- state_order[FALSE,]
      state_lagged_nrow <- 0
      for (var_idx in 1:n_vars) {
        var_name <- list_vars[[var_idx]]
        state_lagged_nrow <- state_lagged_nrow + 1
        state_lagged [state_lagged_nrow,] <- state_order [row_of_var[[var_name]],]
        if (tau[[var_idx]] > 0)
          state_lagged [state_lagged_nrow, "var_names"] <- paste0 (var_name, "_lag0")
      }
      #
      # Duplicate rows for lagged variables
      # 
      tau_max <- max (tau)
      for (tau_idx in 1:tau_max) {
        for (var_idx in 1:n_vars) {
          tau_of_var <- tau[[var_idx]];
          if (tau_idx <= tau_of_var) {
            var_name <- list_vars[[var_idx]]
            state_lagged_nrow <- state_lagged_nrow + 1
            state_lagged [state_lagged_nrow,] <- state_order [row_of_var[[var_name]],]
            lag <- tau_idx * delta_tau[[var_idx]];
            state_lagged [state_lagged_nrow, "var_names"] <- paste0 (var_name, "_lag", lag)
          }
        }
      }
      state_order <- state_lagged
    }
  }
  #
  # Lag the true_edges if supplied: input X, Y, lagZ => ouput X_lagZ, Y
  #
  if (!is.null (true_edges)) {
    if (nrow (true_edges) > 0) {
      factor_cols <- sapply (true_edges, is.factor)
      true_edges[factor_cols] <- lapply (true_edges[factor_cols], as.character)
      
      for (i in 1:nrow (true_edges)) {
        orig_node <- true_edges [i, 1]
        orig_node_idx <- which (list_vars == orig_node)
        if (length(orig_node_idx) <= 0)
          stop (paste0 ("Error: Variable ", orig_node, " from true edges not found in input data" ) )
        if (tau [[orig_node_idx]] >= 1) 
          true_edges [i, 1] <- paste0 (true_edges [i, 1], "_lag", true_edges [i, 3])
        
        dest_node <- true_edges [i, 2]
        dest_node_idx <- which (list_vars == dest_node)
        if (length(dest_node_idx) <= 0)
          stop (paste0 ("Error: Variable ", dest_node, " from true edges not found in input data" ) )
        if (tau [[dest_node_idx]] >= 1) 
          true_edges [i, 2] <- paste0 (true_edges [i, 2], "_lag0")
      }
      true_edges <- true_edges[,-3]
    }
  }
  #
  # Lag the black_box if supplied: input is X, Y, lagZ => ouput X_lagZ, Y
  #
  if (!is.null(black_box)) {
    if (nrow(black_box) > 0) {
      factor_cols <- sapply (black_box, is.factor)
      black_box[factor_cols] <- lapply (black_box[factor_cols], as.character)
      
      for (i in 1:nrow (black_box)) {
        orig_node <- black_box [i, 1]
        orig_node_idx <- which (list_vars == orig_node)
        if (length(orig_node_idx) <= 0)
          stop (paste0 ("Error: Variable ", orig_node, " from black box not found in input data" ) )
        if (tau [[orig_node_idx]] >= 1) 
          black_box [i, 1] <- paste0 (black_box [i, 1], "_lag", black_box [i, 3])
        
        dest_node <- black_box [i, 2]
        dest_node_idx <- which (list_vars == dest_node)
        if (length(dest_node_idx) <= 0)
          stop (paste0 ("Error: Variable ", dest_node, " from black box not found in input data" ) )
        if (tau [[dest_node_idx]] >= 1) 
          black_box [i, 2] <- paste0 (black_box [i, 2], "_lag0")
      }
      black_box <- black_box[,-3]
    }
  }
  #
  # Lag the dataset
  #
  ret_data <- data.frame ()
  arg_list <- list ("tau"=tau, "delta_tau"=delta_tau, "movavg"=movavg, 
                    "keep_max_data"=keep_max_data)
  cpp_input <- list ("input_data"=input_data)
  if (base::requireNamespace ( "Rcpp", quietly=TRUE) ) {
    rescpp <- lagData (cpp_input, arg_list)
    ret_data <- as.data.frame(rescpp$output_data)
  }
  #
  # Check that variables in true edges and black box are valid after lagging
  #
  error_flag <- FALSE
  cols_lagged <- colnames (ret_data)
  if (!is.null (true_edges))
    if (nrow (true_edges) > 0)
        for (i in 1:nrow (true_edges)) {
          if ( ! (true_edges[[i,1]] %in% cols_lagged) ) {
            print (paste0 ("Error: variable ", true_edges[[i,1]], " in true_edge does not exist in lagged varaiables") )
            error_flag <- TRUE
          }
          if ( ! (true_edges[[i,2]] %in% cols_lagged) ) {
            print (paste0 ("Error: variable ", true_edges[[i,2]], " in true_edge does not exist in lagged varaiables") )
            error_flag <- TRUE
          }
        }
  if (!is.null (black_box))
    if (nrow (black_box) > 0)
        for (i in 1:nrow (black_box)) {
          if ( ! (black_box[[i,1]] %in% cols_lagged) ) {
            print (paste0 ("Error: variable ", black_box[[i,1]], " in black_box does not exist in lagged varaiables") )
            error_flag <- TRUE
          }
          if ( ! (black_box[[i,2]] %in% cols_lagged) ) {
            print (paste0 ("Error: variable ", black_box[[i,2]], " in black_box does not exist in lagged varaiables") )
            error_flag <- TRUE
          }
        }
  if (error_flag)
    stop ("Error: Mismatch between input data and complementary parameters")
  
  ret <- list (input_data=ret_data, tau=tau, delta_tau=delta_tau, 
               state_order=state_order, true_edges=true_edges, 
               black_box=black_box)
  return (ret)
}

#-----------------------------------------------------------------------------
# tmiic.precompute_lags_layers_and_shifts
#-----------------------------------------------------------------------------
# @description
# Utility function to precompute lags, layers and shifts of nodes in the 
# lagged network
#
# @param tmiic.res [a tmiic object] The object returned by miic's 
# execution in temporal mode.
#
# @return a dataframe with lagged nodes as row name and 3 columns: 
# \itemize{
#  \item \emph{lags:} the lag of each lagged node
#  \item \emph{corresp_nodes:} the corresponding non lagged node  
#  \item \emph{shifts:} the shift to apply to find the next lagged node
#  }
#-----------------------------------------------------------------------------
tmiic.precompute_lags_layers_and_shifts <- function (tmiic.res) {
  list_nodes_not_lagged = tmiic.res$tmiic_specific[["nodes_not_lagged"]]
  n_nodes_not_lagged = length (list_nodes_not_lagged)
  list_taus <- tmiic.res$tmiic_specific[["tau"]]
  list_delta_taus <- tmiic.res$tmiic_specific[["delta_tau"]]
  #
  # Identify lag and layer of each node
  #
  list_lags <- rep(0, n_nodes_not_lagged)
  list_nodes_lagged <- c()
  list_corresp_nodes <- c()
  i = 1
  for (node_idx in 1:n_nodes_not_lagged) {
    node_name <- list_nodes_not_lagged[[node_idx]]
    list_corresp_nodes[[i]] <- node_name
    
    if (list_taus[[node_idx]] >= 1)
      node_name <- paste0 (node_name, "_lag0")
    list_nodes_lagged [[i]] <- node_name
    i <- i + 1
  }
  
  tau_max <- max (list_taus)
  for (tau_idx in 1:tau_max) {
    for (node_idx in 1:n_nodes_not_lagged) {
      tau_of_var <- list_taus[[node_idx]];
      if (tau_idx <= tau_of_var) {
        node_name <- list_nodes_not_lagged[[node_idx]]
        list_corresp_nodes[[i]] <- node_name
        
        lag <- tau_idx * list_delta_taus[[node_idx]];
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
  for (tau_idx in 1:(tau_max+1) ) {
    for (node_idx in 1:n_nodes_not_lagged) {
      tau_of_var <- list_taus[[node_idx]];
      if (tau_idx <= tau_of_var)
        list_shifts <- append (list_shifts, n_nodes_shifts)
      else if (!end_reached[[node_idx]]) {
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

#-----------------------------------------------------------------------------
# tmiic.combine_lag
#-----------------------------------------------------------------------------
# @description
# Utility function to combine lags  when flattening the network.
#
# @param df [a dataframe] A non empty dataframe with the edges to combine
#-----------------------------------------------------------------------------
tmiic.combine_lag <- function (df) {
  #
  # Reverse inverted edges (orient == -2) and duplicate lags of bidrectional
  # temporal edges (lag != 0) 
  # NB: for non lag 0 edges, such cases are more than likely errors 
  # and will generate negative lags however showing them will allow the user
  # to identify possible issues
  #
  for (idx in 1:nrow(df) ) { 
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

#-----------------------------------------------------------------------------
# tmiic.combine_orient
#-----------------------------------------------------------------------------
# @description
# Utility function to combine edges orientations when flattening the network.
#
# @param df [a dataframe] The dataframe with the edges to combine
# @param col_name [a string] The orientation column
#-----------------------------------------------------------------------------
tmiic.combine_orient <- function (df, col_name) {
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
# tmiic.combine_probas
#-----------------------------------------------------------------------------
# @description
# Utility function to combine edges probabilities when flattening the network.
# Depending on the combined edges orientation, chose the appropriate max, min
# or mean probabilities to compute the combined edge probabilities
#
# @param df [a dataframe] The dataframe with the edges to combine
# @param comb_orient [a integer] The orientation of the combined edge
#-----------------------------------------------------------------------------
tmiic.combine_probas <- function (df, comb_orient) {
  valid_probas <- grepl (';', df$proba, fixed=TRUE)
  df <- df[valid_probas,]
  if (nrow (df) <= 0)
    return (NA)
  #
  # We set probas like if we have node X <= node Y
  #
  for ( idx in 1:nrow(df) )
    if (df[idx,"x"] > df[idx,"y"]) {
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
#' Note that only the summary data frame is flattened.
#' 
#' @param tmiic.res [a tmiic object] The object returned by miic's 
#' execution in temporal mode.
#' 
#' @param flatten_mode [a string]. Optional, default value "compact".
#' Possible values are \emph{"compact"}, \emph{"combine"}, \emph{"unique"}, 
#' \emph{"drop"}:
#' \itemize{
#' \item When \emph{flatten_mode} = \emph{"compact"}, the default. Nodes 
#'   and edges are converted into a flattened version preserving all 
#'   the initial information.\cr
#'   i.e.: X_lag1->Y_lag0, X_lag0<-Y_lag2 become respectively X->Y lag=1, 
#'   X<-Y lag=2. 
#' \item When \emph{flatten_mode} = \emph{"combine"}, one edge will be kept
#'   per couple of nodes. The info_shifted will be the highest of the 
#'   summarized edges whilst the lag and orientation of the summarized 
#'   edge will be an agregation.\cr 
#'   i.e.: X_lag2->Y_lag0, X_lag0<-Y_lag1 will become X<->Y lag=1,2 with
#'   the info_shifted of X_lag2->Y_lag0 if info_shifted of 
#'   X_lag2->Y_lag0 > X_lag0<-Y_lag1.
#' \item When \emph{flatten_mode} = \emph{"unique"}, only the edges having the
#'   highest info_shifted for a couple of nodes are kept in the flattened 
#'   network. If several edges between the sames nodes have the same
#'   info_shifted, then the edge kept is the one with the minimum lag.\cr 
#'   i.e.: X_lag1->Y_lag0, X_lag0<-Y_lag2 with info_shifted of 
#'   X_lag1->Y_lag0 > X_lag0<-Y_lag2 become X->Y lag=1.
#' \item When \emph{flatten_mode} = \emph{"drop"}, only the edges having the 
#'   highest info_shifted for a couple of nodes are kept in the flattened 
#'   network. If several edges between the sames nodes have the same
#'   info_shifted, then the edge kept is the one with the minimum lag.\cr
#'   i.e. :  X_lag1->Y_lag0, X_lag0<-Y_lag2 with info_shifted of 
#'   X_lag1->Y_lag0 > X_lag0<-Y_lag2 become X->Y. The lag information is 
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
#' @return [a tmiic object] The returned tmiic object is the one received 
#' as input where the summary dataframe has been flattened.
#'     
#' @export
#-----------------------------------------------------------------------------
tmiic.flatten_network <- function (tmiic.res, flatten_mode="compact", 
                                   keep_edges_on_same_node=TRUE) {
  if (! tmiic.res$tmiic_specific[["graph_type"]] == "raw")
    stop ("Error: The miic network must be raw")
  #
  # Keep only edges found by miic
  #
  df_edges <- tmiic.res$all.edges.summary[tmiic.res$all.edges.summary$type %in% c('P', 'TP', 'FP'), ]
  if (nrow(df_edges) <= 0) {
    if (flatten_mode != "drop")
      df_edges$lag = numeric(0)
    tmiic.res$all.edges.summary <- df_edges
    tmiic.res$tmiic_specific[["graph_type"]] = flatten_mode
    return (tmiic.res)
  }
  #
  # Precompute lag and layer of each node
  #
  df_precomputed <- tmiic.precompute_lags_layers_and_shifts (tmiic.res)
  #
  # First step, perform flatten_mode="compact": 
  # from summary, remove lag info from nodes names and put it into a lag column
  #
  df_edges$lag <- -1
  for ( edge_idx in 1:nrow(df_edges) ) {
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
    else {
      df_edges [edge_idx, c("x","y","lag")] <- c(node_y, node_x, -lag)
      if (abs (one_edge$infOrt) == 2)
         df_edges [edge_idx,"infOrt"] <- -one_edge$infOrt
      if ( !is.na (one_edge$trueOrt ) ) 
        if (abs (one_edge$trueOrt ) == 2)
          df_edges [edge_idx,"trueOrt"] <- -one_edge$trueOrt 
      if ( !is.na (one_edge$proba ) ) {
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
  if (nrow(df_edges) <= 0) {
    if (flatten_mode == "drop")
      df_edges$lag <- NULL
    tmiic.res$all.edges.summary <- df_edges
    tmiic.res$tmiic_specific[["graph_type"]] = flatten_mode
    return (tmiic.res)
  }
  #
  # "compact" mode is done
  #
  if (flatten_mode != "compact") {
    #
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
    for ( xy_idx in 1:nrow(df_xy) ) {
      ref_x <- df_xy[xy_idx,"x"]
      ref_y <- df_xy[xy_idx,"y"]
      cond_same_edges = ( ( (df_edges[["x"]] == ref_x) & (df_edges[["y"]] == ref_y) )
                        | ( (df_edges[["x"]] == ref_y) & (df_edges[["y"]] == ref_x) ) )
      df_same <- df_edges[cond_same_edges,]
      
      if (nrow (df_same) > 1) {
        if (flatten_mode == "combine") {
          #
          # Combine lag, orient and proba 
          #
          df_same$new_lag <- tmiic.combine_lag (df_same)
          comb_infOrt <- tmiic.combine_orient (df_same, "infOrt")
          df_same$proba <- tmiic.combine_probas (df_same, comb_infOrt)
          df_same$trueOrt <- tmiic.combine_orient (df_same, "trueOrt")
          df_same$infOrt <- comb_infOrt
          #
          # Orientations and probas have been computed for x <= y, 
          # so force x <= y on all rows
          #
          if (ref_x <= ref_y) {
            df_same[,"x"] <- ref_x
            df_same[,"y"] <- ref_y
          }
          else {
            df_same[, "x"] <- ref_y
            df_same[, "y"] <- ref_x
          }
        }
        
        max_info <- max (df_same[["info_shifted"]])
        df_same <- df_same[ (df_same[["info_shifted"]] == max_info),]
      }
      if (nrow(df_same) > 1) {
        min_lag <- min (df_same[["lag"]])
        df_same <- df_same[ (df_same[["lag"]] == min_lag),]
      }
      if ("new_lag" %in% colnames(df_same) ) {
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
  else {
    #
    # For contextual variable, we clean the lag info
    #
    is_contextual <- tmiic.res$tmiic_specific[["is_contextual"]]
    if (!is.null(is_contextual)) {
      list_nodes_not_lagged = tmiic.res$tmiic_specific[["nodes_not_lagged"]]
      for ( edge_idx in 1:nrow(df_edges) ) {
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
  tmiic.res$all.edges.summary <- df_edges
  tmiic.res$tmiic_specific[["graph_type"]] = flatten_mode
  return (tmiic.res)
}

#-----------------------------------------------------------------------------
# tmiic.repeat_edges_over_history
#-----------------------------------------------------------------------------
#' tmiic.repeat_edges_over_history
#'
#' @description
#' Duplicates edges found by miic over the history assuming stationnarity
#' 
#' @details 
#' In temporal mode, the network returned by miic contains only edges
#' with at least one contemporaneous node (lag0). To improve the visual 
#' aspect when plotting, this function duplicates the edges over the history.  
#' i.e: assuming that we used tau=3, the edge X_lag0-X_lag1 will be 
#' copied as X_lag1-X_lag2 and X_lag2-X_lag3.\cr
#' Note that only the summary data frame is modified.
#' 
#' @param tmiic.res [a tmiic object] The object returned by miic's 
#' execution in temporal mode.
#' 
#' @return [a tmiic object] The tmiic object with a modified summary 
#'     
#' @export
#-----------------------------------------------------------------------------
tmiic.repeat_edges_over_history <- function (tmiic.res) {
  if (! tmiic.res$tmiic_specific[["graph_type"]] == "raw")
    stop ("Error: The miic network must be raw")
  #
  # Consider only edges found by miic  type = "P", "TP", "FP"
  #
  df_edges <- tmiic.res$all.edges.summary[tmiic.res$all.edges.summary$type %in% c('P', 'TP', 'FP'), ]
  if (nrow(df_edges) <= 0) {
    tmiic.res$all.edges.summary <- df_edges
    tmiic.res$tmiic_specific[["graph_type"]] <- "lagged"
    return (tmiic.res)
  }
  #
  # Precompute lag, layer and shift of each node
  #
  df_precomp <- tmiic.precompute_lags_layers_and_shifts (tmiic.res)
  
  list_taus <- tmiic.res$tmiic_specific[["tau"]]
  list_nodes_not_lagged <- tmiic.res$tmiic_specific[["nodes_not_lagged"]]
  #
  # Duplicate the edges over all layers of history
  #
  n_edges <- nrow(df_edges)
  for (edge_idx in 1:n_edges) {
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
    tau_x <- list_taus [[which (list_nodes_not_lagged == node_x_base)]]
    tau_y <- list_taus [[which (list_nodes_not_lagged == node_y_base)]]
    same_lag_needed = TRUE
    if (tau_x <= 0)
      same_lag_needed = FALSE
    if (tau_y <= 0)
      same_lag_needed = FALSE
    #
    # Duplication of the edge
    #
    while (TRUE) {
      #
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
      if (same_lag_needed) {
        new_lag = df_precomp [node_x_pos, "lags"] - df_precomp [node_y_pos, "lags"]
        while (sav_lag != new_lag) {
          if (sav_lag < new_lag) {
            node_y_shift = df_precomp [node_y_pos, "shifts"]
            if (node_y_shift <= 0) {
              same_lag_impossible = TRUE
              break
            }
            node_y_pos = node_y_pos + node_y_shift;
          }
          else { # sav_lag > new_lag
            node_x_shift = df_precomp [node_x_pos, "shifts"]
            if (node_x_shift <= 0) {
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
  tmiic.res$all.edges.summary <- df_edges
  tmiic.res$tmiic_specific[["graph_type"]] <- "lagged"
  return (tmiic.res)
}

