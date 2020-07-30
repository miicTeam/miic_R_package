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
#' tmiic.transform_data_for_miic
#'
#' @description
#' Reorganizes the data using the history to create lagged nodes and  
#' samples in a format usable by miic
#'
#' @details 
#' The function slices the input data according to the \emph{tau}
#' and \emph{delta_tau} parameters. Data are expected to be received in a 
#' dataframe with variables as columns and timeseries/timesteps as rows. 
#' The timestep information must be suplied in the first column and, 
#' for each timeseries, be in ascending order.
#' 
#' The number of variables is increased and renamed on \emph{tau} 
#' / \emph{delta_tau} layers.\cr 
#' i.e. with \emph{tau}=6 and \emph{delta_tau}=3 : node1, node2 => 
#' node1_lag0, node2_lag0, node1_lag3, node2_lag3, node1_lag6, node2_lag6.
#' 
#' Every timestep (until number of timesteps - \emph{tau} * \emph{delta_tau}) 
#' is converted into a sample in the lagged data. Exemple with tau=6 and
#' delta_tau=3:
#' 
#' \tabular{ccccccc}{
#' Timestep \tab  Node & value  \tab  Node & value  \tab => \tab Sample \tab  Node & value  \tab  Node & value \cr
#'   t-6    \tab node1_val(t-6) \tab node2_val(t-6) \tab => \tab   i    \tab node1_lag6_val \tab node2_lag6_val\cr
#'   t-3    \tab node1_val(t-3) \tab node2_val(t-3) \tab => \tab   i    \tab node1_lag3_val \tab node2_lag3_val\cr
#'    t     \tab  node1_val(t)  \tab  node2_val(t)  \tab => \tab   i    \tab node1_lag0_val \tab node2_lag0_val\cr
#'   \cr    \tab                \tab                \tab    \tab        \tab                \tab               \cr
#'   t-7    \tab node1_val(t-7) \tab node2_val(t-7) \tab => \tab   i'   \tab node1_lag6_val \tab node2_lag6_val\cr
#'   t-4    \tab node1_val(t-4) \tab node2_val(t-4) \tab => \tab   i'   \tab node1_lag3_val \tab node2_lag3_val\cr
#'   t-1    \tab node1_val(t-1) \tab node2_val(t-1) \tab => \tab   i'   \tab node1_lag0_val \tab node2_lag0_val\cr
#'   \cr    \tab                \tab                \tab    \tab        \tab                \tab               \cr
#'   t-8    \tab node1_val(t-8) \tab node2_val(t-8) \tab => \tab   i"   \tab node1_lag6_val \tab node2_lag6_val\cr
#'   t-5    \tab node1_val(t-5) \tab node2_val(t-5) \tab => \tab   i"   \tab node1_lag3_val \tab node2_lag3_val\cr
#'   t-2    \tab node1_val(t-2) \tab node2_val(t-2) \tab => \tab   i"   \tab node1_lag0_val \tab node2_lag0_val\cr
#'   \cr    \tab                \tab                \tab    \tab        \tab                \tab               \cr
#'   ...    \tab .............. \tab .............. \tab => \tab ...... \tab .............. \tab ............. \cr
#' }
#' until number of timesteps - \emph{tau} * \emph{delta_tau} is reached. 
#' The same process is applied to all input timeseries.
#'
#' @param data [a dataframe] 
#' A dataframe with the time series with variables as columns and
#' timeseries/timesteps as rows. The timestep information must be suplied 
#' in the first column and, for each timeseries, in an ascending order.
#'
#' @param tau [an int > 0] A strictly positive int defining the max lag.
#' Note that if \emph{delta_tau} is also supplied, \emph{tau} must be a
#' multiple of \emph{delta_tau}.
#' 
#' @param categoryOrder [a data frame] Optional, NULL by default. 
#' A data frame giving information about how to order the various states of 
#' categorical variables. In temporal mode, this data frame will be lagged 
#' as the input data on \emph{tau} / \emph{delta_tau} layers.
#' 
#' @param movavg [an integer] Optional, -1 by default.\cr
#' When \emph{movavg} is supplied (integer > 1), a moving average 
#' operation is applied to the time series.\cr
#' 
#' @param delta_tau [an integer] Optional, 1 by default.\cr
#' When \emph{delta_tau} is supplied (integer > 1), the samples will be 
#' construted using 1 timestep every \emph{delta_tau} timesteps starting 
#' from the last.\cr 
#' i.e.: on 1000 timesteps with  \emph{tau} = 14 and \emph{delta_tau} = 7, 
#' the timesteps kept for the samples conversion will be 1000, 993, 986 
#' for the first sample, the next sample will use 999, 992, 985 and so on.\cr
#' 
#' @return a list with two elements:
#' \itemize{
#'  \item \emph{inputData:} the samples generated from the timeseries
#'  \item \emph{categoryOrer:} the lagged catagoryOrder 
#'  }
#'                                    
#' @export
#' @useDynLib miic
#-----------------------------------------------------------------------------
tmiic.transform_data_for_miic <- function (data, tau, categoryOrder=NULL, 
                                           movavg=-1, delta_tau=1)
  {
  DEBUG <- FALSE
  
  n_rows <- nrow(data)
  n_nodes <- ncol(data) - 1
  list_nodes <- colnames(data)[-1]

  if (DEBUG)
    {
    print ("tmiic.transform_data_for_miic:")
    print (paste ("Nb rows   :", n_rows, sep="") )
    print (paste ("Nb nodes  :", n_nodes, sep="") )
    print (paste ("Tau       :", tau, sep="") )
    print (paste ("movavg    :", movavg, sep="") )
    print (paste ("delta_tau :", delta_tau, sep="") )
    print ("List nodes:")
    print (list_nodes)
    print ("")
    print ("Input df on tau timesteps:")
    print (data[1:tau,])
    print ("")
    print ("Input df: last rows:")
    print (data[(n_rows - tau):n_rows,])
    }
  #
  # Lag the nodes
  #
  list_nodes_lagged <- list()
  for (tau_idx in seq(0,tau,by=delta_tau) )
    for (node_idx in 1:n_nodes)
      list_nodes_lagged <- append (list_nodes_lagged, paste (list_nodes[[node_idx]], "_lag", tau_idx, sep="" ) )
  n_nodes_lagged <- length(list_nodes_lagged)
  if (DEBUG)
    {
    print ("Lagged nodes list")
    print (list_nodes_lagged)
    print (paste ("n_nodes_lagged=", n_nodes_lagged, sep="") )
    }
  #
  # Lag the categoryOrder if supplied
  #
  if (!is.null(categoryOrder))
    {
    categories_lagged <- categoryOrder [FALSE,]
    categ_col1 <- colnames(categoryOrder)[[1]]
    for (tau_idx in seq(0,tau,by=delta_tau) )
      {
      for (old_categ_idx in 1:nrow(categoryOrder))
        {
        lagged_categ_idx <- nrow(categories_lagged) + 1
        categories_lagged [lagged_categ_idx,] <- categoryOrder [old_categ_idx,]
        
        node_name <- categoryOrder [old_categ_idx, categ_col1]
        node_name_lagged <- paste (node_name, "_lag", tau_idx, sep="" )
        categories_lagged [lagged_categ_idx, categ_col1] <- node_name_lagged
        }
      }
    categoryOrder <- categories_lagged
    if (DEBUG)
      {
      print ("Lagged categoryOrder")
      print (categoryOrder)
      }
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
  for (row_idx in 1:n_rows)
    {
    timestep = data[row_idx,1]
    if (timestep < previous_timestep)
      {
      if (DEBUG)
        print (paste ("Timeseries found between ", previous_row_idx, " and ", row_idx-1),
               sep="")
      one_timeseries <- data[previous_row_idx:(row_idx-1), list_nodes]
      df_ret <- miic:::tmiic.lag_one_timeseries (one_timeseries, list_nodes_lagged,
                                                 tau, movavg, delta_tau)
      df_lagged <- rbind (df_lagged, df_ret)
      previous_row_idx <- row_idx
      }
    previous_timestep <- timestep
    }
  if (DEBUG)
    print (paste ("Timeseries found between ", previous_row_idx, " and ", row_idx),
           sep="")
  one_timeseries <- data[previous_row_idx:row_idx, list_nodes]
  df_ret <- miic:::tmiic.lag_one_timeseries (one_timeseries, list_nodes_lagged,
                                             tau, movavg, delta_tau)
  df_lagged <- rbind (df_lagged, df_ret)

  if (DEBUG)
    {
    print ("tmiic.transform_data_for_miic return:")
    print ("Returned df on tau timesteps:")
    print (df_lagged[1:tau,])
    print ("")
    print ("Returned df:  last rows:")
    n_rows_lagged = nrow(df_lagged)
    print (df_lagged[(n_rows_lagged - tau):n_rows_lagged,])
    }
  #
  # returns the dataframe as miic expects and categoryOrder (if supplied) lagged
  #
  return (list (inputData=df_lagged, categoryOrder=categoryOrder) )
  }
  
#-----------------------------------------------------------------------------
# tmiic.lag_one_timeseries
#-----------------------------------------------------------------------------
#' tmiic.lag_one_timeseries
#'
#' @description
#' Reorganizes the data of one timeseries using the history to create lagged
#' samples in a format usable by miic
#'
#' @details 
#' The function slices the input data according to  \emph{tau} and 
#' \emph{delta_tau} parameters. Data are expected to be received in a dataframe 
#' with variables as columns and timesteps as rows. 
#' 
#' Every timestep (until number of timesteps - \emph{tau} * \emph{delta_tau}) 
#' is converted into a sample in the lagged graph. Exemple with tau=6 and
#' delta_tau=3:
#' 
#' \tabular{ccccccc}{
#' Timestep \tab  Node & value  \tab  Node & value  \tab => \tab Sample \tab  Node & value  \tab  Node & value \cr
#'   t-6    \tab node1_val(t-6) \tab node2_val(t-6) \tab => \tab   i    \tab node1_lag6_val \tab node2_lag6_val\cr
#'   t-3    \tab node1_val(t-3) \tab node2_val(t-3) \tab => \tab   i    \tab node1_lag3_val \tab node2_lag3_val\cr
#'    t     \tab  node1_val(t)  \tab  node2_val(t)  \tab => \tab   i    \tab node1_lag0_val \tab node2_lag0_val\cr
#'   \cr    \tab                \tab                \tab    \tab        \tab                \tab               \cr
#'   t-7    \tab node1_val(t-7) \tab node2_val(t-7) \tab => \tab   i'   \tab node1_lag6_val \tab node2_lag6_val\cr
#'   t-4    \tab node1_val(t-4) \tab node2_val(t-4) \tab => \tab   i'   \tab node1_lag3_val \tab node2_lag3_val\cr
#'   t-1    \tab node1_val(t-1) \tab node2_val(t-1) \tab => \tab   i'   \tab node1_lag0_val \tab node2_lag0_val\cr
#'   \cr    \tab                \tab                \tab    \tab        \tab                \tab               \cr
#'   t-8    \tab node1_val(t-8) \tab node2_val(t-8) \tab => \tab   i"   \tab node1_lag6_val \tab node2_lag6_val\cr
#'   t-5    \tab node1_val(t-5) \tab node2_val(t-5) \tab => \tab   i"   \tab node1_lag3_val \tab node2_lag3_val\cr
#'   t-2    \tab node1_val(t-2) \tab node2_val(t-2) \tab => \tab   i"   \tab node1_lag0_val \tab node2_lag0_val\cr
#'   \cr    \tab                \tab                \tab    \tab        \tab                \tab               \cr
#'   ...    \tab .............. \tab .............. \tab => \tab ...... \tab .............. \tab ............. \cr
#' }
#' until number of timesteps - \emph{tau} * \emph{delta_tau} is reached. 
#'
#' @param data [a dataframe] 
#' A dataframe containing the time series with variables as columns and
#' timesteps as rows. 
#'
#' @param tau [an int > 0] A strictly positive int defining the max lag.
#' Note that if \emph{delta_tau} is also supplied, \emph{tau} must be a
#' multiple of \emph{delta_tau}.
#' 
#' @param list_nodes_lagged [a list] 
#' The list of variables lagged over \emph{tau} / \emph{delta_tau}
#' 
#' @param movavg [an integer] Optional, -1 by default.\cr
#' When \emph{movavg} is supplied (integer > 1), a moving average 
#' operation is applied to the time series.
#' 
#' @param delta_tau [an integer] Optional, 1 by default.\cr
#' When \emph{delta_tau} is supplied (integer > 1), the samples will be 
#' construted using 1 timestep every \emph{delta_tau} timesteps starting 
#' from the last.\cr 
#' i.e.: on 1000 timesteps with  \emph{tau} = 14 and \emph{delta_tau} = 7, 
#' the timesteps kept for the samples conversion will be 1000, 993, 986 
#' for the first sample, the next sample will use 999, 992, 985 and so on.\cr
#' 
#' 
#' @return a dataframe with the samples generated form the timeseries
#'                                    
#' @useDynLib miic
#-----------------------------------------------------------------------------
tmiic.lag_one_timeseries <- function (df_timeseries, list_nodes_lagged, 
                                      tau, movavg, delta_tau)
  {
  n_timesteps <- nrow(df_timeseries)
  n_nodes <- ncol(df_timeseries) 
  list_nodes <- colnames(df_timeseries)
  n_nodes_lagged <- length(list_nodes_lagged)
    
  DEBUG <- FALSE
  if (DEBUG)
    {
    print ("tmiic.lag_one_timeseries:")
    print (paste ("n_timesteps :", n_timesteps, sep="") )
    print (paste ("Nb nodes   :", n_nodes, sep="") )
    print (paste ("Nb nodes lagged  :", n_nodes_lagged, sep="") )
    print (paste ("Tau        :", tau, sep="") )
    print (paste ("movavg     :", movavg, sep="") )
    print (paste ("delta_tau  :", delta_tau, sep="") )
    print ("")
    print ("Input df on tau timesteps:")
    print (df_timeseries[1:tau,])
    print ("")
    print ("Input df: last rows:")
    print (df_timeseries[(n_timesteps - tau):n_timesteps,])
    }
  #
  # Init df to be returned
  #
  df_lagged <- data.frame (matrix(ncol=n_nodes_lagged, nrow=0), 
                           stringsAsFactors=FALSE)
  colnames (df_lagged) <- list_nodes_lagged
  #
  # Apply moving average if requested
  #
  if (movavg > 1)
    {
    if (DEBUG)
      {
      print ("movavg:")
      print ("Original df:")
      print (df_timeseries[1:movavg,])
      print (df_timeseries[(n_timesteps-movavg):n_timesteps,])
      }
    df_movavg <- data.frame (matrix(ncol=n_nodes, nrow=n_timesteps), 
                             stringsAsFactors=FALSE)
    colnames (df_movavg) <- list_nodes
    for (node_idx in 1:n_nodes)
      df_movavg[, node_idx] <- data.table::frollmean (df_timeseries[,node_idx], movavg)
    
    if (DEBUG)
      {
      print ("After raw movavg:")
      print (dim(df_movavg))
      print (df_movavg[1:movavg,])
      print (df_movavg[(n_timesteps-movavg):n_timesteps,])
      }
    #
    # Do not keep the first timesteps as NA due to the moving average
    #
    df_timeseries <- df_movavg[movavg:n_timesteps,]
    n_timesteps <- nrow(df_timeseries)
    if (DEBUG)
      {
      print ("After movavg:")
      print (dim(df_timeseries))
      print (df_timeseries[1:movavg,])
      print (df_timeseries[(n_timesteps-movavg):n_timesteps,])
      }
    }
  #
  # Loop over timesteps to create lagged samples
  #
  if (DEBUG)
    print (paste ("idx min:", tau+1, sep="") )
  for (row_idx in seq(from=n_timesteps, to=tau+1, by=-1) )
    {
    new_line <- nrow(df_lagged)+1
    for (node_idx in 1:n_nodes)
      for (tau_idx in seq(0,tau,by=delta_tau) )
        {
        lagged_node_idx <- (tau_idx %/% delta_tau) * n_nodes + node_idx
        df_lagged [new_line,lagged_node_idx] <- df_timeseries[row_idx-tau_idx,node_idx]
        }
    }
 
  if (DEBUG)
    {
    print ("")
    print (paste ("Check: Input data over first tau*3", sep="") )
    print ( df_timeseries [1:(tau*3),] )
    print ("Check, output over last tau")
    print (df_lagged [(nrow(df_lagged)-tau):nrow(df_lagged),])
    print ("")
    print ("Check: Input data over last tau*3")
    print ( df_timeseries [(n_timesteps-tau*3):n_timesteps,] )
    print ("Check: Output data over tau")
    print (df_lagged [1:tau,])
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
#' Note that only the adjMatrix and summary data frames are flatened.
#' 
#' @param miic_return [a miic object] The object returned by miic's 
#' execution.
#' 
#' @param flatten_mode [a string]. Optional, default value "normal".
#' Possible values are \emph{"normal"}, \emph{"combine"}, \emph{"unique"}, 
#' \emph{"drop"}:
#' \itemize{
#' \item When \emph{flatten_mode} = \emph{"normal"}, the default, edges are 
#'   simply converted into a flattened version.\cr
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
#' @param keep_edges_on_same_node [a boolean] Optional, TRUE by default.
#' When TRUE, the edges like X_lag0-X_lag1 are kept during flatenning
#' (it becomes an X-X edge). When FALSE, only edges having different nodes 
#' are kept in the flatten network.
#' 
#' @return [a list] The returned list is the one received as input where
#' the adjacency and summary dataframes has been flattened.\cr
#' Note that the adjacency matrix does not contain any more the orientation
#' as this information make less sense after flattening. Orientations are 
#' replaced by lag(s) in the adjacency matrix.
#'     
#' @importFrom magrittr "%>%"                             
#' @export
#' @useDynLib miic
#-----------------------------------------------------------------------------
tmiic.flatten_network <- function (miic_return, flatten_mode="normal", 
                                   keep_edges_on_same_node=TRUE)
  {
  DEBUG <- FALSE
  if (DEBUG)
    {
    print ("tmiic.flatten_network:")
    print ("adjacency matrix:")
    print (miic_return$adjMatrix)
    print ("summary:")
    print (miic_return$all.edges.summary %>% dplyr::select(x,y,type,infOrt,sign,log_confidence,proba) )
    print (paste ("flatten_mode=", flatten_mode) )
    }
  #
  # Construct the list of nodes not lagged 
  #
  list_nodes_lagged = colnames (miic_return$adjMatrix)
  list_nodes_not_lagged = list ()
  for (node_idx in 1:length (list_nodes_lagged) )
    {
    #
    # As long we find "_lag0" at the end of the node name
    #
    one_node = list_nodes_lagged[[node_idx]]
    if ( endsWith(one_node, "_lag0") )
      {
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
  
  if (DEBUG)
    {
    print ("list of node not lagged:")
    print (list_nodes_not_lagged)
    print (paste ("number of nodes not lagged=", n_nodes) )
    }
  #
  # First step, perform flatten_mode="normal" : 
  # - on adjacency matrix : construct a n_nodes x n_nodes with lags in each cell
  # - on summary: remove lag info from nodes names and put it into a lag column
  #
  df_adj <- array ( data="", dim=c(n_nodes ,n_nodes),
                    dimnames=list(list_nodes_not_lagged, list_nodes_not_lagged) )
  df_edges <- miic_return$all.edges.summary
  df_edges$lag <- -1
  for (edge_idx in 1:nrow(df_edges) )
    {
    #
    # Get nodes the edge wthout "_lagX"
    #
    one_edge <- df_edges[edge_idx,]
    #
    # Edges with type != "P" are not true edges
    #
    if (one_edge$type != "P")
      {
      if (DEBUG)
        {
        print ( paste("Edge ", one_edge$x, "-", one_edge$y, ", type=", one_edge$type,
                      " => type!='P' => dropped", sep="") )
        }
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
    if ( !is.na(pos_lag_x[1]) )
      {
      tau_idx_x <- stringr::str_replace(node_x, ".*_lag", "")
      tau_idx_x <- strtoi (tau_idx_x)
      }
    pos_lag_y <- stringr::str_locate(node_y, "_lag")
    tau_idx_y <- 0
    if ( !is.na(pos_lag_y[1]) )
      {
      tau_idx_y <- stringr::str_replace(node_y, ".*_lag", "")
      tau_idx_y <- strtoi (tau_idx_y)
      }
    #
    # It temporal, ensure to order from oldest to newest
    #
    lag = tau_idx_x - tau_idx_y
    if (lag < 0)
      {
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
    if ( (!keep_edges_on_same_node) & (node_x == node_y) )
      {
      if (DEBUG)
        {
        print ( paste("Edge ", one_edge$x, "-", one_edge$y, ", type=", one_edge$type,
                      " => node_x==node_y' => dropped", sep="") )
        }
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
  for (i in 1:ncol (df_adj) )
    {
    for (j in i:ncol (df_adj) )
      {
      lags_str <- df_adj[i,j]
      if (nchar (lags_str) > 0)
        if ( grepl (',', lags_str, fixed = TRUE) )
          {
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
  if (DEBUG)
    {
    print ("End step 1: lag moved into its own column")
    print ("df_adj:")
    print (df_adj)
    print ("df_edges:")
    print (df_edges %>% dplyr::select(x,y,lag,type,infOrt,sign,log_confidence,proba))
    }
  #
  # If we want only one edge per couple of nodes 
  #
  if (flatten_mode != "normal")
    {
    # At first, we ensure X string < Y string and we create an extra
    # column containing x-y so we have an id to perform a duplicate check
    #
    for (edge_idx in 1:nrow(df_edges)) 
      {
      one_edge <- df_edges[edge_idx,]
      if (one_edge$x < one_edge$y)
        {
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
    if (DEBUG) 
      {
      print ("df_edges after ordering x > y:")
      print (df_edges %>% dplyr::select(x,y,lag,type,infOrt,trueOrt,sign,log_confidence,proba) )
      }
    #
    # Keep the rows having max log_confidence when grouped by on x, y
    # If some rows between same nodes have the same log_confidence 
    # We addd a second step keeping the edges with the minimum lag 
    #
    df_group <- df_edges %>% dplyr::group_by(x,y) %>% dplyr::top_n(n=1, wt=log_confidence)
    if (DEBUG) 
      {
      print ("df_group non unique:")
      print (df_group %>% dplyr::select(x,y,lag,type,infOrt,trueOrt,sign,log_confidence,proba) )
      }
    df_group <- df_group %>% dplyr::group_by(x,y) %>% dplyr::top_n(n=1, wt=-lag)
    df_group <- as.data.frame (df_group)
    if (DEBUG) 
      {
      print ("df_group:")
      print (df_group %>% dplyr::select(x,y,lag,type,infOrt,trueOrt,sign,log_confidence,proba) )
      }
    #
    # If the flatten_mode is "unique" or "drop", nothing more to do for now
    #
    # If the flatten_mode is "combine", combine orientations
    # 
    if (flatten_mode == "combine")
      {
      #
      # Update the orientions (in case the grouped edges have different orientations)
      # and put in the lag column the list of lags
      #
      for (edge_idx in 1:nrow(df_group) )
        {
        one_edge <- df_group[edge_idx,]
        node_x <- one_edge$x
        node_y <- one_edge$y
        lag <- one_edge$lag
        if (DEBUG)
          {
          print ("combine: processing grouped edge:")
          print (one_edge %>% dplyr::select(x,y,lag,type,infOrt,trueOrt,sign,log_confidence,proba) )
          }
        #
        # Update the list of lags (value has been computed in lag matrix)
        #
        df_group[edge_idx,]$lag <- df_adj[node_x, node_y]
        #
        # Select the other edges between same nodes for orientation update
        #
        cond_for_orient <- ( (df_edges[["x"]] == node_x) 
                           & (df_edges[["y"]] == node_y) 
                           & (df_edges[["lag"]] != lag) )
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
        if ( (abs(one_edge$infOrt) >= 2) & (!is.na(one_edge$proba)) )
          {
          edge_probas <- strsplit (one_edge$proba, ';')[[1]]
          if (one_edge$infOrt > 0)
            {
            head_proba <- edge_probas[[1]]
            tail_proba <- edge_probas[[2]]
            }
          else
            {
            head_proba <- edge_probas[[2]]
            tail_proba <- edge_probas[[1]]
            }
          }
        if (DEBUG)
          print (paste ("grouped edge, head=", head_proba, " tail=", tail_proba, sep="") )
        
        df_other <- df_edges[cond_for_orient,]
        for (other_idx in 1:nrow(df_other) )
          {
          other_edge <- df_other[other_idx,]
          if (DEBUG)
            {
            print ("combine: processing other_edge edge:")
            print (other_edge %>% dplyr::select(x,y,lag,type,infOrt,trueOrt,sign,log_confidence,proba) )
            }
          #
          # Update infOrt and trueOrt
          #
          for (col_to_update in c("infOrt", "trueOrt"))
            {
            other_edge_orient <- other_edge[[col_to_update]]
            group_edge_orient <-  df_group[edge_idx,][[col_to_update]]
            #
            # If other edge is na has no orient, nothing to do
            #
            if ( (is.na(other_edge_orient)) | (other_edge_orient == 1) )
              {
              if (DEBUG)
                print (paste ("Update ", col_to_update, " orientation: ", 
                       node_x, "-", node_y, " lag=", lag, 
                       " orient=", group_edge_orient, 
                       ", other edge orientation=", other_edge_orient, 
                       " => na or unoriented, nothing to do", sep="") )
              next            
              }
            #
            # If grouped edge is na or has no orient, we can always update
            #
            if (  (is.na(group_edge_orient)) | (group_edge_orient == 1) )
              {
              if (DEBUG)
                print (paste ("Update ", col_to_update, " orientation: ", 
                       node_x, "-", node_y, " lag=", lag,
                       " orient=", group_edge_orient, 
                       ", other edge orientation=", other_edge_orient, 
                       " => update to ", other_edge_orient, sep="") )
              df_group[edge_idx,][[col_to_update]] <- other_edge_orient
              next            
              }
            #
            # The 2 edges have an orientation (unidirectional or birectional)
            #
            # If the other edge is bidirectional or opposite to the grouped one 
            # -> update to bidirectional
            #
            if ( (other_edge_orient == 6) | (other_edge_orient == -group_edge_orient) )
              {
              if (DEBUG)
                print (paste ("Update ", col_to_update, " orientation: ", 
                       node_x, "-", node_y, " lag=", lag,
                       " orient=", group_edge_orient, 
                       ", other edge orientation=", other_edge_orient, 
                       " => update to birectional=6", sep="") )
              df_group[edge_idx,][[col_to_update]] <- 6
              next            
              }
            if (DEBUG)
              print (paste ("Update ", col_to_update, " orientation: ", 
                     node_x, "-", node_y, " lag=", lag,
                     " orient=", group_edge_orient, 
                     ", other edge orientation=", other_edge_orient, 
                     " => nothing done", sep="") )
            }
          #
          # Keep max proba of tail/queue
          #
          other_edge_head_proba <- 0
          other_edge_tail_proba <- 0
          if ( (abs(other_edge$infOrt) >= 2) & (!is.na(other_edge$proba)) )
            {
            other_edge_probas <- strsplit (other_edge$proba, ';')[[1]]
            if (other_edge$infOrt > 0)
              {
              other_edge_head_proba <- other_edge_probas[[1]]
              other_edge_tail_proba <- other_edge_probas[[2]]
              }
            else
              {
              other_edge_head_proba <- other_edge_probas[[2]]
              other_edge_tail_proba <- other_edge_probas[[1]]
              }
            if (DEBUG)
              print (paste ("other edge, probas=", other_edge_probas, " orient=", other_edge$infOrt,
                            " head=", other_edge_head_proba, " tail=", other_edge_tail_proba, sep="") )
            }
          if (other_edge_head_proba > head_proba)
            head_proba <- other_edge_head_proba
          if (other_edge_tail_proba > tail_proba)
            tail_proba <- other_edge_tail_proba
          if (DEBUG)
            print (paste ("afer update new head=", head_proba, 
                          " new tail=", tail_proba, sep="") )
          }
        #
        # Update the grouped edge with the new probas
        #
        if (df_group[edge_idx,]$infOrt >= 0)
          df_group[edge_idx,]$proba <- paste(head_proba, tail_proba, sep=";")
        else
          df_group[edge_idx,]$proba <- paste(tail_proba, head_proba, sep=";")
        if (DEBUG)
          print (paste ("updated edge, probas=", df_group[edge_idx,]$proba, sep="") )
        }
      }
      
    df_edges <- df_group
    if (DEBUG)
      {
      print (paste ("after flatten_mode == ", flatten_mode) )
      print (df_edges %>% dplyr::select(x,y,lag,type,infOrt,trueOrt,sign,log_confidence,proba) )
      }
    }
  #
  # If we do not want to keep info about lag at all
  #
  if (flatten_mode == "drop")
    {
    if (DEBUG)
      {
      df_edges$lag <- NULL
      print ("after flatten_mode == 'drop':")
      print (df_edges %>% dplyr::select(x,y,type,infOrt,sign,log_confidence,proba))
      }
    }
  #
  # returns the list of dataframes where network has been flattened
  #
  miic_return$adjMatrix <- df_adj
  miic_return$all.edges.summary <- df_edges
  if (DEBUG)
    {
    print ("Returned values:")
    print ("adjMatrix:")
    print (miic_return$adjMatrix)
    print ("Summary:")
    print (miic_return$all.edges.summary %>% dplyr::select(x,y,lag,type,infOrt,sign,log_confidence,proba))
    }
  return (miic_return)
  }

