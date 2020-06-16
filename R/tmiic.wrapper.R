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
# - 15 june 2020 : add subtiming and moving average
#*****************************************************************************
library(dplyr)

#-----------------------------------------------------------------------------
# tmiic.transform_data_for_miic
#-----------------------------------------------------------------------------
#' tmiic.transform_data_for_miic
#'
#' @description
#' Reorganizes the data using the history to create lagged nodes and extra 
#' samples in a format usable by miic
#'
#' @details 
#' The function slices the input data according to the lag max argument 
#' \emph{tau}. Data are expected to be received in a 3 dimensional array 
#' [n_samples * n_nodes * n_time]. History is assumed to be time ordered from 
#' the oldest (first rows) to lastest (ending rows).
#' 
#' The number of nodes is increased and renamed on \emph{tau} layers.\cr 
#' i.e. with \emph{tau}=2: node1, node2 => node1_lag0, node2_lag0, node1_lag1, 
#' node2_lag1, node1_lag2, node2_lag2.
#' 
#-----------------------------------------------------------------------------
#' Every timestep (until number of timesteps - \emph{tau}) is converted into 
#' a sample in the lagged graph:
#' 
#' \tabular{ccccccc}{
#' Timestep \tab  Node & value  \tab  Node & value  \tab => \tab Sample \tab  Node & value  \tab  Node & value \cr
#'   t-2    \tab node1_val(t-2) \tab node2_val(t-2) \tab => \tab   i    \tab node1_lag2_val \tab node2_lag2_val\cr
#'   t-1    \tab node1_val(t-1) \tab node2_val(t-1) \tab => \tab   i    \tab node1_lag1_val \tab node2_lag1_val\cr
#'    t     \tab  node1_val(t)  \tab  node2_val(t)  \tab => \tab   i    \tab node1_lag0_val \tab node2_lag0_val\cr
#'   \cr    \tab                \tab                \tab    \tab        \tab                \tab               \cr
#'   t-3    \tab node1_val(t-3) \tab node2_val(t-3) \tab => \tab   i'   \tab node1_lag2_val \tab node2_lag2_val\cr
#'   t-2    \tab node1_val(t-2) \tab node2_val(t-2) \tab => \tab   i'   \tab node1_lag1_val \tab node2_lag1_val\cr
#'   t-1    \tab node1_val(t-1) \tab node2_val(t-1) \tab => \tab   i'   \tab node1_lag0_val \tab node2_lag0_val\cr
#'   \cr    \tab                \tab                \tab    \tab        \tab                \tab               \cr
#'   t-4    \tab node1_val(t-4) \tab node2_val(t-4) \tab => \tab   i"   \tab node1_lag2_val \tab node2_lag2_val\cr
#'   t-3    \tab node1_val(t-3) \tab node2_val(t-3) \tab => \tab   i"   \tab node1_lag1_val \tab node2_lag1_val\cr
#'   t-2    \tab node1_val(t-2) \tab node2_val(t-2) \tab => \tab   i"   \tab node1_lag0_val \tab node2_lag0_val\cr
#'   \cr    \tab                \tab                \tab    \tab        \tab                \tab               \cr
#'   ...    \tab .............. \tab .............. \tab => \tab ...... \tab .............. \tab ............. \cr
#' }
#' until number of timesteps - \emph{tau} is reached. The same process is applied
#' to all input samples.
#'
#' @param data_tab [a 2D or 3D array] 
#' An array of the time series for all nodes and samples of dimensions
#' [n_samples * n_nodes * n_time] for 3D or [n_nodes * n_time] for 2D.
#' When a 2D array is supplied, the number of samples is assumed to be 1.
#'
#' @param tau [an int > 0] A strictly positive int defining the max lag.
#' 
#' @param categoryOrder [a data frame] Optional, NULL by default. 
#' A data frame giving information about how to order the various states of 
#' categorical variables. This data frame will be lagged as the input 
#' data on \emph{tau} timesteps.
#' 
#' @param movavg [an integer] Optional, -1 by default.\cr
#' When \emph{movavg} is supplied (integer > 1), a moving average 
#' operation is applied to the time series.\cr
#' Note that when both moving average and subtiming are applied,
#' the moving average is performed before the subtiming.
#' 
#' @param subtiming [an integer] Optional, -1 by default.\cr
#' When \emph{subtiming} is supplied (integer > 1), the time series will be 
#' subtimed using 1 timestep every \emph{subtiming} timesteps starting from
#' the last.\cr 
#' i.e.: on 1000 timesteps with \emph{subtiming} = 7, the  timesteps kept
#' will be 1000, 993, 986, ..., 13, 6.\cr
#' Note that when both moving average and subtiming are applied,
#' the moving average is performed before the subtiming.
#' 
#' @param bootstrap [an int] Optional, default=-1.\cr
#' Experimental. When -1, no bootstraping is performed. 
#' When > 0, select randomly \emph{bootstrap} lagged samples (the samples 
#' obtained after transformation of the input samples over \emph{tau} 
#' timesteps).\cr
#' As normal when using bootstrapping, the \emph{bootstrap} value can be
#' greater than the number of lagged samples as a lagged sample can be  
#' selected more than once.
#' 
#' @return a 2D array of dimensions [ n_samples * (n_time -  \emph{tau}) ), 
#'                                    n_nodes * ( \emph{tau}+1) ]
#-----------------------------------------------------------------------------
tmiic.transform_data_for_miic <- function (data_tab, tau, categoryOrder=NULL, 
                                           movavg=-1, subtiming=-1, bootstrap=-1)
  {
  DEBUG <- FALSE
  if (DEBUG)
    print ("tmiic.transform_data_for_miic:")
  #
  # If input is a 2D array, convert it to 3D with 1 sample
  #
  dim_data <- dim(data_tab)
  if (length(dim_data) == 2)
    {
    tmp_dim_names <- dimnames(data_tab)
    dim(data_tab) <- c(1, dim_data)
    dimnames(data_tab) <- c(c(1), tmp_dim_names) 
    if (DEBUG)
      {
      print ("Change 2D into 3D:")
      print (dim_data)
      print ( dim(data_tab) )
      }
    dim_data <- dim(data_tab)
    }
  n_samples <- dim_data[[1]]
  n_nodes <- dim_data[[2]]
  n_time <- dim_data[[3]]
  list_nodes <- dimnames(data_tab)[[2]]

  # file_trace <- file("trace.txt")
  # writeLines (paste ("tmiic.transform_data_for_miic called at ", Sys.time(), "\n",
  #                    "Nb samples=", n_samples, "\n",
  #                    "Nb nodes=", n_nodes, "\n",
  #                    "Nb timesteps=", n_time, "\n",
  #                    "Tau max=", tau, "\n",
  #                    "Mov avg=", movavg, "\n",
  #                    "Subtiming=", subtiming, "\n",
  #                    "Bootstrap=", bootstrap, 
  #                    sep=""), file_trace)
  # close (file_trace)
  
  if (DEBUG)
    {
    print (paste ("Nb samples  :", n_samples, sep="") )
    print (paste ("Nb nodes    :", n_nodes, sep="") )
    print (paste ("Nb timesteps:", n_time, sep="") )
    print (paste ("Tau         :", tau, sep="") )
    print (paste ("Moving avg  :", movavg, sep="") )
    print (paste ("Subtiming   :", subtiming, sep="") )
    print (paste ("Bootstrap   :", bootstrap, sep="") )
    print ("")
    print ("Input df: sample 1 over tau timesteps:")
    print (data_tab[1,,1:tau])
    print ("")
    print ("Input df: sample 1 over last tau timesteps:")
    print (data_tab[1,,(n_time - tau + 1):n_time])
    }
  #
  # If moving average is requested
  #
  if (movavg > 1)
    {
    if (DEBUG)
      {
      print ("movavg:")
      print ("Original tab:")
      print (data_tab[1,,1:movavg])
      print (data_tab[1,,(n_time-movavg):n_time])
      }
    tab_movavg <- array ( data=NA, dim=c (n_samples, n_nodes, n_time),
            dimnames=list(seq(1,n_samples), list_nodes, seq(1,n_time)) )
    for (sample_idx in 1:n_samples)
      for (node_idx in 1:n_nodes)
        {
        tab_movavg[sample_idx,node_idx,] <- frollmean (data_tab[sample_idx,node_idx,], movavg)
        }
    if (DEBUG)
      {
      print ("Tab after raw movavg:")
      print (dim(tab_movavg))
      print (tab_movavg[1,,1:movavg])
      print (tab_movavg[1,,(n_time-movavg):n_time])
      }
    #
    # Do not keep the first timesteps as NA due to the moving average
    #
    data_tab <- tab_movavg[,,movavg:n_time]
    #
    # Switch again to 3D when only one sample 
    #
    dim_data <- dim(data_tab)
    if (length(dim_data) == 2)
      {
      tmp_dim_names <- dimnames(data_tab)
      dim(data_tab) <- c(1, dim_data)
      dimnames(data_tab) <- c(c(1), tmp_dim_names) 
      if (DEBUG)
        {
        print ("Change 2D into 3D:")
        print (dim_data)
        print ( dim(data_tab) )
        }
      dim_data <- dim(data_tab)
      }
    #
    # Rename the timestep from 1 to n_time - movavg + 1
    #
    n_time <- dim_data[[3]]
    dimnames(data_tab)[[3]] <- seq (1,n_time)  
    if (DEBUG)
      {
      print ("Tab after movavg:")
      print (dim(data_tab))
      print (paste ("Nb timesteps:", n_time, sep="") )
      print (data_tab[1,,1:movavg])
      print (data_tab[1,,(n_time-movavg):n_time])
      }
    }
  #
  # If subtiming is requested : keep only keep 1 every x timesteps
  #
  if (subtiming > 1)
    {
    if (DEBUG)
      {
      print ("Sub timing:")
      print ("Original tab, sample 1:")
      print (data_tab[1,,1:subtiming])
      print (data_tab[1,,(n_time-subtiming):n_time])
      }
    #
    # Subtiming keeping the last timestep 
    # (=> skip some timesteps to end on the last one) 
    #
    nb_skip <- n_time - (n_time %/% subtiming) * subtiming
    if (nb_skip == 0)
      nb_skip <- subtiming
    timesteps_kept <- seq.int (from=nb_skip, to=n_time, by=subtiming) 
    data_tab <- data_tab[,,timesteps_kept]
    #
    # Switch again to 3D when only one sample 
    #
    dim_data <- dim(data_tab)
    if (length(dim_data) == 2)
      {
      tmp_dim_names <- dimnames(data_tab)
      dim(data_tab) <- c(1, dim_data)
      dimnames(data_tab) <- c(c(1), tmp_dim_names) 
      if (DEBUG)
        {
        print ("Change 2D into 3D:")
        print (dim_data)
        print ( dim(data_tab) )
        }
      dim_data <- dim(data_tab)
      }
    n_time <- dim_data[[3]]
    if (DEBUG)
      {
      print ("Tab after subtiming, sample 1:")
      print (dim(data_tab))
      print (paste ("Nb timesteps:", n_time, sep="") )
      print (data_tab[1,,1:2])
      print (data_tab[1,,(n_time-1):n_time])
      }
    #
    # Update dimnames so we go from 1 to n_time subtimed
    #
    dimnames(data_tab)[[3]] <- seq (1,n_time)  
    if (DEBUG)
      {
      print ("Tab after subtiming and colnames update, sample 1:")
      print (dim(data_tab))
      print (paste ("Nb timesteps:", n_time, sep="") )
      print (data_tab[1,,1:2])
      print (data_tab[1,,(n_time-1):n_time])
      }
    }
  # if (DEBUG)
  #   tmiic.plot_one_sample (data_tab, 1, title="Post ", filename="./Sample1_Post.png")
  #
  # Add the lagged nodes names to the nodes list 
  #
  new_list_nodes <- list()
  for (tau_idx in 0:tau)
    {
    for (node_idx in 1:n_nodes)
      {
      new_list_nodes <- append (new_list_nodes, paste (list_nodes[[node_idx]], "_lag", tau_idx, sep="" ) )
      if (DEBUG)
        {
        print (paste ("new nodes=", new_list_nodes[[length(new_list_nodes)]], sep="") )
        }
      }
    }
  new_n_nodes <- n_nodes * (tau + 1)
  if (DEBUG)
    {
    print (paste ("new nb nodes=", new_n_nodes, sep="") )
    }
  #
  # Transform the other parameter to fit with the new list of nodes
  #
  if (!is.null(categoryOrder))
    {
    new_categories <- categoryOrder [FALSE,]
    categ_col1 <- colnames(categoryOrder)[[1]]
    for (tau_idx in 0:tau)
      {
      for (old_categ_idx in 1:nrow(categoryOrder))
        {
        new_categ_idx <- tau_idx * nrow(categoryOrder) + old_categ_idx
        new_categories [new_categ_idx,] <- categoryOrder [old_categ_idx,]
        new_categories [new_categ_idx, categ_col1] <- paste (categoryOrder [old_categ_idx, categ_col1], 
                                                      "_lag", tau_idx, sep="" ) 
        if (DEBUG)
          {
          print (paste ("old categ=", categoryOrder [old_categ_idx, categ_col1],
                        " new categ=", new_categories [new_categ_idx, categ_col1], sep="") )
          }
        }
      }
    if (DEBUG)
      {
      print ("New categ df=")
      print (paste ("new nb nodes=", new_n_nodes, sep="") )
      }
    categoryOrder <- new_categories
    }
  #
  # Create the new array with dimensions = 
  # [ n_samples * (ntime %/% tau_plus_1), new_n_nodes
  #
  tau_plus_1 <- tau + 1
  n_samples_from_history = n_time - tau
  
  tab_lagged <- array ( data=NA, dim=c (n_samples * n_samples_from_history, new_n_nodes),
                        dimnames=list(seq(1,n_samples * n_samples_from_history), new_list_nodes) )
  if (DEBUG)
    {
    print ("create empty array of dim: ")
    print ( dim(tab_lagged) )
    }
  #
  # Reallocate data into new array
  #
  for (new_sample_idx in 1:(n_samples_from_history * n_samples) )
    {
    for (new_node_idx in 1:new_n_nodes)
      {
      tau_idx <- (new_node_idx-1) %/% n_nodes
      old_node_idx <- ( (new_node_idx-1) %% n_nodes ) + 1
      
      old_sample_idx <- ( (new_sample_idx-1) %/% n_samples_from_history) + 1
      old_time_idx <- ( (new_sample_idx-1)  %% n_samples_from_history) + tau_idx + 1
      old_time_rev_idx <- n_time + 1 - old_time_idx
      if ( (DEBUG) & (   (new_sample_idx <= 3) | (new_sample_idx == n_samples_from_history) 
                       | (new_sample_idx == 148) | (new_sample_idx == n_samples_from_history * n_samples) ) )
        {
        print("")
        print (paste ("new_sample_idx=", new_sample_idx, " new_node_idx=", new_node_idx, 
                      " node=", new_list_nodes[[new_node_idx]], sep="") )        
        print (paste ("=> old_node_idx=", old_node_idx, " node=", list_nodes[[old_node_idx]], 
                      " tau_idx=", tau_idx, sep="") ) 
        print (paste ("=> old_sample_idx=", old_sample_idx, " old_time_idx=", old_time_idx,
                      " => old_time_rev_idx=", old_time_rev_idx, sep="") ) 
        }
      tab_lagged[new_sample_idx, new_node_idx] <- data_tab [old_sample_idx, old_node_idx, old_time_rev_idx]
      }
    }
  
  if (DEBUG)
    {
    print ("")
    print ("Check sample 1 at t : Input data over last tau")
    print ( data_tab [1, , (n_time-tau):n_time] )
    print ("Check sample 1: Output data")
    print (tab_lagged [1,])
    print ("")
    print ("Check sample 1 a t-1 : Input data over last tau+1")
    print ( data_tab [1, , (n_time-tau-1):n_time] )
    print ("Check sample 1: Output data for new sample 2")
    print (tab_lagged [2,])
    print ("")
    print (paste ("Check sample 1: Input data over first tau", sep="") )
    print ( data_tab [1, , 1:(tau*2)] )
    print ("Check sample1: Output data, last new sample from sample 1")
    print (tab_lagged [n_samples_from_history,])
    print ("")
    print (paste ("Check last sample: Input data over last tau", sep="") )
    print ( data_tab [n_samples, , (n_time-tau):n_time] )
    print ("Check sample1: Output data")
    print (tab_lagged [(n_samples_from_history * (n_samples-1) + 1),])
    print ("")
    print (paste ("Check last sample: Input data over first tau", sep="") )
    print ( data_tab [n_samples, , 1:(tau*2)] )
    print ("Check sample1: Output data of last new sample")
    print (tab_lagged [n_samples_from_history * n_samples,])
    print ("")
    print ("Return df for sample 1 and last:")
    print (tab_lagged[1,])
    print (tab_lagged[(n_samples_from_history * n_samples),])
    }
  #
  # If Bootstrapping
  #
  if (bootstrap > 0)
    {
    tab_boost <- array ( data=NA, dim=c (bootstrap, new_n_nodes),
                          dimnames=list(seq(bootstrap), new_list_nodes) )
    if (DEBUG)
      {
      print ("Bootstrapping:")
      print ("tab samples avt:")
      print (dim(tab_lagged))
      print (paste ("nb samples orig:", n_samples) )
      print ("create empty array of dim: ")
      print ( dim(tab_boost) )
      }
    #
    # To fill the bootstrapped tab, take randomly samples until arry filled
    #
    for (boost_idx in 1:bootstrap)
      {
      rand_idx = sample(1:(n_samples * n_samples_from_history), 1)
      tab_boost[boost_idx,] <- tab_lagged[rand_idx,]
      }
    tab_lagged <- tab_boost
    if (DEBUG)
      {
      print ("Tab samples ap:")
      print (dim(tab_lagged))
      print (tab_lagged[1,])
      print (tab_lagged[boost_idx,])
      }
    }
  #
  # returns the dataframe as miic expects and categoryOrder (if supplied) lagged
  #
  return (list (inputData=as.data.frame(tab_lagged, stringsAsFactors = FALSE),
                categoryOrder=categoryOrder) )
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
#' (X_lag0, X_lag1, ...) over tau timesteps. This function flatten the 
#' network depending of the \emph{flatten_mode} parameter.\cr
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
#'   network.\cr 
#'   i.e.: X_lag1->Y_lag0, X_lag2<-Y_lag0 with log_confidence of 
#'   X_lag1->Y_lag0 > X_lag2<-Y_lag0 become X->Y lag=1.
#' \item When \emph{flatten_mode} = \emph{"drop"}, only the edges having the 
#'   highest log_confidence for a couple of nodes are kept in the flattened 
#'   network.\cr
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
    print (miic_return$all.edges.summary %>% select(x,y,type,infOrt,sign,log_confidence) )
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
    pos_lag_x <- str_locate(node_x, "_lag")
    tau_idx_x <- 0
    if ( !is.na(pos_lag_x[1]) )
      {
      tau_idx_x <- str_remove(node_x, ".*_lag")
      tau_idx_x <- strtoi (tau_idx_x)
      }
    pos_lag_y <- str_locate(node_y, "_lag")
    tau_idx_y <- 0
    if ( !is.na(pos_lag_y[1]) )
      {
      tau_idx_y <- str_remove(node_y, ".*_lag")
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
    print (df_edges %>% select(x,y,lag,type,infOrt,sign,log_confidence))
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
      print (df_edges %>% select(x,y,lag,type,infOrt,trueOrt,sign,log_confidence) )
      }
    #
    # Keep the rows having max log_confidence when grouped by on x, y
    #
    df_group <- df_edges %>% group_by(x,y) %>% top_n(n=1, wt = log_confidence)
    df_group <- as.data.frame (df_group)
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
        log_c <- one_edge$log_confidence
        #
        # Update the list of lags (value has been computed in lag matrix)
        #
        df_group[edge_idx,]$lag <- df_adj[node_x, node_y]
        #
        # Update orientation
        #
        cond_for_orient <- ( (df_edges[["x"]] == node_x) 
                           & (df_edges[["y"]] == node_y) 
                           & (df_edges[["log_confidence"]] != log_c) )
        #
        # If edge was unique (no different lag between the nodes), nothing to do
        #
        if (sum(cond_for_orient) == 0) 
          next
        #
        # Edge was having two or more lags
        #
        df_other <- df_edges[cond_for_orient,]
        for (other_idx in 1:nrow(df_other) )
          {
          other_edge <- df_other[other_idx,]
          #
          # Loop on infOrt and trueOrt
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
                       node_x, "-", node_y, " log=", log_c, 
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
                       node_x, "-", node_y, " log=", log_c,
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
                       node_x, "-", node_y, " log=", log_c,
                       " orient=", group_edge_orient, 
                       ", other edge orientation=", other_edge_orient, 
                       " => update to birectional=6", sep="") )
              df_group[edge_idx,][[col_to_update]] <- 6
              next            
              }
            if (DEBUG)
              print (paste ("Update ", col_to_update, " orientation: ", 
                     node_x, "-", node_y, " log=", log_c,
                     " orient=", group_edge_orient, 
                     ", other edge orientation=", other_edge_orient, 
                     " => nothing done", sep="") )
            }
          }
        }
      }
      
    df_edges <- df_group
    if (DEBUG)
      {
      print (paste ("after flatten_mode == ", flatten_mode) )
      print (df_edges %>% select(x,y,lag,type,infOrt,trueOrt,sign,log_confidence) )
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
      print (df_edges %>% select(x,y,type,infOrt,sign,log_confidence))
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
    print (miic_return$all.edges.summary %>% select(x,y,lag,type,infOrt,sign,log_confidence))
    }
  return (miic_return)
  }

