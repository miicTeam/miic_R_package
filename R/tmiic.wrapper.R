#*****************************************************************************
# Filename   : tmiic.wrapper.R                   Creation date: 24 march 2020
#
# Description: Data transformation of time series for miic
#
# Author     : Franck SIMON (fsimon.informaticien@wanadoo.fr)
#
# Changes history:
# - 24 march 2020 : initial version
#*****************************************************************************

#-----------------------------------------------------------------------------
# tmiic.transform_data_for_miic
#-----------------------------------------------------------------------------
#' tmiic.transform_data_for_miic
#'
#' @description
#' Reorganizes the data using the history to create lagged nodes and samples 
#' in a format usable by miic
#'
#' @details 
#' The function slices the input data according to the lag max argument tau :
#' Data are expected to be received in a 3 dimensional array [n_samples * n_nodes * n_time]
#' History us assumed to be time ordered from the oldest (first rows) to lastest (ending rows)
#' 
#' The number of nodes is increased and renamed on tau layers, ie with tau=2:
#' node1, node2 => node1_lag0, node2_lag0, node1_lag1, node2_lag1, node1_lag2, node2_lag2
#' 
#' Every timestep (until n_time - tau) is converted into a sample in the lagged graph:
#' 
#' Timestep Node & value   Node & value   => Sample  Node & value   Node & value
#'   t-2    node1_val(t-2) node2_val(t-2) =>    i    node1_lag2_val node2_lag2_val
#'   t-1    node1_val(t-1) node2_val(t-1) =>    i    node1_lag1_val node2_lag1_val 
#'   t      node1_val(t)   node2_val(t)   =>    i    node1_lag0_val node2_lag0_val 
#'   
#'   t-3    node1_val(t-3) node2_val(t-3) =>    i'   node1_lag2_val node2_lag2_val
#'   t-2    node1_val(t-2) node2_val(t-2) =>    i'   node1_lag1_val node2_lag1_val
#'   t-1    node1_val(t-1) node2_val(t-1) =>    i    node1_lag0_val node2_lag0_val 
#'   
#'   t-4    node1_val(t-4) node2_val(t-4) =>    i"   node1_lag2_val node2_lag2_val
#'   t-3    node1_val(t-3) node2_val(t-3) =>    i"   node1_lag1_val node2_lag1_val
#'   t-2    node1_val(t-2) node2_val(t-2) =>    i"   node1_lag0_val node2_lag0_val
#'   ...    .............. .............. ..    ..   .............. ..............
#' until n_time-tau is reached. The same process is applied to all input samples
#'
#' @param data_tab [a 2D or 3D array] 
#' An array of the time series for all nodes and samples of dimensions
#' [n_samples * n_nodes * n_time] for 3D or [n_nodes * n_time] for 2D.
#' When a 2D array is supplied, the number of samples is assumed to be 1.
#'
#' @param tau [an int]
#' An int representating the max lag
#' 
#' @return a 2D array of dimensions [ n_samples * (n_time - tau) ), 
#'                                    n_nodes * (tau+1) ]
#'                                  
tmiic.transform_data_for_miic <- function (data_tab, tau)
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

  if (DEBUG)
    {
    print (paste ("Nb samples  :", n_samples, sep="") )
    print (paste ("Nb nodes    :", n_nodes, sep="") )
    print (paste ("Nb timesteps:", n_time, sep="") )
    print (paste ("Tau max     :", tau, sep="") )
    print ("")
    print ("Input df: sample 1 over tau timesteps:")
    print (data_tab[1,,1:tau])
    print ("")
    print ("Input df: sample 1 over last tau timesteps:")
    print (data_tab[1,,(n_time - tau + 1):n_time])
    }
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
  
  #--------------------------------  
  # Test only keep 1 sample every x
  #--------------------------------  
  # print ("tab samples avt:")
  # print (dim(tab_lagged))
  # print (paste ("nb samples orig:", n_samples) )
  # new_nb_samples <- (n_samples * n_samples_from_history) %/% 10
  # new_tab_lagged2 <- array ( data=NA, dim=c (new_nb_samples, new_n_nodes),
  #                     dimnames=list(seq(1,new_nb_samples), new_list_nodes) )
  # for (i in 1:new_nb_samples)
  #   {
  #   new_tab_lagged2 [i,] <- tab_lagged[i * 10,]
  #   }
  # tab_lagged <- new_tab_lagged2
  # print ("tab samples ap:")
  # print (dim(tab_lagged))
  # print (tab_lagged[1,])
  # print (tab_lagged[new_nb_samples,])
  #--------------------------------  
  # End test only keep 1 sample every x
  #--------------------------------  
  #
  # return a dataframe as miic expects 
  #
  return ( as.data.frame(tab_lagged, stringsAsFactors = FALSE) )
  }
