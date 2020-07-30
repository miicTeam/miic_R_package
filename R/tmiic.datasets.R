#*****************************************************************************
# Filename   : tmiic.datasets.R                   Creation date: 24 march 2020
#
# Description: time series data sets generation
#
# Author     : Franck SIMON (fsimon.informaticien@wanadoo.fr)
#
# Changes history:
# - 24 march 2020 : initial version
#*****************************************************************************
#
# 5 predefined models are defined :
#
# Model 0 (debug)
#
# X1(t) = sample * 10000 + node * 1000 + timestep
# X2(t) = sample * 10000 + node * 1000 + timestep
# X3(t) = sample * 10000 + node * 1000 + timestep
# X4(t) = sample * 10000 + node * 1000 + timestep
# X5(t) = sample * 10000 + node * 1000 + timestep
#
# Model 1
#
# X1(t) = 0.20 * X1(t-1) + noise1(t)
# X2(t) = 0.60 * X2(t-1) - 0.287 * f( X5(t-1) ) + 0.287 * f( X4(t-1) ) + noise2(t)
# X3(t) = noise3(t) 
# X4(t) = 0.40 * X4(t-1) - 0.287 * f( X1(t-1) ) + 0.287 * f( X5(t-1) ) + noise4(t)
# X5(t) = 0.60 * X5(t-1) - 0.287 * f( X2(t-1) ) + noise5(t)
#
# Model 2
#
# X1(t) <- 0.40 * X1(t-1) + 0.287 * f( X5(t-1) ) + 0.287 * f( X4(t-1) ) + noise1(t)
# X2(t) <- 0.20 * X2(t-1) - 0.287 * f( X4(t-1) ) + noise2(t)
# X3(t) <- noise3(t)
# X4(t) <- 0.80 * X4(t-1) + 0.287 * f( X5(t-1) ) - 0.287 * f( X3(t-2) ) + noise4(t)
# X5(t) <- 0.60 * X5(t-1)  + noise5(t)
#
# Model 3
#
# X1(t) <- - 0.287 * f( X4(t-1) ) + noise1(t)
# X2(t) <-   0.40 * X2(t-1) - 0.287 * f( X4(t-1) ) + 0.287 * f( X1(t-1) ) +  noise2(t)
# X3(t) <-   0.90 * X3(t-1) - 0.287 * f( X1(t-2) ) + noise3(t)
# X4(t) <-   0.90 * X4(t-1) + noise4(t)
# X5(t) <-   0.80 * X5(t-1) - 0.287 * f( X3(t-2) )  + noise5(t)
#
# Model 4
#
# X1(t) <- 0.40 * X1(t-1) + noise1(t)
# X2(t) <- 0.20 * X2(t-1) - 0.287 * f( X1(t-1) ) + noise2(t)
# X3(t) <- noise3(t)
# X4(t) <- 0.80 * X4(t-1) + 0.287 * f( X1(t-1) ) - 0.287 * f( X2(t-1) ) + noise4(t)
# X5(t) <- 0.60 * X5(t-1) + 0.287 * f( X1(t-2) ) - 0.287 * f( X4(t) ) + noise5(t)
#
# Model 5 : test isolated edges orientation using time
#
# X1(t) <- 0.287 * X2(t-1) + noise1(t)
# X2(t) <- noise2(t)
# X3(t) <- 0.287 * X3(t-1) + noise3(t)
#
# Model 6 : orientation of non lagged edges with 4 points
#
# X1(t) <- noise1(t)
# X2(t) <- - 0.287 * X1(t-1) + noise2(t)
# X3(t) <- + 0.287 * X2(t)   - 0.287 * X4(t-1) + noise3(t)
# X4(t) <- noise4(t)
#
#
# Model 7 : latent variable on X3 (erased by new white noise after generation)
#
# X1(t) <- 0.40 * X1(t-1) + noise1[t]
# X2(t) <- 0.30 * X2(t-1) - 0.287 * f( X1(t-1) ) + 0.287 * f( X3(t-1) ) + noise2[t]
# X3(t) <- noise3[t]
# X4(t) <- 0.20 * X4(t-1) + 0.287 * f( X5(t-2) ) - 0.287 * f( X3(t-1) ) + noise4[t]
# X5(t) <- 0.60 * X5(t-1)  + noise5[t]
#
# noises ~ N(0, 1)
# f1 (x) = x
# f2 (x) = ( 1 - 4 * e[-(x^2)/2] ) * x
# f3 (x) = ( 1 - 4 * x^3 * e[-(x^2)/2] ) * x
#
#
# True graphs
#
tmiic.DF_TRUE_EDGES_MODEL_1 = data.frame (
  orig     = c(1   , 2   , 4   , 5   , 1     , 2     , 5     , 4     , 5     ), 
  dest     = c(1   , 2   , 4   , 5   , 4     , 5     , 2     , 2     , 4     ),
  lag      = c(1   , 1   , 1   , 1   , 1     , 1     , 1     , 1     , 1     ),
  strength = c(+0.2, +0.6, +0.4, +0.6, -0.287, -0.287, -0.287, +0.287, +0.287)  ,
  stringsAsFactors=FALSE)

tmiic.DF_TRUE_EDGES_MODEL_2 = data.frame (
  orig     = c(1   , 2   , 4   , 5   , 3     , 4     , 4     , 5     , 5    ), 
  dest     = c(1   , 2   , 4   , 5   , 4     , 1     , 2     , 1     , 4    ),
  lag      = c(1   , 1   , 1   , 1   , 2     , 1     , 1     , 1     , 1    ),
  strength = c(+0.4, +0.2, +0.8, +0.6, -0.287, +0.287, -0.287, +0.287,+0.287),
  stringsAsFactors=FALSE)

tmiic.DF_TRUE_EDGES_MODEL_3 = data.frame (
  orig     = c(2,    3,    4,    5,    1     , 1     , 3     , 4     , 4     ), 
  dest     = c(2,    3,    4,    5,    2     , 3     , 5     , 1     , 2     ),
  lag      = c(1   , 1   , 1   , 1   , 1     , 2     , 2     , 1     , 1     ),
  strength = c(+0.4, +0.9, +0.9, +0.8, +0.287, -0.287, -0.287, -0.287, -0.287),
  stringsAsFactors=FALSE)

tmiic.DF_TRUE_EDGES_MODEL_4 = data.frame (
  orig     = c(1   , 2,    4,    5,    1  ,    1  ,    2     , 1     , 4     ), 
  dest     = c(1   , 2,    4,    5,    2  ,    4  ,    4     , 5     , 5     ),
  lag      = c(1   , 1   , 1   , 1   , 1     , 1     , 1     , 2     , 0     ),
  strength = c(+0.4, +0.2, +0.8, +0.6, -0.287, +0.287, -0.287, +0.287, -0.287),
  stringsAsFactors=FALSE)

tmiic.DF_TRUE_EDGES_MODEL_5 = data.frame (
  orig     = c(2     , 3     ), 
  dest     = c(1     , 3     ),
  lag      = c(1     , 1     ),
  strength = c(+0.287, +0.287),
  stringsAsFactors=FALSE)

tmiic.DF_TRUE_EDGES_MODEL_6 = data.frame (
  orig     = c(1     , 2     , 4     ), 
  dest     = c(2     , 3     , 3     ),
  lag      = c(1     , 0     , 1     ),
  strength = c(-0.287, +0.287, -0.287),
  stringsAsFactors=FALSE)

tmiic.DF_TRUE_EDGES_MODEL_7 = data.frame (
  orig     = c(1   , 2   , 1     , 3     , 4   , 3     , 5     , 5   ), 
  dest     = c(1   , 2   , 2     , 2     , 4   , 4     , 4     , 5   ),
  lag      = c(1   , 1   , 1     , 1     , 1   , 1     , 2     , 1   ),
  strength = c(+0.4, +0.3, -0.287, +0.287, +0.2, -0.287, +0.287, +0.6),
  stringsAsFactors=FALSE)

#-----------------------------------------------------------------------------
# tmiic.f1
#-----------------------------------------------------------------------------
#' tmiic.f1
#'
#' @description 
#' basic linear function used for samples generation f(x) = x
#'
#' @param x [a number]
#'
#' @return the number received as argument 
#'
#' @export
#' @useDynLib miic
#-----------------------------------------------------------------------------
tmiic.f1 <- function (x) 
  {
  return (x)
  }

#-----------------------------------------------------------------------------
# tmiic.f2
#-----------------------------------------------------------------------------
#' tmiic.f2
#'
#' @description 
#' non linear function used for samples generation 
#' \eqn{f(x) = ( 1 - 4 * exp[-(x^2)/2] ) * x}
#'
#' @param x [a number]
#'
#' @return the number computed by f(x)
#' 
#' @export
#' @useDynLib miic
#-----------------------------------------------------------------------------
tmiic.f2 <- function (x) 
  {
  return ( ( 1 - 4 * exp(-(x^2)/2) ) * x )
  }

#-----------------------------------------------------------------------------
# tmiic.f3
#-----------------------------------------------------------------------------
#' tmiic.f3
#'
#' @description 
#' non linear function used for samples generation 
#' \eqn{f(x) = ( 1 - 4 * x^3 * exp[-(x^2)/2] ) * x}
#'
#' @param x [a number]
#'
#' @return the number computed by f(x)
#' 
#' @export
#' @useDynLib miic
#-----------------------------------------------------------------------------
tmiic.f3 <- function (x) 
  {
  return ( ( 1 - 4 * x^3 * exp(-(x^2)/2) ) * x )
  }

#-----------------------------------------------------------------------------
# Constant to find the model by its index
#-----------------------------------------------------------------------------
tmiic.LIST_MODELS = list (model1 = tmiic.DF_TRUE_EDGES_MODEL_1, 
                          model2 = tmiic.DF_TRUE_EDGES_MODEL_2, 
                          model3 = tmiic.DF_TRUE_EDGES_MODEL_3, 
                          model4 = tmiic.DF_TRUE_EDGES_MODEL_4, 
                          model5 = tmiic.DF_TRUE_EDGES_MODEL_5, 
                          model6 = tmiic.DF_TRUE_EDGES_MODEL_6, 
                          model7 = tmiic.DF_TRUE_EDGES_MODEL_7)

#-----------------------------------------------------------------------------
# tmiic.call_funct
#-----------------------------------------------------------------------------
#' tmiic.call_funct
#'
#' @description 
#' Utility function to call a function supplied as parameter. 
#'
#' @param funct [a function] The function to use. 
#' Typically, the function will be \emph{tmiic.f1}, \emph{tmiic.f2} or 
#' \emph{tmiic.f3} but any function accepting a float as parameter and 
#' returning a float can be used.
#'  
#' @return the number computed by the funct function
#-----------------------------------------------------------------------------
tmiic.call_funct <- function(funct, x)
  { 
  funct(x)
  }

#-----------------------------------------------------------------------------
# tmiic.generate_predefined_dataset
#-----------------------------------------------------------------------------
#' tmiic.generate_predefined_dataset
#'
#' @description 
#' Generate n_timeseries samples for x nodes over n_timesteps timesteps
#' using a predefined model
#'
#' @param model_idx [an integer betwen 0 and 7] The model for the data 
#' generation:
#' \itemize{
#' \item Model 0 is for debug (intialazing nodes with easily identiable values)
#' \item Model 1 uses only t-1 information. 
#' \item Models 2 and 3 use both t-1 and t-2.
#' \item Model 4 uses t-1, t-2 and contemporaneous edge
#' \item Model 5 isolated temporal edges oriented using time
#' \item Model 6 contemporaneous edge oriented using the 4 points condition
#' \item Model 7 latent node unveiled by edge between X2 and X4
#' }
#' @param funct [a function] the function used for the generation.\cr
#' It can be any function computing a float from a float or one of the
#' predefined ones:
#' \itemize{
#' \item \emph{tmiic.f1}: \eqn{f(x) = x}
#' \item \emph{tmiic.f2}: \eqn{f(x) = ( 1 - 4 * exp[-(x^2)/2] ) * x}
#' \item \emph{tmiic.f3}: \eqn{f(x) = ( 1 - 4 * x^3 * exp[-(x^2)/2] ) * x}
#' }
#' @param n_timeseries [an integer] The number of samples to generate
#' @param n_timesteps [an integer] The number of timesteps of the time series 
#' generated
#' @param seed [an integer] Optional, NULL by default. The seed to use when 
#' generating the dataset
#'
#' @return a list with two entries :
#' \itemize{
#' \item samples as an array of dimensions \emph{n_timeseries} * 
#' \emph{n_nodes} * \emph{n_timesteps}
#' \item a dataframe with the true edes of the network
#' }
#' 
#' @export
#' @useDynLib miic
#-----------------------------------------------------------------------------
tmiic.generate_predefined_dataset <- function (model_idx, funct, n_timeseries, 
                                               n_timesteps, seed=NULL) 
  {
  DEBUG <- FALSE
  if (DEBUG)
    {
    print ("tmiic.generate_predefined_dataset:")
    print (paste ("model_idx=", model_idx, sep="") )
    print (paste ("n_timeseries=", n_timeseries, sep="") )
    print (paste ("n_timesteps=", n_timesteps, sep="") )
    }
  
  if ( (model_idx < 0) | (model_idx > length(tmiic.LIST_MODELS)) )
    {
    stop (paste ("Predefined model", model_idx, "doesn't exist") )
    }
  #
  # All the predefined models have 5 nodes
  #
  list_nodes = list ("X1", "X2", "X3", "X4", "X5")
  n_nodes = length(list_nodes)
  #
  # Model 0 is debug
  #
  if (model_idx == 0)
    {
    true_edges <- data.frame ( orig=character(), dest=character(),
                               lag=integer(), strength=double(), 
                               stringsAsFactors = FALSE)
    df_data <- data.frame (matrix(ncol=n_nodes+1, nrow=n_timeseries*n_timesteps), 
                             stringsAsFactors=FALSE)
    colnames (df_data) <- c("timestep", list_nodes)
  
    for(timseseries_idx in 1:n_timeseries)
      {
      for(timestep_idx in 1:n_timesteps)
        {
        new_line_idx <- (timseseries_idx-1) * n_timesteps + timestep_idx
        df_data [new_line_idx, 1] <- timestep_idx
        for(node_idx in 1:n_nodes)
          {
          df_data [new_line_idx, node_idx+1] <- timseseries_idx * 10000 + node_idx * 1000 + timestep_idx
          }
        }
      }
    }
  else
    {
    #
    # Model >= 1 are predefined models
    #
    true_edges <- tmiic.LIST_MODELS[[model_idx]]
    df_data <- tmiic.generate_dataset (true_edges, funct, list_nodes, 
                                       n_timeseries, n_timesteps, seed) 
    #
    # For model 7, erase 3rd node data which is the latent node and modify the 
    # true edges dataframe to remove edges starting from node3 and add a latent
    # edge between node 2 and 4 (indicated by strengh = 0)
    #
    if (model_idx == 7)
      {
      df_data[,4] <- rnorm (n_timeseries * n_timesteps)
      true_edges <- true_edges[true_edges$orig != 3,]
      true_edges [nrow(true_edges) + 1,] = c(2, 4, 0, 0)
      }
    }
  if (DEBUG)
    {
    print ("")
    print ("Returned true edges:")
    print (true_edges)
    print ("Returned array dimension:")
    print (dim(df_data))
    print ("Returned data for first sample=1, 5 first timestep:")
    print (df_data[1:5,])
    print (paste ("Returned array on last sample=", n_timeseries, ", 5 last timesteps", sep="") )
    last_row = nrow(df_data)
    print (df_data[(last_row-5):last_row,])
    }
  return ( list(df_data, true_edges) )
  }

#-----------------------------------------------------------------------------
# tmiic.generate_dataset
#-----------------------------------------------------------------------------
#' tmiic.generate_dataset
#'
#' @description 
#' Generate \emph{n_timeseries} samples for \emph{n_nodes} nodes over 
#' \emph{n_timesteps} timesteps
#'
#' @details 
#' The function uses the \emph{funct} and the \emph{true_edges} parameters 
#' to generate temporal dataset.
#'  
#' i.e.: if in the true edges, the node1 is:
#' \itemize{
#' \item connected to itself at t-1 with strengh 0.5 
#' \item connected to the node2 at t-2 with strength 0.2
#' } 
#' the generation will compute for each timestep of node 1:\cr
#' node1(t) <- 0.5 * funct ( node1[t-1] ) + 0.2 * funct ( node2[t-2] ) + white noise
#' 
#' Values of the first iterations are discarded, so the temporal data
#' generated are not affected by the initial values.
#' 
#' @param true_edges [a dataframe] The true_edges dataframe must at least 
#' contains the columns orig, dest, lag and strengh. The orig and dest 
#' columns are the indexes of the edges's nodes. Lag is a postive or 0 integer.
#' Strengh is a real number positive or negative indicating the impact of the
#' orig node to the dest one.\cr
#' CAUTION: if a 0 lag is used, the order of the edges in the input dataframe
#' is important. The values of orig node(s) of the lag 0 edge must be completly
#' computed before the dest node of the lag 0 edge is computed.\cr
#' i.e. : if we have edges X3(t-1) -> X2(t) and X2(t) -> X1(t)
#' the edge X3(t-1) -> X2(t) must appear first in the true_edges parameter.
#' @param funct [a function] The function to apply when an edge exist between
#' two nodes
#' @param list_nodes [a list] The list of nodes in the dataset
#' @param n_timeseries [an integer] The number of samples to generate
#' @param n_timesteps [an integer] The number of timesteps of the time series 
#' generated
#' @param seed [an integer] Optional, NULL by default. The seed to use when 
#' generating the dataset
#' 
#' @return an array of dimensions n_timeseries * n_nodes * n_timesteps
#' 
#' @export
#' @useDynLib miic
#-----------------------------------------------------------------------------
tmiic.generate_dataset <- function (true_edges, funct, list_nodes, n_timeseries, 
                                    n_timesteps, seed=NULL) 
  {
  DEBUG <- FALSE
  if (DEBUG)
    {
    print ("tmiic.generate_dataset:")
    print ("true edges=")
    print (true_edges)
    print ("list_nodes=")
    print (list_nodes)
    print (paste ("n_timeseries=", n_timeseries, sep="") )
    print (paste ("n_timesteps=", n_timesteps, sep="") )
    }
  n_nodes = length (list_nodes)
  max_lag = max(true_edges$lag)
  #
  # We will generate more timesteps than requested to discard the first ones
  # so the generated data are not affected by the intial values.
  #
  n_time_generation = n_timesteps * 2
  # 
  # Set seed if specified and init an array with white noise
  #
  if ( !is.null (seed) )
    set.seed (seed)
  df_random <- array ( data=rnorm(n_timeseries * n_nodes * n_time_generation),
                       dim=c(n_timeseries, n_nodes, n_time_generation),
                       dimnames=list ( seq(1,n_timeseries), list_nodes, seq(1,n_time_generation)
                     )               )
  #
  # Generate n_timeseries time series
  #
  df_data <- data.frame (matrix(ncol=n_nodes+1, nrow=0), stringsAsFactors=FALSE)
  colnames (df_data) <- c("timestep", list_nodes)
  
  for(timseseries_idx in 1:n_timeseries)
    {
    # Generate the time serie
    #
    df_tmp <- data.frame (matrix(ncol=n_nodes, nrow=n_time_generation), 
                          stringsAsFactors=FALSE)
    for (timestep_idx in 1:n_time_generation)
      {
      if ( (DEBUG) & (timseseries_idx == 1)  & (timestep_idx == max_lag+1) )
        {
        print ("Nodes vals before random init=") 
        print (df_tmp[timestep_idx,])
        print ("noise at t=") 
        print (df_random[timseseries_idx,,timestep_idx] ) 
        }
      df_tmp[timestep_idx,] <- df_random[timseseries_idx,,timestep_idx]
      #
      # Until max_lag, we just initialize the tmp df with the random values
      #
      if (timestep_idx <= max_lag)
        next          
      #
      # Above max_lag, we are able to compute the nodes values using their history
      #
      if ( (DEBUG) & (timseseries_idx == 1) & (timestep_idx == max_lag+1) )
        {
        print (paste ("timestep_idx=", timestep_idx, sep="") ) 
        print ("t-2=") 
        print (df_tmp[(timestep_idx-2),] ) 
        print ("t-1=") 
        print (df_tmp[(timestep_idx-1),] ) 
        print ("nodes vals at t (init with noise)=") 
        print (df_tmp[timestep_idx,] ) 
        }
      
      for (edge_idx in 1:nrow(true_edges) )
        {
        one_edge <- true_edges[edge_idx,]
        data_orig_idx <- timestep_idx - one_edge$lag
        #
        # Function is applied only when edge is from a different node
        #
        if (one_edge$orig == one_edge$dest)
          fct_res <- df_tmp[data_orig_idx, one_edge$orig]
        else
          fct_res <- tmiic.call_funct (funct, df_tmp[data_orig_idx, one_edge$orig])
        
        if ( (DEBUG) & (timseseries_idx == 1)  & (timestep_idx == max_lag+1) )
          {
          print ("edge=") 
          print (one_edge) 
          print (paste ("old val of dest node ", list_nodes[[one_edge$dest]], "=",
                        df_tmp[timestep_idx, one_edge$dest], sep="") )
          print (paste ("use val of orig node ", list_nodes[[one_edge$orig]], 
                        " lag ", one_edge$lag, "=",
                        df_tmp[data_orig_idx, one_edge$orig], sep="") )
          print (paste ("function result=", fct_res, sep="") )
          print (paste ("function result * strength=", one_edge$strength * fct_res, sep="") )
          }
        
        df_tmp[timestep_idx, one_edge$dest] <- df_tmp[timestep_idx, one_edge$dest] + one_edge$strength * fct_res

        if ( (DEBUG) & (timseseries_idx == 1) & (timestep_idx == max_lag+1) )
          {
          print (paste ("new val at t of dest node ", list_nodes[[one_edge$dest]], "=",
                        df_tmp[timestep_idx, one_edge$dest], sep="") ) 
          }
        }

      if ( (DEBUG) & (timseseries_idx == 1) & (timestep_idx == max_lag+1) )
        {
        print ("new val at t=") 
        print (df_tmp[timestep_idx,] ) 
        }
      }
    #
    # One timeseries has been fully computed in df_tmp, add last timesteps into data
    #
    if (DEBUG)
      {
      print ("Tmp array dimension")
      print ( dim(df_tmp) )
      print ("Tmp Data around future 1st row ")
      print (df_tmp[(n_time_generation-n_timesteps):(n_time_generation-n_timesteps+1),])
      dbg_row_data = nrow(df_data)      
      print ("Data array dimension")
      print ( dim(df_data) )
      print ("End row of data")
      print (df_data[dbg_row_data,])
      }
    df_tmp <- df_tmp[(n_time_generation-n_timesteps+1):n_time_generation,]
    df_tmp <- cbind (timestep=seq(1,n_timesteps), df_tmp)
    df_data <- rbind (df_data, df_tmp)
    if (DEBUG)
      {
      print ("New Tmp array dimension")
      print ( dim(df_tmp) )
      print ("Tmp Data 1st row ")
      print (df_tmp[1,])
      print ("Data array new dimension")
      print ( dim(df_data) )
      print ("Data around add of rows")
      print (df_data[dbg_row_data:(dbg_row_data+1),])
      }
    }
  row.names(df_data) <- NULL
  if (DEBUG)
    {
    print ("Returned array dimension")
    print ( dim(df_data) )
    print ("Returned array on first sample=1, 5 first timesteps")
    print (df_data [1:5,])
    print (paste ("Returned array on last sample=", n_timeseries, ", 5 last timesteps", sep="") )
    last_row <- nrow(df_data)
    print (df_data [(last_row-5):last_row,])
    }
  return (df_data)
  }

#-----------------------------------------------------------------------------
# tmiic.plot_one_timeseries
#-----------------------------------------------------------------------------
#' tmiic.plot_one_timeseries
#' 
#' @description 
#' Plot one sample of a temporal serie
#'
#' @param df_data [an array of dimensions n_timeseries * n_nodes * n_timesteps]
#' The array containing the samples
#' @param timseseries_idx [an integer] The sample index to plot
#' @param title [a string] Optional, NULL by default. The title of the plot
#' @param filename [a string] Optional, NULL by default. If supplied, the plot
#' is saved in this file
#'
#' @return None
#' 
#' @export
#' @useDynLib miic
#-----------------------------------------------------------------------------
tmiic.plot_one_timeseries <- function (df_data, timseseries_idx, 
                                       title=NULL, filename=NULL)
  {
  DEBUG <- FALSE
  if (DEBUG)
    {
    print ("tmiic.plot_one_timeseries")
    print ("df_data")
    print (head (df_data))
    print (paste ("timeseries to plot:", timseseries_idx, sep="") )
    }
  #
  # Find the timeseries
  #
  n_rows <- nrow(df_data)
  n_timeseries_found <- 1
  previous_row_idx <- 1
  previous_timestep <- -Inf
  for (row_idx in 1:n_rows)
    {
    timestep = df_data[row_idx,1]
    if (timestep < previous_timestep)
      {
      if (n_timeseries_found == timseseries_idx)
        {
        row_idx <- row_idx - 1
        break
        }
      n_timeseries_found <- n_timeseries_found + 1
      previous_row_idx <- row_idx
      }
    previous_timestep <- timestep
    }
  if (DEBUG)
    print (paste ("Timeseries search result: ", n_timeseries_found, sep="") )
  
  if (n_timeseries_found != timseseries_idx)
    {
    print (paste ("tmiic.plot_one_timeseries: timeseries ", timseseries_idx,
           " not found, nothing to plot...", sep="") )
    return ()
    }
  if (DEBUG)
    print (paste ("Timeseries found between ", previous_row_idx, " and ", row_idx, sep="") )
  
  list_nodes <- colnames(df_data)[-1]
  df_plot <- df_data[previous_row_idx:row_idx, list_nodes]
  if (DEBUG)
    {
    print ("head df to plot")
    print (head (df_plot))
    }
  #
  # Plot part
  #
  if (! is.null(filename) )
    {
    png (filename=filename)
    }
  plot.ts (df_plot)
  if (! is.null(title) )
    {
    title (title)
    }
  if (! is.null(filename) )
    {
    dev.off()
    }
  }

