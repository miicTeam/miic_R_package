#*****************************************************************************
# Filename   : tmiic.scores.R                   Creation date: 24 march 2020
#
# Description: Methods to evaluate the reconstruction of timed graphs
#
# Author     : Franck SIMON (fsimon.informaticien@wanadoo.fr)
#
# Changes history:
# - 24 march 2020 : initial version
#*****************************************************************************

#-----------------------------------------------------------------------------
# tmiic.compute_scores
#-----------------------------------------------------------------------------
#' tmiic.compute_scores
#'
#' @description 
#' Compute confusion matrix, precision, recall and f-score on the 
#' reconstructed graph. Values are computed both for unoriented and 
#' oriendted graphs
#'
#' @param df_edges [a dataframe] The list of edges computed by tmiic
#' (the summary dataframe of tmiic's return)
#' @param df_true_edges [a dataframe]  The list of true edges
#' @param list_nodes [a list] The list of nodes in the network
#' @param n_layers [an integer] The number of layers in the temporal graph
#' 
#' @returns [a list] The returned list contains two data frames containing 
#' both the confusion matrix, the precision, recall and f-score. The first 
#' dataframe is for unoriented edges and the second about oriented ones.
#'
#' @export
#' @useDynLib miic
#-----------------------------------------------------------------------------
tmiic.compute_scores <- function (df_edges, df_true_edges, list_nodes, n_layers) 
  {
  DEBUG <- FALSE
  if (DEBUG)
    {
    print ("tmiic.compute_scores:")
    print ("List nodes:")
    print (list_nodes)
    print (paste ("n_layers:", n_layers, sep="") )
    print (paste ("Input edges df (class ", class(df_edges), ") limited to type 'P'=", sep="") )
    print (df_edges[df_edges$type=='P',] 
             %>%  dplyr::select (x,y,type,infOrt,sign,info,info_cond,cplx,Nxy_ai,log_confidence))
    print (paste ("True edges df (class ", class(df_true_edges), ")", sep="") )
    print (df_true_edges)
    }
  n_nodes <- length(list_nodes)
  #
  # If true edge list is not supplied (ie debug model), ends
  #
  df_empty <- data.frame ( matrix(ncol = 7, nrow = 0), stringsAsFactors=FALSE)
  colnames (df_empty) <- c("tp", "fp", "fn", "tn", "precision", "recall", "fscore")
  if (is.null(df_true_edges))
    return (list (df_empty, df_empty) )
  if (nrow(df_true_edges) <= 0)
    return (list (df_empty, df_empty) )
  #
  # Evaluate true/false positive
  #
  fp <- 0
  tp <- 0
  fp_orient <- 0
  tp_orient <- 0
  for (row_idx in 1:nrow(df_edges) )
    {
    one_edge <- df_edges[row_idx,]
    #
    # if type != P, this is not a real edge in fact
    #
    if (one_edge$type != "P")
      {
      if (DEBUG)
        {
        print (paste(one_edge$x, "-", one_edge$y, ", type=", one_edge$type, 
                     " => type!='P', don't consider edge", sep="") )
        }
      next  
      }
    node_x <- one_edge$x
    node_y <- one_edge$y
    orient <- one_edge$infOrt
    #
    # Get the lag for each the vertices of the edge
    #
    pos_lag_x <- stringr::str_locate(node_x, "_lag")
    lag_x <- 0
    if ( !is.na(pos_lag_x[1]) )
      {
      lag_x <- stringr::str_remove(node_x, ".*_lag")
      lag_x <- strtoi (lag_x)
      }
    pos_lag_y <- stringr::str_locate(node_y, "_lag")
    lag_y <- 0
    if ( !is.na(pos_lag_y[1]) )
      {
      lag_y <- stringr::str_remove(node_y, ".*_lag")
      lag_y <- strtoi (lag_y)
      }
    #
    # It temporal, ensure to order from oldest to newest
    #
    diff = lag_x - lag_y
    if (diff < 0)
      {
      diff <- - diff
      
      temp_var <- lag_x
      lag_x <- lag_y
      lag_y <- temp_var
      
      temp_var <- node_x
      node_x <- node_y
      node_y <- temp_var
      
      if ( (orient == 2) | (orient == -2) )
        orient = -orient
      }
    #
    # Remove "_lag" information from the nodes names
    # find the index of the nodes in the nodes list
    #
    node_oldest <- gsub("_lag.*","",node_x)
    node_newest <- gsub("_lag.*","",node_y)
    idx_node_oldest <- which (list_nodes == node_oldest)[[1]]
    idx_node_newest <- which (list_nodes == node_newest)[[1]]
    #
    # Search this edge in the true edge df 
    #
    cond_old_new = (  (df_true_edges[["orig"]] == idx_node_oldest) 
                    & (df_true_edges[["dest"]] == idx_node_newest) 
                    & (df_true_edges[["lag"]] == diff) )
    #
    # If edge is not temporal, we must look also the opposite way
    #
    cond_new_old = 0
    if (diff == 0)
      cond_new_old = (  (df_true_edges[["orig"]] == idx_node_newest) 
                     & (df_true_edges[["dest"]] == idx_node_oldest) 
                     & (df_true_edges[["lag"]] == 0) )
    
    if ( (sum(cond_old_new) == 0) && (sum(cond_new_old) == 0) )
      {
      #
      # If edge is not in the true edges list => false pos
      #
      fp <- fp + 1
      
      if (orient == 1)
        {
        #
        # Edge not in the true edges list and was not oriented => Nothing to count for oriented
        #
        if (DEBUG)
          print (paste (one_edge$x, "-", one_edge$y, " type=", one_edge$type, " ort=", one_edge$infOrt, 
                        " => ", node_oldest, "-", node_newest, " lag=", diff, 
                        "  => not found and unoriented => fp + 1, fp_orient unchanged", sep="") )
        next
        }
      #
      # If the edge not in the true edges list is oriented => false pos also for oriented
      #
      if (orient == 6) # bidirectional
        {
        fp_orient <- fp_orient + 2
        if (DEBUG)
          print (paste (one_edge$x, "-", one_edge$y, " type=", one_edge$type, " ort=", one_edge$infOrt, 
                        " => ", node_oldest, "-", node_newest, " lag=", diff, 
                        " => not found and oriented => fp + 1, fp_orient + 2", sep="") )
        }
      else
        {
        fp_orient <- fp_orient + 1
        if (DEBUG)
          print (paste (one_edge$x, "-", one_edge$y, " type=", one_edge$type, " ort=", one_edge$infOrt, 
                        " => ", node_oldest, "-", node_newest, " lag=", diff, 
                        " => not found and oriented => fp + 1, fp_orient + 1", sep="") )
        }
      next
      }
    #
    # The edge is in the true edge list => at least true pos for non oriented
    #
    tp <- tp + 1
    if (orient == 1)
      {
      #
      # Edge computed has no orientation, nothing to count for oriented
      #
      if (DEBUG)
        print (paste (one_edge$x, "-", one_edge$y, " type=", one_edge$type, " ort=", one_edge$infOrt, 
                      " => ", node_oldest, "-", node_newest, " lag=", diff, 
                      " => found and unoriented => tp + 1, tp_orient unchanged",  sep="") )
      next
      }
    #
    # As latent variable are identified in the true edges by strength = 0
    # Set new contion taking excluding latent variable (because they are unoriented)
    #
    cond_old_new_not_latent = ( (df_true_edges[["orig"]] == idx_node_oldest) 
                              & (df_true_edges[["dest"]] == idx_node_newest) 
                              & (df_true_edges[["lag"]] == diff) 
                              & (df_true_edges[["strength"]] != 0) )
    #
    # If edge is not temporal, we must look also the opposite way
    #
    cond_new_old_not_latent = 0
    if (diff == 0)
      cond_new_old_not_latent = ( (df_true_edges[["orig"]] == idx_node_newest) 
                                & (df_true_edges[["dest"]] == idx_node_oldest) 
                                & (df_true_edges[["lag"]] == 0)  
                                & (df_true_edges[["strength"]] != 0) )
    
    #
    # Edge computed has an orientation, check if correct
    #
    if (orient == 6)
      {
      # 
      # Edge commputed is birectional => we must find 2 true edges not latent
      #
      if (sum(cond_old_new_not_latent) + sum(cond_new_old_not_latent) == 2)
        {
        #
        # Computed edge (implies not temporal) found in the 2 ways => + 2 true positive
        #
        if (DEBUG)
          print (paste (one_edge$x, "-", one_edge$y, " type=", one_edge$type, " ort=", one_edge$infOrt, 
                        " => ", node_oldest, "-", node_newest, " lag=", diff, 
                        " => found 2 times (bidirectional)",
                        " => tp + 1, tp_orient + 2", sep="") )
        tp_orient <- tp_orient + 2
        next
        }
      #
      # Computed edge found in 1 way only => + 1 true positive, + 1 false positive
      #
      if (sum(cond_old_new_not_latent) + sum(cond_new_old_not_latent) == 1)
        {
        if (DEBUG)
          print (paste (one_edge$x, "-", one_edge$y, " type=", one_edge$type, " ort=", one_edge$infOrt, 
                        " => ", node_oldest, "-", node_newest, " lag=", diff, 
                        " => found 1 time (unidirectional)",
                        " => tp + 1, tp_orient + 1, fp_orient + 1", sep="") )
        tp_orient <- tp_orient + 1
        fp_orient <- fp_orient + 1
        next
        }
      #
      # Computed edge not found in non latent  => + 2 false positive
      #
      if (DEBUG)
        print (paste (one_edge$x, "-", one_edge$y, " type=", one_edge$type, " ort=", one_edge$infOrt, 
                      " => ", node_oldest, "-", node_newest, " lag=", diff, 
                      " => found but latent (not oriented)",
                      " => tp + 1, fp_orient + 2", sep="") )
      fp_orient <- fp_orient + 2
      next
      }
    #
    # Computed edge is unidirectional (orient == 2)
    # and has been found including latent
    # => look if still good when latent are excluded
    #
    if (sum(cond_old_new_not_latent) + sum(cond_new_old_not_latent) == 0)
      {
      if (DEBUG)
        print (paste (one_edge$x, "-", one_edge$y, " type=", one_edge$type, " ort=", one_edge$infOrt, 
                      " => ", node_oldest, "-", node_newest, " lag=", diff, 
                      " => found but latent (not oriented) => tp + 1, fp_orient + 1", sep="") )
      fp_orient <- fp_orient + 1
      next
      }
    #
    # Computed edge is unidirectional (orient == 2)
    # and has been found excluding latent
    # => look if correctly oriented
    #
    if (sum(cond_old_new_not_latent) > 0)
      {
      if (DEBUG)
        print (paste (one_edge$x, "-", one_edge$y, " type=", one_edge$type, " ort=", one_edge$infOrt, 
                      " => ", node_oldest, "-", node_newest, " lag=", diff, 
                      " => found and correctly oriented => tp + 1, tp_orient + 1", sep="") )
      tp_orient <- tp_orient + 1
      next
      }
    #
    # Edge found but badly oriented
    #
    if (DEBUG)
      print (paste (one_edge$x, "-", one_edge$y, " type=", one_edge$type, " ort=", one_edge$infOrt, 
                    " => ", node_oldest, "-", node_newest, " lag=", diff, 
                    " => found and badly oriented => tp + 1, fp_orient + 1", sep="") )
    fp_orient <- fp_orient + 1
    }
  #
  # Compute the number of true edges
  #
  nb_true_edges_total <- nrow(df_true_edges)
  nb_true_edges_orient <- nrow(df_true_edges[df_true_edges$strength != 0,])
  nb_true_edges_non_orient <- nb_true_edges_total
  for (row_idx in 1:nrow(df_true_edges) )
    {
    # Exclude lagged edges as they appear always once
    one_edge <- df_true_edges[row_idx,]
    if (one_edge$lag != 0)
      next
    # Exclude latent edges as they appear always once
    if (one_edge$strength != 0) 
      next
    cond_opposite = (  (df_true_edges[["orig"]] == one_edge$dest) 
                     & (df_true_edges[["dest"]] == one_edge$orig) 
                     & (df_true_edges[["lag"]] == 0) )
    #
    # If we have opposite edges at lag=0 between the same edge,
    # then we must decrease the number of unoriented by 1
    # => remove 0.5 because we will find a 2 matches
    #
    if (sum(cond_opposite) > 0)
      nb_true_edges_non_orient <- nb_true_edges_non_orient - 0.5
    }
  
  tp_rate <- tp / nb_true_edges_non_orient
  tp_orient_rate <- tp_orient / nb_true_edges_orient
  #
  # Compute nb of no edges 
  #
  nb_possible_edges_for_one_lag0_node <- (n_layers-1) + n_layers * (n_nodes - 1)
  nb_possible_edges_with_one_node_lag0 <- nb_possible_edges_for_one_lag0_node * n_nodes
  nb_possible_edges_non_orient <- nb_possible_edges_with_one_node_lag0 - (n_nodes - 1)
  nb_possible_edges_orient <- nb_possible_edges_non_orient * 2
  
  nb_no_edges_non_orient <- nb_possible_edges_non_orient - nb_true_edges_non_orient
  nb_no_edges_orient <- nb_possible_edges_orient - nb_true_edges_orient
  fp_rate <- fp / nb_no_edges_non_orient 
  fp_orient_rate <- fp_orient / nb_no_edges_orient
  #
  # Evaluate false negative
  #
  fn <- 0
  fn_orient <- 0
  for (row_idx in 1:nrow(df_true_edges) )
    {
    one_edge <- df_true_edges[row_idx,]
    #
    # Get nodes the edge
    #
    node_orig <- list_nodes[[one_edge$orig]] 
    node_dest <- list_nodes[[one_edge$dest]]
    lag <- one_edge$lag
    strength <- one_edge$strength
    node_orig <- paste(node_orig, "_lag", lag, sep="")
    node_dest <- paste(node_dest, "_lag0", sep="")
    #
    # Search if this edge has been found
    #
    cond_old_new = (  (df_edges[["x"]] == node_orig) 
                    & (df_edges[["y"]] == node_dest)
                    & (df_edges[["type"]] == "P") )
    cond_new_old = (  (df_edges[["x"]] == node_dest) 
                    & (df_edges[["y"]] == node_orig)
                    & (df_edges[["type"]] == "P") )
    #
    # If edge is not in the computed edges list => false neg
    #
    if ( (sum(cond_old_new) == 0) & (sum(cond_new_old) == 0) )
      {
      fn <- fn +1
      if (strength != 0)
        {
        fn_orient <- fn_orient + 1
        if (DEBUG)
          print (paste (one_edge$orig, "-", one_edge$dest, " lag=", lag, " strength=", strength, 
                        " => ", node_orig, "-", node_dest, 
                        " not found in computed edges => fn + 1, fn_orient + 1", sep="") )
        }
      else
        {
        if (DEBUG)
          print (paste (one_edge$orig, "-", one_edge$dest, " lag=", lag, " strength=", strength, 
                        " => ", node_orig, "-", node_dest, 
                        " not found in computed edges => fn + 1, not fn_orient because latent", sep="") )
        }
      next
      }
    #
    # We have an edge matching, it is not a FN for unoriented
    #
    if (DEBUG)
      print (paste (one_edge$orig, "-", one_edge$dest, " lag=", lag, 
                    " => ", node_orig, "-", node_dest, 
                    " => found (not an unoriented FN)", sep="") )
    #
    # If true edge is a latent variable, it is unoriented 
    # => Nothing to do about orientation (can not be a FN oriented)
    #
    if (strength == 0)
      {
      if (DEBUG)
        print (paste (one_edge$orig, "-", one_edge$dest, " lag=", lag, " strength=", strength, 
                      " => ", node_orig, "-", node_dest, 
                      " => found and true edge is latent => noting to do on FN oriented", sep="") )
      next
      }
    #    
    # We have an edge matching the true edge and the true edge is oriented
    # => look for orientation
    #
    # Normally we should not find two edges with the nodes (we can find the edge with 
    # the same starting and ending nodes or the opposite but normally, not the two)
    #
    if ( (sum(cond_old_new) > 0) & (sum(cond_new_old) > 0) )
      {
      message ("Warning: tmiic compute score, 2 edges found : should not happen")
      if (  ( (df_edges[cond_old_new,]$infOrt != 6) & (df_edges[cond_old_new,]$infOrt != 2) )
          & ( (df_edges[cond_new_old,]$infOrt != 6) & (df_edges[cond_new_old,]$infOrt != -2) ) )
        {
        fn_orient <- fn_orient + 1
        if (DEBUG)
          print (paste (one_edge$orig, "-", one_edge$dest, " lag=", lag, 
                        " => ", node_orig, "-", node_dest, 
                        " => found 2 times but wrong orient => fn_orient + 1", sep="") )

        }
      next
      }
    #
    # We know that one edge and only one has been found
    #
    if ( (sum(cond_old_new) > 0) )
      {
      orient <- df_edges[cond_old_new,]$infOrt
      if ( (orient != 6) & (orient != 2) )
        {
        fn_orient <- fn_orient + 1
        if (DEBUG)
          print (paste (one_edge$orig, "-", one_edge$dest, " lag=", lag, 
                        " => ", node_orig, "-", node_dest, 
                        " found but wrong orient=", orient, " => fn_orient + 1", sep="") )
        }
      else
        {
        if (DEBUG)
          print (paste (one_edge$orig, "-", one_edge$dest, " lag=", lag, 
                        " => ", node_orig, "-", node_dest, 
                        " => found with good orient=", orient, 
                        " => Not a FN oriented", sep="") )
        }
      next
      }
    #
    # We know there is an edge and it is not the one with the same starting end ending node
    #
    orient <- df_edges[cond_new_old,]$infOrt
    if ( (orient != 6) & (orient != -2) )
      {
      fn_orient <- fn_orient + 1
      if (DEBUG)
        print (paste (one_edge$orig, "-", one_edge$dest, " lag=", lag, 
                      " => ", node_orig, "-", node_dest, 
                      " => found but wrong orient=", orient, " => fn_orient + 1", sep="") )
      }
    else
      {
      if (DEBUG)
        print (paste (one_edge$orig, "-", one_edge$dest, " lag=", lag, 
                      " => ", node_orig, "-", node_dest, 
                      " => found with good orient=", orient, " => Not a FN oriented", sep="") )
      }
    }
  fn_rate <- fn / nb_true_edges_non_orient
  fn_orient_rate <- fn_orient / nb_true_edges_orient
  #
  # True negatif
  #
  tn <- nb_possible_edges_non_orient - tp - fp - fn
  tn_rate <- tn / nb_no_edges_non_orient

  tn_orient <- nb_possible_edges_orient - tp_orient - fp_orient - fn_orient
  tn_orient_rate <- tn_orient / nb_no_edges_orient
  #
  # Precision, Recall and F-Score
  #
  val_precision <- 0
  if (tp + fp > 0)
    {
    val_precision <- tp / (tp + fp)
    }
  val_recall <- 0
  if (tp + fn > 0)
    {
    val_recall <- tp / (tp + fn)
    }
  val_fscore <- 0
  if ( (2 * tp) + fp + fn > 0 )
    {
    val_fscore <- (2 * tp) / ( (2 * tp) + fp + fn)
    }

  val_precision_orient <- 0
  if (tp_orient + fp_orient > 0)
    {
    val_precision_orient <- tp_orient / (tp_orient + fp_orient)
    }
  val_recall_orient <- 0
  if (tp_orient + fn_orient > 0)
    {
    val_recall_orient <- tp_orient / (tp_orient + fn_orient)
    }
  val_fscore_orient <- 0
  if ( (2 * tp_orient) + fp_orient + fn_orient > 0 )
    {
    val_fscore_orient <- (2 * tp_orient) / ( (2 * tp_orient) + fp_orient + fn_orient)
    }
  #
  # Trace and return results
  #
  ret <- data.frame (tp=tp_rate, fp=fp_rate, fn=fn_rate, tn=tn_rate, 
                     precision=val_precision, recall=val_recall, fscore=val_fscore,
                     stringsAsFactors=FALSE)
  ret_orient <- data.frame (tp=tp_orient_rate, fp=fp_orient_rate, fn=fn_orient_rate, tn=tn_orient_rate, 
                           precision=val_precision_orient, recall=val_recall_orient, fscore=val_fscore_orient,
                           stringsAsFactors=FALSE)
  if (DEBUG)
    {
    print (paste ("tp ", tp, " fp ", fp, " fn ", fn, " tn ", tn, sep="") )
    print (paste ("nb true edges ", nb_true_edges_non_orient, " nb no edges ", nb_no_edges_non_orient, 
            " nb tot edges possible ", nb_possible_edges_non_orient, sep="") )
    print (paste ("tp rate ", tp_rate, " fp rate ", fp_rate, 
                  " fn rate ", fn_rate, " tn rate ", tn_rate, 
                  " prec ", val_precision, " recall ", val_recall, " fscore ", val_fscore, 
                  sep="") )
    print ("")
    print (paste ( "tp orient ", tp_orient, " fp orient ", fp_orient, 
                  " fn orient ", fn_orient, " tn orient ", tn_orient, sep="") )
    print (paste ("nb true edges ", nb_true_edges_orient, " nb no edges orient ", nb_no_edges_orient, 
                  " nb tot edges possible orient ", nb_possible_edges_orient, sep="") )
    print (paste ("tp orient rate ", tp_orient_rate, " fp orient rate ", fp_orient_rate, 
                  " fn orient rate ", fn_orient_rate, " tn orient rate ", tn_orient_rate, 
                  " prec orient ", val_precision_orient, " recall orient ", val_recall_orient, 
                  " fscore orient ", val_fscore_orient, sep="") )
    print ("")
    print ("Returned df:")
    print (ret)
    print ("Returned orient df:")
    print (ret_orient)
    }
  return ( list (ret, ret_orient) )
  }

#-----------------------------------------------------------------------------
# tmiic.compute_list_scores
#-----------------------------------------------------------------------------
#' tmiic.compute_list_scores
#'
#' @description 
#' Compute scores on a list of tmiic results
#'
#' @param list_res [a list] the list of tmiic's results. It can be a list of
#' full results or only summary
#' @param df_true_edges [a dataframe] The list of the true edges
#' @param list_nodes [a list] The list of nodes in the network
#' @param n_layers [an integer] The number of layers in the temporal graph
#' 
#' @returns [a list] the returned list contains two dataframes containing 
#' both the list of scores for each tmiic result. The first dataframe is for
#' unoriented edges and the second about oriented ones.
#'
#' @export
#' @useDynLib miic
#-----------------------------------------------------------------------------
tmiic.compute_list_scores <- function (list_res, df_true_edges, list_nodes, n_layers) 
  {
  DEBUG <- FALSE
  if (DEBUG)
    {
    print ("tmiic.compute_list_scores:")
    print (paste ("input results list length:", length(list_res), sep="") )
    print ("input df_true_edges:")
    print (df_true_edges)
    }
  #
  # prepare df to  store performance evalutation
  #
  df_no_result <- data.frame (tp=NaN, fp=NaN, fn=NaN, tn=NaN, 
                              precision=NaN, recall=NaN, fscore=NaN, stringsAsFactors=FALSE)
  
  df_scores_unoriented <- data.frame ( matrix(ncol = 7, nrow = 0) )
  colnames (df_scores_unoriented) <- c("tp", "fp", "fn", "tn", "precision", "recall", "fscore")
  
  df_scores_oriented <- data.frame ( matrix(ncol = 7, nrow = 0) )
  colnames (df_scores_oriented) <- c("tp", "fp", "fn", "tn", "precision", "recall", "fscore")
  #
  # If true edge list is not supplied (ie debug model), ends
  #
  if (is.null(df_true_edges))
    return (list (df_scores_unoriented, df_scores_oriented) )
  if (nrow(df_true_edges) <= 0)
    return (list (df_scores_unoriented, df_scores_oriented) )
  #
  # Iterate over results
  # 
  for (list_idx in 1:length(list_res) )
    {
    one_res <- list_res[[list_idx]]
    if (all(is.na (one_res)))
      {
      df_scores_unoriented[nrow(df_scores_unoriented) + 1,] <- df_no_result
      df_scores_oriented  [nrow(df_scores_oriented)   + 1,] <- df_no_result
      }
    else
      {
      if (! is.data.frame(one_res))
        one_res <- one_res$all.edges.summary
      tmp_ret <- tmiic.compute_scores (one_res, df_true_edges, list_nodes, n_layers)
      df_scores_unoriented[nrow(df_scores_unoriented) + 1,] <- tmp_ret[[1]]
      df_scores_oriented  [nrow(df_scores_oriented)   + 1,] <- tmp_ret[[2]]
      }
    }
  
  if (DEBUG)
    {
    print ("Returned scores non oriented:")
    print (df_scores_unoriented)
    print ("Returned scores oriented:")
    print (df_scores_oriented)
    }
  return ( list(df_scores_unoriented, df_scores_oriented) )
  }

#-----------------------------------------------------------------------------
# tmiic.plot_scores
#-----------------------------------------------------------------------------
#' tmiic.plot_scores
#'
#' @description 
#' plot evolution of precision, recall and f-score or TP, FP, FN, TN rates
#'
#' @param scores [a dataframe] A dataframe of scores. Expected columns are
#' tp, fp, tn, fn, precision, recall and fscore  whilst each row correspond 
#' to a run of tmiic.
#' @param title [a string] The tille of the plot
#' @param list_labels [a list] The list of labels on the X axis. Typically,
#' these labels are the settings used for the different runs.
#' @param type [a string ] Optional, \emph{"PRFS"} by default. 
#' Represents the type of the plot: \emph{"PRFS"} plots Precision, Recall 
#' and F-Score whilst the other type available \emph{"TPFP"} plots 
#' True/False positive rates. 
#' @param plot_fntn [a boolean] Optional, FALSE by default. Applies only to  
#' type \emph{"TPFP"} plot. If TRUE, displays also True and False negative
#' rates.
#' @param filename [a string] Optional, NULL by default. The file name where
#' to save the plot,
#' 
#' @return None
#'
#' @export
#' @useDynLib miic
#-----------------------------------------------------------------------------
tmiic.plot_scores <- function(scores, title, list_labels, 
                              type="PRFS", plot_fntn=FALSE, filename=NULL)
  {
  DEBUG <- FALSE
  if (DEBUG)
    {
    print ("tmiic.plot_scores:")
    print ("Input scores matrix df:")
    print (scores)
   }
  #
  # If df score is not supplied or empty (ie debug model), ends
  #
  if (is.null(scores))
    return ()
  if (nrow(scores) <= 0)
    return ()
  #
  # Define graphic output
  #
  if (! is.null(filename) )
    {
    png (filename=filename, res=100)
    }
  #
  # Plot score
  #
  if (type == "PRFS")
    {
    plot (scores$precision, type="l", xlab="Number of samples", ylab="Rate", 
          main=title, cex.axis = 0.7, xaxt = 'n', ylim=c(0,1), col="blue")
    lines (scores$recall, col="red")
    lines (scores$fscore, col="black")
    legend(x=1, y=1, legend=c("Precision", "Recall", "F-Score"), col=c("blue","red","black"), lty=c(1,1,1), cex=0.5)
    }
  else  
    {
    plot (scores$tp, type="l", xlab="Number of samples", ylab="Rate", 
          main=title, cex.axis = 0.7, xaxt = 'n', ylim=c(0,1), col="blue")
    lines (scores$fp, col="red")
    if (plot_fntn)
      {
      lines (scores$fn, col="orange")
      lines (scores$tn, col="black")
      legend(x=1, y=1, legend=c("TP", "FP", "FN", "TN"), col=c("blue","red","orange", "black"), 
             lty=c(1,1), cex=0.5)
      }
    else
      {
      legend(x=1, y=1, legend=c("TP", "FP"), col=c("blue","red"), lty=c(1,1), cex=0.5)
      }
    }
  axis (side=1, cex.axis=0.7, at=seq ( 1:length(list_labels) ), labels=list_labels)
  
  if (! is.null(filename) )
    {
    dev.off()
    }
  }

