#' MIIC, causal network learning algorithm including latent variables
#'
#' @description MIIC (Multivariate Information based Inductive Causation) combines
#' constraint-based and information-theoretic approaches to disentangle direct
#' from indirect effects amongst correlated variables, including cause-effect
#' relationships and the effect of unobserved latent causes.
#'
#' @details In standard mode, starting from a complete graph, the method iteratively removes
#' dispensable edges, by uncovering significant information contributions from
#' indirect paths, and assesses edge-specific confidences from randomization of
#' available data. The remaining edges are then oriented based on the signature
#' of causality in observational data.
#'
#' In temporal mode, miic reorganizes the dataset using the \emph{n_layers} and
#' \emph{delta_t} parameters to transform the time steps into lagged samples.
#' As starting point, a lagged graph is created with only edges having at
#' least one node laying on the last time step.
#' Then, miic standard algorithm is applied to remove dispensable edges.
#' The remaining edges are then oriented by using the temporality and the
#' signature of causality in observational data. The use of temporal mode
#' is exposed in Simon \emph{et al.}, eLife reviewed preprint 2024.
#'
#' The method relies on an information theoretic based (conditional) independence
#' test which is described in (Verny \emph{et al.}, PLoS Comp. Bio. 2017),
#' (Cabeli \emph{et al.}, PLoS Comp. Bio. 2020). It deals with both categorical
#' and continuous variables by performing optimal context-dependent discretization.
#' As such, the input data frame may contain both numerical columns which will be
#' treated as continuous, or character / factor columns which will be treated
#' as categorical. For further details on the optimal discretization method and
#' the conditional independence test, see the function discretizeMutual.
#' The user may also choose to run miic with scheme presented in
#' (Li \emph{et al.}, NeurIPS 2019) to improve the end result's interpretability
#' by ensuring consistent separating set during the skeleton iterations.
#'
#' @seealso \code{\link{discretizeMutual}} for optimal discretization and
#' (conditional) independence test.
#'
#' @references
#' \itemize{
#' \item{Simon et al., \emph{eLife reviewed preprint}  https://www.biorxiv.org/content/10.1101/2024.02.06.579177v1.abstract }
#' \item{Ribeiro-Dantas et al., \emph{iScience 2024}  https://arxiv.org/abs/2303.06423 }
#' \item{Cabeli et al., \emph{PLoS Comp. Bio. 2020.}  https://doi.org/10.1371/journal.pcbi.1007866 }
#' \item{Li et al., \emph{NeurIPS 2019} http://papers.nips.cc/paper/9573-constraint-based-causal-structure-learning-with-consistent-separating-sets.pdf }
#' \item{Verny et al., \emph{PLoS Comp. Bio. 2017.}  https://doi.org/10.1371/journal.pcbi.1005662 }
#' }
#'
#' @param input_data [a data frame, required]
#'
#' A n*d data frame (n samples, d variables) that contains the observational data.
#'
#' In standard mode, each column corresponds to one variable and each row is a
#' sample that gives the values for all the observed variables.
#' The column names correspond to the names of the observed variables.
#' Numeric columns with at least 5 distinct values will be treated as continuous
#' by default whilst numeric columns with less than 5 distinct values, factors
#' and characters will be considered as categorical.
#'
#' In temporal mode, the expected data frame layout is variables as columns
#' and time series/time steps as rows.
#' The time step information must be supplied in the first column and,
#' for each time series, be consecutive and in ascending order (increment of 1).
#' Multiple trajectories can be provided, miic will consider that a new trajectory
#' starts each time a smaller time step than the one of the previous row is encountered.
#'
#' @param state_order [a data frame, optional, NULL by default]
#'
#' A data frame providing extra information for variables. It must have d rows
#' where d is the number of input variables and possible columns are described
#' below. For optional columns, if they are not provided or contain missing
#' values, default values suitable for \emph{input_data} will be used.
#'
#' \emph{"var_names"} (required) contains the name of each variable as specified
#' by colnames(input_data). In temporal mode, the time steps column should
#' not be mentioned in the variables list.
#'
#' \emph{"var_type"} (optional) contains a binary value that specifies if each
#' variable is to be considered as discrete (0) or continuous (1).
#'
#' \emph{"levels_increasing_order"} (optional) contains a single character string
#' with all of the unique levels of the ordinal variable in increasing order,
#' delimited by comma ','. It will be used during the post-processing to compute
#' the sign of an edge using Spearman's rank correlation. If a variable is
#' continuous or is categorical but not ordinal, this column should be NA.
#'
#' \emph{"is_contextual"} (optional) contains a binary value that specifies if a
#' variable is to be considered as a contextual variable (1) or not (0).
#' Contextual variables cannot be the child node of any other variable (cannot
#' have edge with arrowhead pointing to them).
#'
#' \emph{"is_consequence"} (optional) contains a binary value that specifies if a
#' variable is to be considered as a consequence variable (1) or not (0).
#' Edges between consequence variables are ignored, consequence variables
#' cannot be the parent node of any other variable and cannot be used as
#' contributors. Edges between a non consequence and consequence variables
#' are pre-oriented toward the consequence.
#'
#' Several other columns are possible in temporal mode:
#'
#' \emph{"n_layers"} (optional) contains an integer value that specifies the
#' number of layers to be considered for the variable.
#' Note that if a \emph{"n_layers"} column is present in the \emph{state_order},
#' its values will overwrite the function parameter.
#'
#' \emph{"delta_t"} (optional) contains an integer value that specifies the number
#' of time steps between each layer for the variable.
#' Note that if a \emph{"delta_t"} column is present in the \emph{state_order},
#' its values will overwrite the function parameter.
#'
#' \emph{"movavg"} (optional) contains an integer value that specifies the size of
#' the moving average window to be applied to the variable.
#' Note that if \emph{"movavg"} column is present in the \emph{state_order},
#' its values will overwrite the function parameter.
#'
#' @param true_edges [a data frame, optional, NULL by default]
#'
#' A data frame containing the edges of the true graph for
#' computing performance after the run.\cr
#' In standard mode, the expected layout is a two columns data frame, each row
#' representing a true edge with in each column, the variable names.
#' Variables names must exist in the \emph{input_data} data frame.\cr
#' In temporal mode, the expected layout is a three columns data frame,
#' with the first two columns being variable names and the third the lag.
#' Variables names must exist in the \emph{input_data} data frame and the lag
#' must be valid in the time unfolded graph. e.g. a row var1, var2, 3 is valid
#' with \emph{n_layers} = 4 + \emph{delta_t} = 1 or
#' \emph{n_layers} = 2 + \emph{delta_t} = 3
#' but not for \emph{n_layers} = 2 + \emph{delta_t} = 2 as there is no matching
#' edge in the time unfolded graph. Please note that the order is important:
#' var1, var2, 3 is interpreted as var1_lag3 - var2_lag0. Please note also that,
#' for contextual variables that are not lagged, the expected value in the
#' third column for the time lag is NA.
#'
#' @param black_box [a data frame, optional, NULL by default]
#'
#' A data frame containing pairs of variables that will be considered
#' as independent during the network reconstruction. In practice, these edges
#' will not be included in the skeleton initialization and cannot be part of
#' the final result.\cr
#' In standard mode, the expected layout is a two columns data frame, each row
#' representing a forbidden edge with in each column, the variable names.
#' Variables names must exist in the \emph{input_data} data frame.\cr
#' In temporal mode, the expected layout is a three columns data frame,
#' with the first two columns being variable names and the third the lag.
#' Variables names must exist in the \emph{input_data} data frame and the lag
#' must be valid in the time unfolded graph. e.g. a row var1, var2, 3 is valid
#' with \emph{n_layers} = 4 + \emph{delta_t} = 1 or
#' \emph{n_layers} = 2 + \emph{delta_t} = 3
#' but not for \emph{n_layers} = 2 + \emph{delta_t} = 2 as there is no matching
#' edge in the time unfolded graph. Please note that the order is important:
#' var1, var2, 3 is interpreted as var1_lag3 - var2_lag0. Please note also that,
#' for contextual variables that are not lagged, the expected value in the
#' third column for the time lag is NA.
#'
#' @param n_threads [a positive integer, optional, 1 by default]
#'
#' When set greater than 1, n_threads parallel threads will be used for computation. Make sure
#' your compiler is compatible with openmp if you wish to use multithreading.
#'
#' @param cplx [a string, optional, "nml" by default, possible values:
#' "nml", "mdl"]
#'
#' In practice, the finite size of the input  dataset requires that
#' the 2-point and 3-point information measures should be \emph{shifted}
#' by a \emph{complexity} term. The finite size corrections can be based on
#' the Minimal Description Length (MDL) criterion.
#' However, the MDL complexity criterion tends to underestimate the
#' relevance of edges connecting variables with many different categories,
#' leading to the removal of false negative edges. To avoid such biases
#' with finite datasets, the (universal) Normalized Maximum Likelihood (NML)
#' criterion can be used (see Affeldt \emph{et al.}, UAI 2015).
#'
#' @param orientation [a boolean value, optional, TRUE by default]
#'
#' The miic network skeleton can be partially directed by orienting
#' edge directions, based on the sign and magnitude of the conditional
#' 3-point information of unshielded triples and, in temporal mode, using time.
#' If set to FALSE, the orientation step is not performed.
#'
#' @param ori_proba_ratio [a floating point between 0 and 1, optional,
#' 1 by default]
#'
#' The threshold when deducing the type of an edge tip (head/tail)
#' from the probability of orientation.
#' For a given edge tip, denote by p the probability of it being a head,
#' the orientation is accepted if (1 - p) / p < \emph{ori_proba_ratio}.
#' 0 means reject all orientations, 1 means accept all orientations.
#'
#' @param ori_consensus_ratio [a floating point between 0 and 1, optional,
#' NULL by default]
#'
#' The threshold when deducing the type of an consensus edge tip (head/tail)
#' from the average probability of orientation.
#' For a given edge tip, denote by p the probability of it being a head,
#' the orientation is accepted if (1 - p) / p < \emph{ori_consensus_ratio}.
#' 0 means reject all orientations, 1 means accept all orientations.
#' If not supplied, the \emph{ori_consensus_ratio} will be initialized with
#' the \emph{ori_proba_ratio} value.
#'
#' @param propagation [a boolean value, optional, FALSE by default]
#'
#' If set to FALSE, the skeleton is partially oriented with only the
#' v-structure orientations. Otherwise, the v-structure orientations are
#' propagated to downstream undirected edges in unshielded triples following
#' the propagation procedure, relying on probabilities (for more details,
#' see Verny \emph{et al.}, PLoS Comp. Bio. 2017).
#'
#' @param latent [a string, optional, "orientation" by default, possible
#' values: "orientation", "no", "yes"]
#'
#' When set to "yes", the network reconstruction is taking into account hidden
#' (latent) variables. When set to "orientation", latent variables are not
#' considered during the skeleton reconstruction but allows bi-directed edges
#' during the orientation.
#' Dependence between two observed variables due to a latent variable is
#' indicated with a '6' in the adjacency matrix and in the network
#' edges.summary and by a bi-directed edge in the (partially) oriented graph.
#'
#' @param n_eff [a positive integer, optional, -1 by default]
#'
#' In standard mode, the n samples given in the \emph{input_data} data frame are
#' expected to be independent. In case of correlated samples such as in
#' Monte Carlo sampling approaches, the effective number of independent samples
#' \emph{n_eff} can be estimated using the decay of the autocorrelation function
#' (Verny et al., PLoS Comp. Bio. 2017). This effective number \emph{n_eff}
#' of independent samples can be provided using this parameter.
#'
#' @param n_shuffles [a positive integer, optional, 0 by default]
#'
#' The number of shufflings of the original dataset in order to evaluate
#' the edge specific confidence ratio of all inferred edges.
#' Default is 0: no confidence cut is applied. If the number of shufflings
#' is set to an integer > 0, the confidence threshold must also be > 0
#' (e.g. \emph{n_shuffles} = 100 and \emph{conf_threshold} = 0.01).
#'
#' @param conf_threshold [a positive floating point, optional, 0 by default]
#'
#' The threshold used to filter the less probable edges following the skeleton
#' step. See Verny \emph{et al.}, PLoS Comp. Bio. 2017.
#' Default is 0: no confidence cut is applied. If the the confidence threshold
#' is set > 0, the number of shufflings must also be > 0
#' (e.g. \emph{n_shuffles} = 100 and \emph{conf_threshold} = 0.01).
#'
#' @param sample_weights [a numeric vector, optional, NULL by default]
#'
#' An vector containing the weight of each observation.  If defined, it must be
#' a vector of floats in the range [0,1] of size equal to the number of samples.
#'
#' @param test_mar [a boolean value, optional, TRUE by default]
#'
#' If set to TRUE, distributions with missing values will be tested with
#' Kullback-Leibler divergence: conditioning variables for the given link
#' \eqn{X - Y}, \eqn{Z} will be considered only if the divergence
#' between the full distribution and the non-missing distribution
#' \eqn{KL(P(X,Y) | P(X,Y)_{!NA})} is low enough (with \eqn{P(X,Y)_{!NA}} as
#' the joint distribution of \eqn{X} and \eqn{Y} on samples which are
#' not missing on Z.
#' This is a way to ensure that data are missing at random for the considered
#' interaction and to avoid selection bias.
#'
#' @param consistent [a string, optional, "no" by default, possible values:
#' "no", "orientation", "skeleton"]
#'
#' If set to "orientation": iterate over skeleton and orientation steps to
#' ensure consistency of the network.
#' If set to "skeleton": iterate over skeleton step to get a consistent skeleton,
#' then orient edges and discard inconsistent orientations to ensure consistency
#' of the network (see Li \emph{et al.}, NeurIPS 2019 for details).
#'
#' @param max_iteration [a positive integer, optional, 100 by default]
#'
#' When the \emph{consistent} parameter is set to "skeleton" or "orientation",
#' the maximum number of iterations allowed when trying to find a consistent
#' graph.
#'
#' @param consensus_threshold [a floating point between 0.5 and 1.0, optional,
#' 0.8 by default]
#'
#' When the \emph{consistent} parameter is set to "skeleton" or "orientation"
#' and when the result graph is inconsistent or is a union of more than
#' one inconsistent graphs, a consensus graph will be produced based on
#' a pool of graphs.
#' If the result graph is inconsistent, then the pool is made of
#' \emph{max_iteration} graphs from the iterations, otherwise it is made of
#' those graphs in the union.
#' In the consensus graph, an edge is present when the proportion of non-zero
#' status in the pool is above the threshold. For example, if the pool contains
#' [A, B, B, 0, 0], where "A", "B" are different status of the edge and "0"
#' indicates the absence of the edge. Then the edge is set to connected ("1")
#' if the proportion of non-zero status (0.6 in the example) is equal to
#' or higher than \emph{consensus_threshold}. (When set to connected,
#' the orientation of the edge will be further determined by the average
#' probability of orientation.)
#'
#' @param negative_info [a boolean value, optional, FALSE by default]
#'
#' For test purpose only. If TRUE, negative shifted mutual information is
#' allowed during the computation when mutual information is inferior to the
#' complexity term.
#' For small dataset with complicated structures, e.g. discrete variables with
#' many levels, allowing for negative shifted mutual information may help
#' identifying weak v-structures related to those discrete variables,
#' as the negative three-point information in those cases will come from
#' the difference between two negative shifted mutual information terms
#' (expected to be negative due to the small sample size).
#' However, under this setting, a v-structure (X -> Z <- Y) in the final graph
#' does not necessarily imply that X is dependent on Y conditioning on Z,
#' As a consequence, the interpretability of the final graph
#' is hindered. In practice, it's advised to keep this parameter as FALSE.
#'
#' @param mode [a string, optional, "S" by default, possible values are
#' "S": Standard (IID samples) or "TS": Temporal Stationary"]
#'
#' When temporal mode is activated, the time information must be provided
#' in the first column of \emph{input_data}. For more details about temporal
#' stationary mode, see Simon \emph{et al.}, eLife reviewed preprint 2024.
#'
#' @param n_layers [an integer, optional, NULL by default, must be >= 2
#' if supplied]
#'
#' Used only in temporal mode, \emph{n_layers} defines the number of layers
#' that will be considered for the variables in the time unfolded graph.
#' The layers will be distant of \emph{delta_t} time steps.
#' If not supplied, the number of layers is estimated from the dynamic of the
#' dataset and the maximum number of nodes \emph{max_nodes} allowed in the
#' final lagged graph.
#'
#' @param delta_t [an integer, optional, NULL by default, must be >= 1
#' if supplied]
#'
#' Used only in temporal mode, \emph{delta_t} defines the number of time steps
#' between each layer.
#' i.e. on 1000 time steps with \emph{n_layers} = 3 and \emph{delta_t} = 7,
#' the time steps kept for the samples conversion will be 1, 8, 15
#' for the first sample, the next sample will use 2, 9, 16 and so on.
#' If not supplied, the number of time steps between layers is estimated
#' from the dynamic of the dataset and the number of layers.
#'
#' @param movavg [an integer, optional, NULL by default, must be >= 2
#' if supplied]
#'
#' Used only in temporal mode. When supplied, a moving average operation is
#' applied to all integer and numeric variables that are not contextual
#' variables.
#'
#' @param keep_max_data [a boolean value, optional, FALSE by default]
#'
#' Used only in temporal mode. If TRUE, rows where some NAs have been
#' introduced during the moving averages and lagging will be kept
#' whilst they will be dropped if FALSE.
#'
#' @param max_nodes [an integer, optional, 50 by default]
#'
#' Used only in temporal mode and if the \emph{n_layers} or \emph{delta_t}
#' parameters are not supplied. \emph{max_nodes} is used as the maximum number
#' of nodes in the final graph to compute \emph{n_layers} and/or \emph{delta_t}.
#' The default is 50 to produce quick runs and can be increased up to 200
#' or 300 on recent computers to produce more precise results.
#'
#' @param verbose [a boolean value, optional, FALSE by default]
#'
#' If TRUE, debugging output is printed.
#'
#' @return A \emph{miic-like} object that contains:
#' \itemize{
#'  \item{\emph{all.edges.summary:} a data frame with information about the relationship between
#'  each pair of variables
#'  \itemize{
#'  \item{ \emph{x:} X node}
#'  \item{ \emph{y:} Y node}
#'  \item{ \emph{type:} contains 'N' if the edge has been removed or 'P' for
#'  retained edges. If the true graph is supplied in the \emph{true_edges}
#'  parameter, 'P' becomes 'TP' (True Positive) or 'FP' (False Positive),
#'  while 'N' becomes 'TN' (True Negative) or 'FN' (False Negative).
#'  Note that, as the \emph{all.edges.summary} does not contain all the
#'  negative edges, edges not present are 'TN'.}
#'  \item{ \emph{ai:} the contributing nodes found by the method which
#'  participate in the mutual information between \emph{x} and \emph{y},
#'  and possibly separate them.}
#'  \item{ \emph{raw_contributions:} describes the share of total mutual
#'  information between \emph{x} and \emph{y} explained by each contributor.}
#'  \item{ \emph{contributions:} describes the share of remaining mutual
#'  information between \emph{x} and \emph{y} explained  by each successive
#'  contributors.}
#'  \item{ \emph{info:} provides the pairwise mutual information times
#'  \emph{Nxyi} for the pair (\emph{x}, \emph{y}).}
#'  \item{ \emph{info_cond:} provides the conditional mutual information
#'  times \emph{Nxy_ai} for the pair (\emph{x}, \emph{y}) when conditioned
#'  on the collected nodes \emph{ai}. It is
#'  equal to the \emph{info} column when \emph{ai} is an empty set.}
#'  \item{ \emph{cplx:} gives the computed complexity between the (\emph{x},
#'  \emph{y}) variables taking into account the contributing nodes \emph{ai}.
#'  Edges that have have more conditional information \emph{info_cond}
#'  than \emph{cplx} are retained in the final graph.}
#'  \item{ \emph{Nxy_ai:} gives the number of complete samples on which the
#'  information and the  complexity have been computed. If the input dataset
#'  has no missing value, the number of samples is the same for all pairs
#'  and corresponds to the total number of samples.}
#'  \item{ \emph{info_shifted:} represents the \emph{info} - \emph{cplx} value.
#'  It is a way to quantify the strength of the edge (\emph{x}, \emph{y}).}
#'  \item{ \emph{infOrt:} the orientation of the edge (\emph{x}, \emph{y}).
#'  It is the same value as in the adjacency matrix at row \emph{x} and
#'  column \emph{y} : 1 for unoriented, 2 for an edge from X to Y,
#'  -2 from Y to X and 6 for bidirectional.}
#'  \item{ \emph{trueOrt:} the orientation of the edge (\emph{x}, \emph{y})
#'  present in the true edges are provided.}
#'  \item{ \emph{isOrtOk:} information about the consistency of the inferred
#'  graph’s orientations with a reference graph is given (if true edges
#'  are provided).
#'  'Y': the orientation is consistent; 'N': the orientation is not consistent
#'  with the PAG (Partial Ancestor Graph) derived from the given true graph.}
#'  \item{ \emph{sign:} the sign of the partial correlation between variables
#'  \emph{x} and \emph{y}, conditioned on the contributing nodes \emph{ai}.}
#'  \item{ \emph{partial_correlation:} value of the partial correlation for the
#'  edge (\emph{x}, \emph{y}) conditioned on the contributing nodes \emph{ai}.}
#'  \item{ \emph{is_causal:} details about the nature of the arrow tip for a
#'  directed edge. A directed edge in a causal graph does not necessarily imply
#'  causation but it does imply that the cause-effect relationship is not the
#'  other way around. An arrow-tip which is itself downstream of another
#'  directed edge suggests stronger causal sense and is marked by a 'Y',
#'  or 'N' otherwise.}
#'  \item{ \emph{proba:} probabilities for the inferred orientation, derived
#'  from the three-point mutual information (cf Affeldt & Isambert, UAI 2015
#'  proceedings) and noted as p(x->y);p(x<-y).}
#'  \item{ \emph{confidence:} this column is computed when the confidence cut
#'  is activated. It represents the ratio between the probability to reject
#'  the edge (\emph{x}, \emph{y}) in the dataset versus the mean probability
#'  to do the same in multiple (user defined) number of randomized datasets.}
#'  }
#'  }
#'
#'  \item{\emph{orientations.prob:} this data frame lists the orientation
#'  probabilities of the two edges of all unshielded triples of the
#'  reconstructed network with the structure: node1 -- mid-node -- node2:
#'  \itemize{
#'  \item{ \emph{node1:} node at the end of the unshielded triplet}
#'  \item{ \emph{p1:} probability of the arrowhead node1 <- mid-node}
#'  \item{ \emph{p2:} probability of the arrowhead node1 -> mid-node}
#'  \item{ \emph{mid-node:} node at the center of the unshielded triplet}
#'  \item{ \emph{p3:} probability of the arrowhead mid-node <- node2}
#'  \item{ \emph{p4:} probability of the arrowhead mid-node -> node2}
#'  \item{ \emph{node2:} node at the end of the unshielded triplet}
#'  \item{ \emph{NI3:} 3 point (conditional) mutual information * N}
#'  \item{ \emph{Conflict:} indicates if there a conflict between the
#'  computed probabilities and the NI3 value}
#'  }
#'  }
#'
#'  \item {\emph{adj_matrix:} the adjacency matrix is a square matrix used to
#'  represent the inferred graph. The entries of the matrix indicate whether
#'  pairs of vertices are adjacent or not in the graph. The matrix can be read
#'  as a (row, column) set of couples where the row represents the source node
#'  and the column the target node. Since miic can reconstruct mixed networks
#'  (including directed, undirected and bidirected edges), we will have a
#'  different digit for each case:
#'  \itemize{
#'  \item{ 1: (\emph{x}, \emph{y}) edge is undirected}
#'  \item{ 2: (\emph{x}, \emph{y}) edge is directed as \emph{x} -> \emph{y} }
#'  \item{ -2: (\emph{x}, \emph{y}) edge is directed as \emph{x} <- \emph{y} }
#'  \item{ 6: (\emph{x}, \emph{y}) edge is bidirected}
#'  }
#'  }
#'
#'  \item {\emph{proba_adj_matrix:} the probability adjacency matrix is
#'  a square  matrix used to represent the orientation probabilities associated
#'  to the edges tips. The value at ("row", "column") is the probability,
#'  for the edge between "row" and "column" nodes, of the edge tip on the "row"
#'  side.  A probability less than 0.5 is an indication of a possible tail
#'  (cause) and a probability greater than 0.5 a possible head (effect).
#'  }
#'
#'  \item {\emph{params:} the list of parameters used for the network
#'  reconstruction. The parameters not supplied are initialized to their default
#'  values. Otherwise, the parameters are checked and corrected if necessary.}
#'
#'  \item {\emph{state_order:} the state order used for the network
#'  reconstruction. If no state order is supplied, it is generated by using
#'  default values. Otherwise, it is the state order checked and corrected
#'  if necessary.}
#'
#'  \item {\emph{black_box:} present only if a black box has been supplied,
#' the black box, checked and corrected if necessary, used for the network
#' reconstruction.}
#'
#'  \item {\emph{true_edges:} present only if the true edges have been supplied,
#' the true edges, checked and corrected if necessary, used for the network
#' evaluation.}
#'
#'  \item {\emph{tmiic:} present only in temporal mode.
#'  Named list containing the full list of edges completed by stationarity,
#'  the lagged state order and, if a black box or true edges have been supplied,
#'  the lagged versions of these inputs.}
#' }
#'
#' @export
#' @useDynLib miic
#' @import Rcpp
#'
#' @examples
#' library(miic)
#'
#' # EXAMPLE HEMATOPOIESIS
#' data(hematoData)
#'
#' # execute MIIC (reconstruct graph)
#' miic.res <- miic(
#'   input_data = hematoData[1:1000,], latent = "yes",
#'   n_shuffles = 10, conf_threshold = 0.001
#' )
#'
#' # plot graph
#' if(require(igraph)) {
#'  plot(miic.res, method="igraph")
#' }
#'
#' \donttest{
#' # write graph to graphml format. Note that to correctly visualize
#' # the network we created the miic style for Cytoscape (http://www.cytoscape.org/).
#'
#' miic.write.network.cytoscape(g = miic.res, file = file.path(tempdir(), "temp"))
#'
#' # EXAMPLE CANCER
#' data(cosmicCancer)
#' data(cosmicCancer_stateOrder)
#' # execute MIIC (reconstruct graph)
#' miic.res <- miic(
#'   input_data = cosmicCancer, state_order = cosmicCancer_stateOrder, latent = "yes",
#'   n_shuffles = 100, conf_threshold = 0.001
#' )
#'
#' # plot graph
#' if(require(igraph)) {
#'  plot(miic.res)
#' }
#'
#' # write graph to graphml format. Note that to correctly visualize
#' # the network we created the miic style for Cytoscape (http://www.cytoscape.org/).
#' miic.write.network.cytoscape(g = miic.res, file = file.path(tempdir(), "temp"))
#'
#' # EXAMPLE OHNOLOGS
#' data(ohno)
#' data(ohno_stateOrder)
#' # execute MIIC (reconstruct graph)
#' miic.res <- miic(
#'   input_data = ohno, latent = "yes", state_order = ohno_stateOrder,
#'   n_shuffles = 100, conf_threshold = 0.001
#' )
#'
#' # plot graph
#' if(require(igraph)) {
#'  plot(miic.res)
#' }
#'
#' # write graph to graphml format. Note that to correctly visualize
#' # the network we created the miic style for Cytoscape (http://www.cytoscape.org/).
#' miic.write.network.cytoscape(g = miic.res, file = file.path(tempdir(), "temp"))
#'
#' # EXAMPLE COVID CASES (time series demo)
#' data(covidCases)
#' # execute MIIC (reconstruct graph in temporal mode)
#' tmiic.res <- miic(input_data = covidCases, mode = "TS", n_layers = 3, delta_t = 1, movavg = 14)
#'
#' # to plot the default graph (compact)
#' if(require(igraph)) {
#'  plot(tmiic.res)
#' }
#'
#' # to plot the raw temporal network Using igraph
#' if(require(igraph)) {
#'   plot(tmiic.res, display="raw")
#' }
#'
#' # to plot the full temporal network Using igraph
#' if(require(igraph)) {
#'   plot(tmiic.res, display="lagged")
#' }
#'
#' }
#'
miic <- function(input_data,
                 state_order = NULL,
                 true_edges = NULL,
                 black_box = NULL,
                 n_threads = 1,
                 cplx = "nml",
                 orientation = TRUE,
                 ori_proba_ratio = 1,
                 ori_consensus_ratio = NULL,
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
                 mode = "S",
                 n_layers = NULL,
                 delta_t = NULL,
                 movavg = NULL,
                 keep_max_data = FALSE,
                 max_nodes = 50,
                 verbose = FALSE)
  {
  if (verbose)
    miic_msg ("Start MIIC...")
  if ( is.null(mode) || ( ! (mode %in% MIIC_VALID_MODES) ) )
    miic_error ("parameters check", "invalid mode ", mode,
      ". Possible modes are S (Standard), TS (Temporal Stationnary).")
  if (mode %in% MIIC_TEMPORAL_MODES)
      miic_msg ("Using temporal mode of MIIC")
  #
  # Check base inputs
  #
  input_data = check_input_data (input_data, mode)
  params = check_parameters (input_data = input_data,
                              n_threads = n_threads,
                              cplx = cplx,
                              orientation = orientation,
                              ori_proba_ratio = ori_proba_ratio,
                              ori_consensus_ratio = ori_consensus_ratio,
                              propagation = propagation,
                              latent = latent,
                              n_eff = n_eff,
                              n_shuffles = n_shuffles,
                              conf_threshold = conf_threshold,
                              sample_weights = sample_weights,
                              test_mar = test_mar,
                              consistent = consistent,
                              max_iteration = max_iteration,
                              consensus_threshold = consensus_threshold,
                              mode = mode,
                              negative_info = negative_info,
                              verbose = verbose)
  state_order = check_state_order (input_data, state_order, params$mode)
  black_box = check_other_df (input_data, state_order,
                              black_box, "black box", params$mode)
  true_edges = check_other_df (input_data, state_order,
                               true_edges, "true edges", params$mode)
  #
  # Extra steps depending on the mode
  #
  if (! (mode %in% MIIC_TEMPORAL_MODES) )
    non_lagged_state_order = NULL
  else
    {
    # Check temporal parameters and state_order
    #
    state_order = tmiic_check_state_order_part1 (state_order)
    list_ret = tmiic_check_parameters (state_order = state_order,
                                       params = params,
                                       n_layers = n_layers,
                                       delta_t = delta_t,
                                       movavg = movavg,
                                       keep_max_data = keep_max_data,
                                       max_nodes = max_nodes)
    params = list_ret$params
    state_order = tmiic_check_state_order_part2 (list_ret$state_order)
    list_ts = tmiic_extract_trajectories (input_data)
    list_ts = tmiic_movavg (list_ts, state_order$movavg,
                            keep_max_data=params$keep_max_data,
                            verbose_level=ifelse (params$verbose, 2, 1) )
    state_order = tmiic_estimate_dynamic (list_ts, state_order,
                            max_nodes=params$max_nodes,
                            verbose_level=ifelse (params$verbose, 2, 1) )
    #
    # Lag data and other inputs accordingly
    #
    non_lagged_state_order = state_order
    non_lagged_true_edges = true_edges
    non_lagged_black_box = black_box
    state_order = tmiic_lag_state_order (non_lagged_state_order)
    true_edges = tmiic_lag_other_df (non_lagged_state_order, true_edges)
    true_edges = tmiic_check_other_df_after_lagging (state_order$var_names,
                                                     true_edges, "true edges")
    black_box = tmiic_lag_other_df (non_lagged_state_order, black_box)
    black_box = tmiic_check_other_df_after_lagging (state_order$var_names,
                                                     black_box, "black box")
    list_ts = tmiic_lag_input_data (list_ts, state_order,
                                    keep_max_data=params$keep_max_data)
    input_data = tmiic_group_trajectories (list_ts)
    #
    # Check number of unique values per variable and review discrete/continuous
    # after lagging as some columns may have less number of unique values
    #
    state_order = tmiic_check_after_lagging (input_data, state_order)
    #
    # Adjust n_eff if delta_t > 1 and no eff supplied by the user
    #
    avg_delta_t = mean (state_order$delta_t[state_order$is_contextual == 0])
    if ( (avg_delta_t > 1) && (params$n_eff == -1) )
      {
      params$n_eff = trunc (nrow (input_data) / avg_delta_t)
      miic_msg ("Note : the n_eff has been set to ", params$n_eff,
                " (nb lagged samples= ", nrow (input_data),
                " / delta_t=", round(avg_delta_t, 2), ").")
      }
    }
  #
  # Convert discrete vars as factors
  #
  for ( i in 1:nrow(state_order) )
    if (state_order[i, "var_type"] == 0)
      input_data[, i] <- factor (input_data[, i])
  #
  # Call C++ reconstruction
  #
  if (verbose)
    miic_msg ("-> Start reconstruction...")
  res <- miic.reconstruct (input_data = input_data,
                           n_threads = params$n_threads,
                           cplx = params$cplx,
                           latent = params$latent,
                           n_eff = params$n_eff,
                           black_box = black_box,
                           n_shuffles = params$n_shuffles,
                           orientation = params$orientation,
                           ori_proba_ratio = params$ori_proba_ratio,
                           propagation = params$propagation,
                           conf_threshold = params$conf_threshold,
                           verbose = params$verbose,
                           is_contextual = state_order$is_contextual,
                           is_consequence = state_order$is_consequence,
                           is_continuous = state_order$var_type,
                           sample_weights = params$sample_weights,
                           test_mar = params$test_mar,
                           consistent = params$consistent,
                           mode = params$mode,
                           n_layers = non_lagged_state_order$n_layers,
                           delta_t = non_lagged_state_order$delta_t,
                           max_iteration = params$max_iteration,
                           negative_info = params$negative_info)
  if (res$interrupted)
    stop("Interupted by user")
  if (verbose)
    miic_msg ("-> End reconstruction...")
  #
  # Post-traitment
  #
  res$all.edges.summary <- summarizeResults (
    observations = input_data,
    results = res,
    true_edges = true_edges,
    state_order = state_order,
    consensus_threshold = params$consensus_threshold,
    ori_consensus_ratio = params$ori_consensus_ratio,
    latent = (params$latent != "no"),
    propagation = params$propagation,
    verbose = params$verbose)

  res$params = params
  if (! (mode %in% MIIC_TEMPORAL_MODES) )
    {
    class(res) <- "miic"
    res$state_order = state_order
    res$black_box = black_box
    res$true_edges = true_edges
    }
  else
    {
    class(res) <- "tmiic"
    #
    # clean state_order structure to remove extra columns used internally
    #
    non_lagged_state_order = non_lagged_state_order[,
      colnames(non_lagged_state_order) %in% STATE_ORDER_TEMPORAL_VALID_COLUMNS]
    res$state_order = non_lagged_state_order
    res$black_box = non_lagged_black_box
    res$true_edges = non_lagged_true_edges

    state_order = state_order[,
        colnames(state_order) %in% STATE_ORDER_TEMPORAL_VALID_COLUMNS]
    #
    # The output of the reconstruction is the "raw" temporal graph, without
    # edges identical by stationarity. To have the "real" temporal graph,
    # we duplicate the edges using the stationary assumption and this "real"
    # graph is stored the "all.edges.stationarity" data frame.
    #
    edges_dup_stat = tmiic_repeat_edges_over_history (res)
    res$tmiic <- list (lagged_state_order = state_order,
                       lagged_black_box = black_box,
                       lagged_true_edges = true_edges,
                       all.edges.stationarity = edges_dup_stat)
    }
  return(res)
  }
