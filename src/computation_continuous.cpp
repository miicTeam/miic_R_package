
#include "computation_continuous.h"

#include <algorithm>  // std::min, std::transform
#include <limits>
#include <numeric>  // std::accumulate
#include <cmath>    // std::pow
#include <Rcpp.h>   // Rcpp::Rcout

#include "linear_allocator.h"
#include "mutual_information.h"
#include "structure.h"

constexpr double kEpsilon = 1.e-5;
constexpr double kPrecision = 1.e-10;
constexpr int kMaxIter = 50;
constexpr int kMaxIterOnU = 3;

namespace miic {
namespace computation {

using std::accumulate;
using std::copy;
using std::min;
using std::max;
using namespace miic::structure;
using namespace miic::utility;

namespace {

//------------------------------------------------------------------------------
// resetCutPoints
//------------------------------------------------------------------------------
// From the number of unique values of each variable and the number of samples,
// determine for each continuous variable the length of the bin and the numbers
// of bin. Then fix a cut every "length of bin"
//
// Inputs:
// - levels: Vector, number of unique values for the variables
// - is_continuous: Vector, flags (1=continuous, 0=discrete) for the variables
// - vars_to_process: Vector of the variables indexes to be processed
// - vars_to_process_begin: int, first index in vars_to_process to consider
// - vars_to_process_end: int, last index in vars_to_process to consider
// - init_nbin: int, number of bins by default
// - n_samples: int, number of samples
// Outputs:
// - cuts: Grid (nb variables to process x nb max cuts), the cuts
//------------------------------------------------------------------------------
void resetCutPoints (const TempVector<int>& levels,
                     const TempVector<int>& is_continuous,
                     const TempVector<int>& vars_to_process,
                     int vars_to_process_begin,
                     int vars_to_process_end,
                     int init_nbin, int n_samples,
                     TempGrid2d<int>& cuts) {
  for (int vars_to_process_idx = vars_to_process_begin;
       vars_to_process_idx < vars_to_process_end;
       ++vars_to_process_idx) {
    if (is_continuous[vars_to_process[vars_to_process_idx]] != 1)
      continue;
    //
    // Determine the length of a bin and the number of bins
    //
    int n_bins = min (init_nbin, levels[vars_to_process[vars_to_process_idx]]);
    int length_bin = n_samples / n_bins;
    if (length_bin < 1) {
      length_bin = 1;
      n_bins = n_samples;
    }
    //
    // Reset cuts: fix a cut every length_bin until n_bins, then add the
    // conventional last cut (= nb samples) and set the remaining positions to 0
    //
    for (int bin_idx = 0; bin_idx < n_bins - 1; ++bin_idx)
      cuts(vars_to_process_idx, bin_idx) = bin_idx * length_bin + length_bin - 1;
    cuts(vars_to_process_idx, n_bins - 1) = n_samples - 1;
    for (size_t bin_idx = n_bins; bin_idx < cuts.n_cols(); ++bin_idx)
      cuts(vars_to_process_idx, bin_idx) = 0;
  }
}

//------------------------------------------------------------------------------
// updateFactors
//------------------------------------------------------------------------------
// Transform data into factors using the cuts supplied in parameters.
// i.e.: for the variable values [0,2,8,4,5,3] with cuts positions at [2,5],
// the computed factors will be [0,0,1,1,1,0].
// Note that when some values are repeated, the updateFactors will not strictly
// respect the cuts positions to ensure that a value is not be associated with
// different factors.
// i.e.: [0,2,8,3,5,3] with cuts positions at [2,5] will not be transformed
// into [0,0,1,0,1,1] as the value 3 would be transformed into 0 and 1.
// updateFactors will not apply the cut as long as it encounters the same value.
// So in this example, the factors will be [0,0,1,0,1,0].
//
// Inputs:
// - data: Grid (variables x samples), the input data (ranked)
// - data_ordered: Grid (variables x samples), the indexes in data of ordered values
// - cuts: Grid (nb variables to process x nb max cuts), the cuts to applied-
// - is_continuous: Vector, flags (1=continuous, 0=discrete) for each variable
// - vars_to_process: Vector of the variables indexes to be processed
// - vars_to_process_begin: int, first index in vars_to_process to consider
// - vars_to_process_end: int, last index in vars_to_process to consider
// Outputs:
// - factors: Grid (nb variables to process x samples), the data transformed into factors
// - cuts_levels: Vector, number of factors (unique values) in each cut
//------------------------------------------------------------------------------
void updateFactors (const TempGrid2d<int>& data,
                    const TempGrid2d<int>& data_ordered,
                    const TempGrid2d<int>& cuts,
                    const TempVector<int>& is_continuous,
                    const TempVector<int>& vars_to_process,
                    int vars_to_process_begin,
                    int vars_to_process_end,
                    TempGrid2d<int>& factors,
                    TempVector<int>& cuts_levels) {
  int n_samples = data_ordered.n_cols();
  for (int vars_to_process_idx = vars_to_process_begin;
       vars_to_process_idx < vars_to_process_end;
       ++vars_to_process_idx) {
    int data_var_idx = vars_to_process[vars_to_process_idx];
    if (is_continuous[data_var_idx] != 1)
      continue;
    //
    // For continuous variables, iterate over the samples values in an ordered way.
    // Affect all the values till the cut position to the value of 'level', then
    // increase the 'level' and proceed in the same way to the next cut position
    //
    int level = 0;
    for (int data_ordered_sample_idx = 0;
         data_ordered_sample_idx < n_samples;
         ++data_ordered_sample_idx) {
      int data_sample_idx = data_ordered (data_var_idx, data_ordered_sample_idx);
      if (data_ordered_sample_idx > cuts(vars_to_process_idx, level)) {
        //
        // When we reach a cut position, look if the sample value has changed,
        // if not, affect all the equal values to the current level
        // When a different value is encountered, increase the level
        //
        int data_sample_prec_idx = data_ordered (data_var_idx, data_ordered_sample_idx-1);
        int value_prec = data (data_var_idx, data_sample_prec_idx);
        int value_iter = data (data_var_idx, data_sample_idx);
        while (value_prec == value_iter) {
          factors (vars_to_process_idx, data_sample_idx) = level;
          ++data_ordered_sample_idx;
          if (data_ordered_sample_idx >= n_samples)
            break;
          data_sample_idx = data_ordered (data_var_idx, data_ordered_sample_idx);
          value_iter = data (data_var_idx, data_sample_idx);
        }
        if (data_ordered_sample_idx >= n_samples)
          break;
        ++level;
      }
      factors (vars_to_process_idx, data_sample_idx) = level;
    }
    cuts_levels[vars_to_process_idx] = level + 1;
  }
}

//------------------------------------------------------------------------------
// reconstructCutCoarse
//------------------------------------------------------------------------------
// Use the results of optimizeCutPoints to establish the cut points list:
// examines the cuts_idx vector starting from its end and follows the chain of
// the best cuts. For each cut, store the data row index coming from cuts_pos.
//
// inputs:
// - cuts_idx: Vector, recorded best cuts, possible values 0...size(cuts_idx).
//   Forms a chain with the best cuts indexes starting from the last position:
//   * only 0s -> stops (one bin)
//   * -k encountered -> stops (two bins ([0 k-1][k ..])
//   * k  encountered -> cut found [.. k-1][k ...]) and continue
// - cuts_pos: Vector, data row index for each cut
// - n_samples: int, number of samples
// output:
// - cuts_ret: Vector, cut points positions (row indexes) returned:
//   [0, cut_pos_a, cut_pos_b, ..., nb samples-1]
//------------------------------------------------------------------------------
template <typename Ccut, typename = IsIntContainer<Ccut>>
void reconstructCutCoarse (const TempVector<int>& cuts_idx,
                           const TempVector<int>& cuts_pos,
                           int n_samples,
                           Ccut& cuts_ret) {
  if (cuts_idx.back() == 0) {
    //
    // only 0s -> stops (one bin)
    //
    cuts_ret[0] = n_samples - 1;
    return;
  }
  //
  // Determine the number of cuts by following the chain of cuts
  // starting from the end
  //
  int ncuts = 0;
  int pop_cut = cuts_idx.back();
  while (pop_cut > 0) {
    ++ncuts;
    pop_cut = cuts_idx[pop_cut - 1];
  }
  if (pop_cut < 0)
    ncuts++;
  //
  // Construct the optimum cuts list by following the chain of cuts
  // starting from the end
  //
  cuts_ret[ncuts] = n_samples - 1;  // conventional last cut
  int cut_idx = ncuts - 1;
  cuts_ret[cut_idx] = cuts_pos.back();  // truly last cut
  --cut_idx;
  pop_cut = cuts_idx.back();
  while (pop_cut > 0 && cut_idx >= 0) {
    cuts_ret[cut_idx] = cuts_pos[pop_cut - 1];
    pop_cut = cuts_idx[pop_cut - 1];
    --cut_idx;
  }
}

//------------------------------------------------------------------------------
// optimizeCutPoints
//------------------------------------------------------------------------------
// Given factors of the target and joint variables, find a set of cut points for
// the target variable by maximizing mutual information over the possible cuts:
//
// Principle:
// coarse <- max (nb_vals_X // maxbin (maxbin=min(50,nb_samples)), 1)
// nb cuts max <-  nb_vals_X // coarse
// for each large bin in bins [ {0}, {0+1}, {0+1+2}, ... , {0+1+2+…+nb cuts max} ]
//   compute info of the large bin
//	 for each possible cut of the large bin
//	 	 compute info of the large_bin with the cut:
//	 	 {0,...,cut} {cut,...,end of large bin}
//		 if better info is found, memorize the cut
// unstack the cuts found starting from the end to get optimal cuts list
//
// Inputs:
// - data_ranked_target: Grid (n samples x 1 col), ranks of the target variable
//   (coming directly from a ranking on the input dataset, not discretized)
// - data_ordered_target: Grid (n samples x 1 col), indexes of ordered values
//   for the target variable. data_ordered_target[i] is the [i]th smallest value
//   i.e.: [0,0,3,4,1] is represented as [0,1,4,2,3] (representing the following
//   matches: 0->0, 1->0, 4->1, 2->3 and 3->4, giving the ordered values: 0,0,1,3,4).
// - data_factors_joint: Vector (n samples), discretized values of the joint
//   variable
// - data_factors_target: Vector (n samples), discretized values of the target
//   variable starting with one single bin
// - nb_factors_joint: int, number of factors (= unique values = levels) of the
//   discretized joint variable
// - nb_factors_target: int, number of factors (= unique values = levels) of the
//   discretized target variable
// - nb_factors_joint_for_cplx: int, number of factors (= unique values = levels)
//   of the discretized joint variable used for complexity. It is the product
//   of number of unique observed levels of all variables except the variable
//   being optimized.
// - nb_factors_target_prev: int, number of factors (= unique values = levels)
//   of the discretized target variable at the previous iteration.
//   Used for the complexity term
// - nb_uniq_vals_target: int, number of unique values in the non discretized
//   target variable (number of unique values after filtering on NAs)
// - weights: Vector (n samples), unique weights for each sample to account for
//   effective N != N, only used if the flag_weights is set to true.
// - flag_weights: bool, set to true if we use weights, false otherwise
// - maxbins: int, maximum number of bins allowed. It controls the resolution
//   of possible cut points. For example, 50 maximum cut points on a 1000
//   samples vector results in possible cut points every 1000/50 = 20 positions.
//   This speeds up significantly the dynamic programming at the cost of finding
//   approximated solutions.
// - cache: shared pointer, cache area storing previous calculations
// - cplx: int, is the choice of complexity : 0 for simple BIC (product of
//   all observed values) and 1 for NML with stochastic complexity and a
//   combinatorial term with the previous number of levels.
//
// Inputs / outputs:
// - cuts_target: Vector, contains the cut points for discretezing the
//   target variables. Modified at the end with the optimal solution.
//------------------------------------------------------------------------------
template <typename Cf0, typename Cf1, typename Ccut, typename =
    void_t<IsIntContainer<Cf0>, IsIntContainer<Cf1>, IsIntContainer<Ccut>>>
void optimizeCutPoints (const TempGrid2d<int>::ConstRow& data_ranked_target,
                        const TempGrid2d<int>::ConstRow& data_ordered_target,
                        const Cf0& data_factors_joint,
                        const Cf1& data_factors_target,
                        int nb_factors_joint,
                        int nb_factors_target,
                        int nb_factors_joint_for_cplx,
                        int nb_factors_target_prev,
                        int nb_uniq_vals_target,
                        const TempVector<double>& weights,
                        bool flag_weights,
                        int maxbins,
                        std::shared_ptr<CtermCache> cache,
                        int cplx,
                        Ccut&& cuts_target) {
  TempAllocatorScope scope;
  // Coarse graining, minimum length of a bin
  int coarse = std::max (ceil (1.0 * nb_uniq_vals_target / maxbins), 1.0);
  // Maximal number of possible cut points on the target
  int n_cuts_max = ceil (1.0 * nb_uniq_vals_target / coarse);
  int n_samples = data_ranked_target.size();
  //
  // For each bin, places "coarse" unique values (unique values from data_ranked)
  // in each bin and counts the number of occurrences of each factor.
  // These counts for each unitary bin, will be used later to compute the
  // entropy on combinations of bins
  // i.e.:
  // data_ranked_target    [0, 0, 0, 0, 0, 1, 2, 3, 0, 1, 2, 0]
  // data_ordered_target   [0, 1, 2, 3, 4, 8,11, 5, 9, 6,10, 7]
  // target discretized    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0] (first optimization, 1 bin)
  // assuming other var    [0, 1, 2, 3, 0, 1, 2, 3, 0, 1, 2, 3] discretized on 2 bin
  // discretized other var [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1]
  // data_factors_joint    [0, 0, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1]
  // now using coarse of 1 and max_bins = 4, we will iterate over sample
  // in an ordered way and the coarse_counts_joint will be:
  //   columns 0 1 (the factors, 0 or 1 possible from data_factors_joint)
  // rows   0: 4 3 (4 factors target=0 and joint=0, 3 factors target=0 and joint=1)
  // (bins) 1: 2 0 (new bin as target ranked value goes from 0 to 1 and coarse=1,
  //                then 2 factors target=0 and joint=0, 0 factors target=0 and joint=1)
  //        2: 0 2 (new bin as target ranked value goes from 1 to 2 and coarse=1,
  //                then 0 factors target=0 and joint=0, 2 factors target=0 and joint=1)
  //        3: 0 1 (new bin as target ranked value goes from 2 to 3 and coarse=1,
  //                then 0 factors target=0 and joint=0, 1 factors target=0 and joint=1)
  //
  // The sample_idx_breaks_per_cut will contains: [7, 9, 11, 12]
  // The cum_weights_per_bin, if "flag_weights" is true would contain:
  // [weights sum for samples [0->6], weights sum for samples[0->8],
  //  weights sum for samples [0->10], weights sum all samples[0->11] ]
  //
  TempVector<int> sample_idx_breaks_per_cut (n_cuts_max);
  TempVector<double> cum_weights_per_bin (n_cuts_max, 0);
  TempGrid2d<int> coarse_counts_target (n_cuts_max, nb_factors_target, 0);
  TempGrid2d<int> coarse_counts_joint (n_cuts_max, nb_factors_joint, 0);

  for (int bin_idx = 0, sample_idx = 0; bin_idx < n_cuts_max; ++bin_idx) {
    if (flag_weights && bin_idx > 0)
      cum_weights_per_bin[bin_idx] = cum_weights_per_bin[bin_idx - 1];

    int count_vals_uniq = 0;  // iterator on not repeated values
    while (count_vals_uniq < coarse && sample_idx < n_samples) {
      int data_factor_idx = data_ordered_target[sample_idx];
      int factor_target = data_factors_target[data_factor_idx];
      int factor_joint = data_factors_joint[data_factor_idx];
      ++coarse_counts_target (bin_idx, factor_target);
      ++coarse_counts_joint (bin_idx, factor_joint);

      if (flag_weights)
        cum_weights_per_bin[bin_idx] += weights[sample_idx];
      if (   (sample_idx + 1 < n_samples)
          && (data_ranked_target[data_ordered_target[sample_idx]] !=
              data_ranked_target[data_ordered_target[sample_idx + 1]]) )
        ++count_vals_uniq;
      ++sample_idx;
    }
    sample_idx_breaks_per_cut[bin_idx] = sample_idx;
  }
  //
  // To find the optimum cuts, we will iterate over an "accumulation" of bins
  // starting from the 1st bin, then {1st+2nd} bin, then {1st+2nd+3rd} bin,
  // ... until we consider all the bins.
  // Inside the loop, starting from the {1st+2nd} iteration, the function will
  // test if the cut between the 1st | 2nd bins is valuable then,
  // on the next iteration {1st+2nd+3rd} bin, cuts between 1st | {2nd+3rd} bins
  // and {1st=2nd} | 3rd bins will be tested and so on.
  //
  // Keep the result of each iteration
  TempVector<double> I (n_cuts_max);
  TempVector<double> Ik (n_cuts_max);
  // Chain of indexes of the cuts (1..n_cuts_max) used to reconstruct the
  // best cuts sequence starting from the end
  TempVector<int> memory_cuts_idx (n_cuts_max);
  // Positions of the cuts (1..n_samples)
  TempVector<int> memory_cuts_pos (n_cuts_max);
  // Store the counts for 1st bin, {1st+2nd} bin, {1st+2nd+3rd} bin, ...
  TempVector<int> cum_counts_target (nb_factors_target);
  TempVector<int> cum_counts_joint (nb_factors_joint);
  //
  // Stochastic complexity used for the simple BIC cost.
  // Its value is 0.5 * number_of_other_levels where number of other levels is
  // the product of number of unique observed levels of all variables except
  // the variable being optimized.
  //
  double k_bic = 0.5 * (nb_factors_joint_for_cplx - 1) * nb_factors_target;

  for (int cum_bin_idx = 0; cum_bin_idx < n_cuts_max; ++cum_bin_idx) {
    // Compute info from the first bin to the "cum_bin_idx" bin
    // (a sort of cumulative info for bins {0}, {0+1}, {0+1+2}, ...)
    //
    int max_sample_idx_of_cum_bin = sample_idx_breaks_per_cut[cum_bin_idx];
    double weight_per_sample = cum_weights_per_bin[cum_bin_idx] / max_sample_idx_of_cum_bin;
    //
    // Compute the cumulative number of samples until the bin "cum_bin_idx"
    // for the target and discretized joint variables. It is cumulative as
    // it reuses the previous values of cum_counts_target at each iteration
    //
    std::transform (begin (cum_counts_target), end (cum_counts_target),
        coarse_counts_target.row_begin (cum_bin_idx), begin (cum_counts_target),
        std::plus<int>());
    std::transform (begin (cum_counts_joint), end (cum_counts_joint),
        coarse_counts_joint.row_begin (cum_bin_idx), begin (cum_counts_joint),
        std::plus<int>());
    //
    // Compute joint entropy
    //
    double Hk_cum_bin_joint = 0, H_cum_bin_joint = 0;
    for (int level = 0; level < nb_factors_joint; level++) {
      int weighted_count = flag_weights
                               ? int(weight_per_sample * cum_counts_joint[level] + 0.5)
                               : cum_counts_joint[level];
      Hk_cum_bin_joint += cache->getH(weighted_count);
    }
    H_cum_bin_joint = Hk_cum_bin_joint;
    //
    // Compute target entropy
    //
    double Hk_cum_bin_target = 0, H_cum_bin_target = 0;
    for (int level = 0; level < nb_factors_target; level++) {
      int weighted_count = flag_weights
                               ? int(weight_per_sample * cum_counts_target[level] + 0.5)
                               : cum_counts_target[level];
      Hk_cum_bin_target -= cache->getH(weighted_count);
      H_cum_bin_target -= cache->getH(weighted_count);

      if (cplx == 0 && cum_counts_target[level] > 0)
        Hk_cum_bin_target -= k_bic * cache->getLog(n_samples);
      else if (cplx == 1)
        Hk_cum_bin_target -= cache->getLogC(weighted_count, nb_factors_joint_for_cplx);
    }
    //
    // Compute mutual info
    //
    I[cum_bin_idx] = H_cum_bin_joint + H_cum_bin_target;
    Ik[cum_bin_idx] = Hk_cum_bin_joint + Hk_cum_bin_target;
    //
    // Evaluate mutual information for partitions of the bins:
    //
    // i.e.: if cum_bin_idx = 2, we just computed the total entropy and
    // information of the group of bins {0+1+2}, and in previous iterations,
    // we computed the best information for bins {0}, {0+1} and {1}. Here,
    // we will determine if introducing a cut after index 0 or 1 is valuable.
    // So we will evaluate if {0} | {1+2} or {0+1} | [2} contains more
    // information than {0+1+2}. As we already computed the best information of
    // bins {0} and {0+1}, we just need to compute entropy and information of
    // the bins {1+2} and {2}.
    //
    TempAllocatorScope scope;
    TempVector<int> cut_counts_joint (cum_counts_joint);    // copy
    TempVector<int> cut_counts_target (cum_counts_target);  // copy
    //
    // Initialize the cuts index and positions. Before trying to find cuts,
    // solution is one single bin from 0 to cum_bin_idx.
    //
    memory_cuts_idx[cum_bin_idx] = 0;
    memory_cuts_pos[cum_bin_idx] = 0;
    //
    // Move cut_idx iterator on possible cuts
    //
    for (int cut_idx = 0; cut_idx < cum_bin_idx; ++cut_idx) {
      int max_sample_idx_of_cut = sample_idx_breaks_per_cut[cut_idx];
      if (flag_weights)
        weight_per_sample = (cum_weights_per_bin[cum_bin_idx]
                        - cum_weights_per_bin[cut_idx])
                        / (max_sample_idx_of_cum_bin - max_sample_idx_of_cut);
      //
      // Compute the sum of the number of samples from the bin "cut_idx"
      // to the bin "cum_bin_idx" for the target variable. It uses a
      // substraction as we start from the cumulative sum and remove
      // the first bin(s)
      //
      std::transform(begin(cut_counts_target), end(cut_counts_target),
          coarse_counts_target.row_begin(cut_idx), begin(cut_counts_target),
          std::minus<int>());
      std::transform(begin(cut_counts_joint), end(cut_counts_joint),
          coarse_counts_joint.row_begin(cut_idx), begin(cut_counts_joint),
          std::minus<int>());
      //
      // Compute the entropy for a bin from cut_idx to cum_bin_idx
      // i.e. : if cum_bin_idx = 3, we evaluated above the information of
      // the large bin {0,1,2,3}. Now, we test the cuts, so if cut_idx = 0,
      // we need to evaluate the information of {0}+{1,2,3}. Information of {0}
      // has already been computed in a previous iteration. So here, we
      // compute the information of the remaining part {1,2,3}
      //
      double Hk_cut_to_cum_joint = 0, H_cut_to_cum_joint = 0;
      for (int level = 0; level < nb_factors_joint; level++) {
        int weighted_count = flag_weights
                                 ? int(weight_per_sample * cut_counts_joint[level] + 0.5)
                                 : cut_counts_joint[level];
        Hk_cut_to_cum_joint += cache->getH(weighted_count);
      }
      H_cut_to_cum_joint = Hk_cut_to_cum_joint;

      double Hk_cut_to_cum_target = 0, H_cut_to_cum_target = 0;
      for (int level = 0; level < nb_factors_target; level++) {
        int weighted_count = flag_weights
                                 ? int(weight_per_sample * cut_counts_target[level] + 0.5)
                                 : cut_counts_target[level];
        Hk_cut_to_cum_target -= cache->getH(weighted_count);
        H_cut_to_cum_target -= cache->getH(weighted_count);

        if (cplx == 0 && cut_counts_target[level] > 0)
          Hk_cut_to_cum_target -= k_bic * cache->getLog(n_samples);
        else if (cplx == 1)
          Hk_cut_to_cum_target -= cache->getLogC(weighted_count, nb_factors_joint_for_cplx);
      }
      //
      // Compute the information value for a bin from cut_idx to cum_bin_idx
      //
      double I_cut_to_cum_bin = H_cut_to_cum_joint + H_cut_to_cum_target;
      double Ik_cut_to_cum_bin = Hk_cut_to_cum_joint + Hk_cut_to_cum_target;
      if (cplx == 1) {
        // Combinatorial approximation
        Ik_cut_to_cum_bin -= cache->getLogChoose (n_cuts_max - 1,
                                         nb_factors_target_prev - 1)
                                         / (nb_factors_target_prev - 1);
      }
      //
      // Compute the total information value if a cut is introduced:
      //   Ik ( [0 -> (possible other cuts) -> cut_idx] ) (computed in a previous iteration)
      // + Ik ( [cut_idx+1 -> cum_bin_idx] )              (computed just above)
      //
      double Ik_cum_bin_with_cut = Ik[cut_idx] + Ik_cut_to_cum_bin;
      if ( (Ik_cum_bin_with_cut - Ik[cum_bin_idx]) > kEpsilon ) {
        //
        // Better information found, we store the new information value
        // for the interval [0 cum_bin_idx] and memorize the cut
        //
        I[cum_bin_idx] = I[cut_idx] + I_cut_to_cum_bin;
        Ik[cum_bin_idx] = Ik_cum_bin_with_cut;
        //
        // Index and position of the (last) optimal cut
        //
        if (memory_cuts_idx[cut_idx + 1] == 0)
          memory_cuts_idx[cum_bin_idx] = -cut_idx - 1;
        else
          memory_cuts_idx[cum_bin_idx] = cut_idx + 1;
        memory_cuts_pos[cum_bin_idx] = max_sample_idx_of_cut - 1;
      }
    }  // inner loop on cut_idx
  }    // outer loop on cum_bin_idx
  reconstructCutCoarse (memory_cuts_idx, memory_cuts_pos, n_samples, cuts_target);
}

//------------------------------------------------------------------------------
// computeIxy
//------------------------------------------------------------------------------
// Initialize Ik(x,y) with equal bin discretization
// Repeat
//   optimize on x Ik(x,y): Hx - Hxy - k_bic
//   optimize on y Ik(x,y): Hy - Hxy - k_bic
// Until convergence
//------------------------------------------------------------------------------
// Inputs:
// - data: grid2d<int> (2 + nb_ui rows * nb_samples_not_NAs cols),
//   ranked data with NA filtered out
// - data_idx: grid2d<int> (2 + nb_ui rows * nb_samples_not_NAs cols),
//   order of values in data
// - is_continuous: vector<int>  (0=discrete, 1=continuous)
//   is_continuous[0] corresponds to x, and is_continuous[1] to y
// - var_idx: vector<int>, indexes of variables to process in data rows
//   var_idx[0] corresponds to x, and var_idx[1] to y
// - levels: Vector, number of unique values for each variable
// - weights: Vector, samples weight
// - flag_sample_weights: 1 if weights are used, 0 otherelse
// - initbins: number of bins to start with, initialized with
//   min (30, cubic root of nb samples rounded to nearest int)
// - maxbins: maximum number of bins (initialized with the number of samples)
// - cplx: int, 0 or 1
// - negative_info: bool
// - cache: point to shared cache
// Outputs:
// - cuts_info:
//------------------------------------------------------------------------------
InfoBlock computeIxy(const TempGrid2d<int>& data,
    const TempGrid2d<int>& data_idx, const TempVector<int>& is_continuous,
    const TempVector<int>& var_idx, const TempVector<int>& levels,
    const TempVector<double>& weights, bool flag_sample_weights, int initbins,
    int maxbins, int cplx, bool negative_info,
    std::shared_ptr<CtermCache> cache,
    std::shared_ptr<CutPointsInfo> cuts_info = nullptr) {
  TempAllocatorScope scope;

  bool save_cuts = cuts_info && !cuts_info->cutpoints.empty();
  int n_samples = data.n_cols();
  double n_eff = accumulate(begin(weights), end(weights), 0.0);

  TempGrid2d<int> datafactors(2, n_samples);
  TempGrid2d<int> cut(2, maxbins);
  TempVector<int> r(2);
  // Initialize discrete factors and n_levels
  for (int l = 0; l < 2; ++l) {
    auto var = var_idx[l];
    if (is_continuous[var] != 1) {
      r[l] = levels[var];
      copy(data.row_begin(var), data.row_end(var), datafactors.row_begin(l));
    }
  }
  TempVector<int> r_temp(3);
  int rxy{0};
  TempVector<int> xy_factors(n_samples);
  // Find the best initial conditions with the same number of bins (equalfreq)
  // on all continuous variables.
  double best_res{std::numeric_limits<double>::lowest()};
  int best_initbins{initbins};
  int n_levels_min{n_samples};
  for (int l = 0; l < 2; ++l) {
    if (is_continuous[var_idx[l]] == 1)
      n_levels_min = min(n_levels_min, levels[var_idx[l]]);
  }
  int n_test_max = min(min(initbins, 20), n_levels_min);
  for (int test_n_bins = 2; test_n_bins < n_test_max; ++test_n_bins) {
    // Initialize cut, datafactors and r
    resetCutPoints(
        levels, is_continuous, var_idx, 0, 2, test_n_bins, n_samples, cut);
    updateFactors (data, data_idx, cut, is_continuous,
                   var_idx, 0, 2, datafactors, r);

    TempAllocatorScope scope;
    rxy = setJointFactors(datafactors, r, TempVector<int>{0, 1}, xy_factors);

    r_temp.assign({r[0], r[1], rxy});
    InfoBlock res_temp = computeMI(datafactors.getRow(0), datafactors.getRow(1),
        xy_factors, r_temp, n_eff, weights, cache, cplx, 0);
    // All set if both variables are discrete
    if (is_continuous[var_idx[0]] == 0 && is_continuous[var_idx[1]] == 0)
      return res_temp;

    double Ik = res_temp.I - res_temp.k;
    if (Ik > best_res) {
      best_initbins = test_n_bins;
      best_res = Ik;
    }
  }
  // Initialize cut, datafactors and r
  resetCutPoints(
      levels, is_continuous, var_idx, 0, 2, best_initbins, n_samples, cut);
  updateFactors (data, data_idx, cut, is_continuous,
                 var_idx, 0, 2, datafactors, r);

  // Run dynamic optimization with the best initial conditions.
  int iter_max = save_cuts ? cuts_info->cutpoints.n_rows() : kMaxIter;
  // Keep the result of each iteration
  TempVector<double> I_list(iter_max);
  TempVector<double> Ik_list(iter_max);
  double Ixy{0}, Ikxy{0};  // to be returned
  for (int step = 0; step < iter_max; ++step) {
    if (is_continuous[var_idx[0]]) {
      TempAllocatorScope scope;
      // max_X{ I(X;Y) }, X starts with a single bin
      optimizeCutPoints(data.getConstRow(var_idx[0]),
          data_idx.getConstRow(var_idx[0]), datafactors.getRow(1),
          TempVector<int>(n_samples, 0), r[1], 1, r[1], r[0],
          levels[var_idx[0]], weights, flag_sample_weights, maxbins, cache,
          cplx, cut.getRow(0));
    }
    if (is_continuous[var_idx[1]]) {
      TempAllocatorScope scope;
      // max_y{ I(x;y) }, Y starts with a single bin
      optimizeCutPoints(data.getConstRow(var_idx[1]),
          data_idx.getConstRow(var_idx[1]), datafactors.getRow(0),
          TempVector<int>(n_samples, 0), r[0], 1, r[0], r[1],
          levels[var_idx[1]], weights, flag_sample_weights, maxbins, cache,
          cplx, cut.getRow(1));
    }
    updateFactors (data, data_idx, cut, is_continuous,
                   var_idx, 0, 2, datafactors, r);

    if (save_cuts) {
      auto& cp = cuts_info->cutpoints;
      copy(cut.row_begin(0), cut.row_end(0), cp.row_begin(step));
      copy(cut.row_begin(1), cut.row_end(1), cp.row_begin(step) + maxbins);
    }

    TempAllocatorScope scope;
    rxy = setJointFactors(datafactors, r, TempVector<int>{0, 1}, xy_factors);
    r_temp.assign({r[0], r[1], rxy});
    InfoBlock res_temp = computeMI(datafactors.getRow(0), datafactors.getRow(1),
        xy_factors, r_temp, n_eff, weights, cache, cplx, 0);
    // Adding combinatorial term
    if (is_continuous[var_idx[0]] && r[0] > 1) {
      int n_cuts_max = min(maxbins, levels[var_idx[0]]);
      if (r[0] < n_cuts_max)
        res_temp.k += cache->getLogChoose(n_cuts_max - 1, r[0] - 1);
    }
    if (is_continuous[var_idx[1]] && r[1] > 1) {
      int n_cuts_max = min(maxbins, levels[var_idx[1]]);
      if (r[1] < n_cuts_max)
        res_temp.k += cache->getLogChoose(n_cuts_max - 1, r[1] - 1);
    }

    I_list[step] = res_temp.I;
    Ik_list[step] = res_temp.I - res_temp.k;
    bool converged{false};
    for (int i = step - 1; i >= 0; --i) {
      if (fabs(Ik_list[step] - Ik_list[i]) < kEpsilon) {
        converged = true;
        Ixy = accumulate(begin(I_list) + i, begin(I_list) + step, 0.0);
        Ikxy = accumulate(begin(Ik_list) + i, begin(Ik_list) + step, 0.0);
        Ikxy /= (step - i);  // average over the periodic cycle
        Ixy /= (step - i);
        break;
      }
    }

    if (save_cuts) cuts_info->n_iterations = step+1;
    if (converged) break;

    Ixy = res_temp.I;
    Ikxy = res_temp.I - res_temp.k;
    // Already optimal if any one of them is discrete
    if (!is_continuous[var_idx[0]] || !is_continuous[var_idx[1]]) break;
  }  // for step

  // I and Ik can always be 0 by choosing 1 bin on either X or Y.
  if (!negative_info && Ikxy < 0) {
    Ixy = 0;
    Ikxy = 0;
  }
  if (save_cuts) {
    cuts_info->Ik = Ikxy / n_eff;
    cuts_info->I = Ixy / n_eff;
    cuts_info->I_equal_freq_max = best_res;
  }

  return InfoBlock(n_samples, Ixy, Ixy - Ikxy);
}

//------------------------------------------------------------------------------
// computeIxyui
//------------------------------------------------------------------------------
InfoBlock computeIxyui(const TempGrid2d<int>& data,
    const TempGrid2d<int>& data_idx, const TempVector<int>& is_continuous,
    const TempVector<int>& var_idx, const TempVector<int>& levels,
    const TempVector<double>& weights, bool flag_sample_weights, int initbins,
    int maxbins, int cplx, bool negative_info,
    std::shared_ptr<CtermCache> cache,
    std::shared_ptr<CutPointsInfo> cuts_info = nullptr) {
  TempAllocatorScope scope;

  bool save_cuts = cuts_info && !cuts_info->cutpoints.empty();
  int n_samples = data.n_cols();
  int n_nodes = var_idx.size();
  int n_ui = n_nodes - 2;
  double n_eff = accumulate(begin(weights), end(weights), 0.0);
  int u_initbins = min(30, max(2, int(0.5 + pow(n_eff, 1.0/(n_nodes)))));
  // allocation factors  x y
  TempGrid2d<int> datafactors(n_nodes, n_samples);
  TempGrid2d<int> cut(n_nodes, maxbins);
  TempGrid2d<int> cut_y_u(n_nodes, maxbins);
  TempGrid2d<int> cut_x_u(n_nodes, maxbins);
  TempVector<int> r(n_nodes);  // n_levels of optimized variables

  // Initialize discrete factors and n_levels
  for (int l = 0; l < n_nodes; ++l) {
    auto var = var_idx[l];
    if (is_continuous[var] != 1) {
      r[l] = levels[var];
      copy(data.row_begin(var), data.row_end(var), datafactors.row_begin(l));
    }
  }
  TempGrid2d<int> uyxfactors(4, n_samples);  // u, uy, ux, uyx
  TempVector<int> ruyx(4);
  // Find the best initial conditions with the same number of bins (equalfreq)
  // on all continuous variables.
  double best_res{std::numeric_limits<double>::lowest()};
  int best_initbins{initbins};
  int n_levels_min{n_samples};
  for (int l = 0; l < n_nodes; ++l) {
    if (is_continuous[var_idx[l]] == 1)
      n_levels_min = min(n_levels_min, levels[var_idx[l]]);
  }
  int n_test_max = min(min(initbins, 20), n_levels_min);

  // FRS 4 jan 2024: remove fix to limit number of joint factors
  // if (std::pow (n_test_max-1, n_ui) >= INT_MAX)
  //   {
  //   n_test_max = std::pow (INT_MAX, 1.0 / n_ui) + 1;
  //   Rcpp::Rcout << "Note: Initial number of bins has been limited to "
  //     << n_test_max-1 << " for " << n_ui << " contributors to avoid overflow\n";
  //   }
  TempVector<int> r_temp(3);
  InfoBlock res_temp;
  for (int test_n_bins = 2; test_n_bins < n_test_max; ++test_n_bins) {
    // Initialize factors, cut and r
    resetCutPoints(levels, is_continuous, var_idx, 0, n_nodes, test_n_bins,
        n_samples, cut);
    updateFactors (data, data_idx, cut, is_continuous,
                   var_idx, 0, n_nodes, datafactors, r);

    setUyxJointFactors(datafactors, r, -1, uyxfactors, ruyx);
    r_temp.assign({r[1], ruyx[2], ruyx[3]});
    res_temp = computeMI(datafactors.getRow(1), uyxfactors.getRow(2),
        uyxfactors.getRow(3), r_temp, n_eff, weights, cache, cplx, 1);
    double Ik_y_xu = res_temp.I - res_temp.k;

    r_temp.assign({r[0], ruyx[1], ruyx[3]});
    res_temp = computeMI(datafactors.getRow(0), uyxfactors.getRow(1),
        uyxfactors.getRow(3), r_temp, n_eff, weights, cache, cplx, 1);
    double Ik_x_yu = res_temp.I - res_temp.k;

    if ((Ik_y_xu + Ik_x_yu) > best_res) {
      best_initbins = test_n_bins;
      best_res = (Ik_y_xu + Ik_x_yu);
    }
  }
  // Initialize X and Y cuts with best_initbins
  resetCutPoints(
      levels, is_continuous, var_idx, 0, 2, best_initbins, n_samples, cut);
  updateFactors (data, data_idx, cut, is_continuous,
                 var_idx, 0, 2, datafactors, r);
  TempVector<int> r_old(r);  // n_levels in the previous iteration

  bool reuse_x_u_cuts = false;
  bool reuse_y_u_cuts = false;

  int iter_max = save_cuts ? cuts_info->cutpoints.n_rows() : kMaxIter;
  // Keep the result of each iteration
  TempVector<double> I_list(iter_max);
  TempVector<double> Ik_list(iter_max);
  double Ixy_ui{0}, Ikxy_ui{0};  // to be returned
  // Run optimization with best initial equal freq.
  for (int step = 0; step < iter_max; ++step) {
    // optimize I(y;xu) over x and u
    // Either reset cutpoints or re-use I(y;u) cutpoints
    if(reuse_y_u_cuts){
      for(int l = 2; l < n_nodes; ++l){
        copy(cut_y_u.row_begin(l), cut_y_u.row_end(l), cut.row_begin(l));
      }
    } else {
      resetCutPoints(
          levels, is_continuous, var_idx, 2, n_nodes, u_initbins,
          n_samples, cut);
    }
    updateFactors (data, data_idx, cut, is_continuous,
                   var_idx, 2, n_nodes, datafactors, r);
    copy(begin(r) + 2, end(r), begin(r_old) + 2);
    for (int count = 0; (count < kMaxIterOnU) && !reuse_y_u_cuts; ++count) {
      for (int l = 2; l < n_nodes; ++l) {
        if (is_continuous[var_idx[l]] == 0) continue;
        // opt u, I(y;xu)
        setUyxJointFactors(datafactors, r_old, l, uyxfactors, ruyx);
        // init variables for the optimization run
        int r0 = ruyx[3];  // xyu
        int r1 = ruyx[2];  // xu
        int sc_levels1 = r_old[1];
        int sc_levels2 = r_old[l];  // old nlevels for combinatorial term
        // Run optimization on U. 2 factors xyu and xu
        optimizeCutPoints(data.getConstRow(var_idx[l]),
            data_idx.getConstRow(var_idx[l]), uyxfactors.getRow(3),
            uyxfactors.getRow(2), r0, r1, sc_levels1, sc_levels2,
            levels[var_idx[l]], weights, flag_sample_weights, maxbins, cache,
            cplx, cut.getRow(l));
      }  // for all Uis
      updateFactors (data, data_idx, cut, is_continuous,
                     var_idx, 2, n_nodes, datafactors, r);
      copy(begin(r) + 2, end(r), begin(r_old) + 2);

      if (n_ui == 1) break;
    }  // Iteration on ui

    setUyxJointFactors(datafactors, r_old, -1, uyxfactors, ruyx);
    r_temp.assign({r_old[1], ruyx[2], ruyx[3]});
    res_temp = computeMI(datafactors.getRow(1), uyxfactors.getRow(2),
        uyxfactors.getRow(3), r_temp, n_eff, weights, cache, cplx, 1);
    double I_y_xu = res_temp.I;  // Before optimization on X.
    double Ik_y_xu = res_temp.I - res_temp.k;
    if (is_continuous[var_idx[0]] && r_old[0] > 1) {
      int n_cuts_max = min(levels[var_idx[0]], maxbins);
      if (r_old[0] < n_cuts_max)
        Ik_y_xu -= cache->getLogChoose(n_cuts_max - 1, r_old[0] - 1);
    }

    if (is_continuous[var_idx[0]] == 1) {
      // I(y;xu), optimize on x
      int r0 = ruyx[1];  // uy
      int r1 = ruyx[0];  // u
      int sc_levels1 = r_old[1];
      int sc_levels2 = r_old[0];
      // Run optimization on X. 2 factors uy and u
      optimizeCutPoints(data.getConstRow(var_idx[0]),
          data_idx.getConstRow(var_idx[0]), uyxfactors.getRow(1),
          uyxfactors.getRow(0), r0, r1, sc_levels1, sc_levels2,
          levels[var_idx[0]], weights, flag_sample_weights, maxbins, cache,
          cplx, cut.getRow(0));
    }

    // optimize I(x;yu) over y and u
    // Either reset cutpoints or re-use I(x;u) cutpoints
    if(reuse_x_u_cuts){
      for(int l = 2; l < n_nodes; ++l){
        copy(cut_x_u.row_begin(l), cut_x_u.row_end(l), cut.row_begin(l));
      }
    } else {
      resetCutPoints(
          levels, is_continuous, var_idx, 2, n_nodes, u_initbins,
          n_samples, cut);
    }
    updateFactors (data, data_idx, cut, is_continuous,
                   var_idx, 2, n_nodes, datafactors, r);
    copy(begin(r) + 2, end(r), begin(r_old) + 2);
    for (int count = 0; (count < kMaxIterOnU) && !reuse_x_u_cuts; ++count) {
      for (int l = 2; l < n_nodes; ++l) {
        if (is_continuous[var_idx[l]] == 0) continue;

        setUyxJointFactors(datafactors, r_old, l, uyxfactors, ruyx);
        // init variables for the optimization run
        int r0 = ruyx[3];  // xyu
        int r1 = ruyx[1];  // yu
        int sc_levels1 = r_old[0];
        int sc_levels2 = r_old[l];
        // Run optimization on U. 2 factors xyu and yu
        optimizeCutPoints(data.getConstRow(var_idx[l]),
            data_idx.getConstRow(var_idx[l]), uyxfactors.getRow(3),
            uyxfactors.getRow(1), r0, r1, sc_levels1, sc_levels2,
            levels[var_idx[l]], weights, flag_sample_weights, maxbins, cache,
            cplx, cut.getRow(l));
      }  // for all Uis
      updateFactors (data, data_idx, cut, is_continuous,
                     var_idx, 2, n_nodes, datafactors, r);
      copy(begin(r) + 2, end(r), begin(r_old) + 2);

      if (n_ui == 1) break;
    }  // Iteration on ui

    setUyxJointFactors(datafactors, r_old, -1, uyxfactors, ruyx);
    r_temp.assign({r_old[0], ruyx[1], ruyx[3]});
    res_temp = computeMI(datafactors.getRow(0), uyxfactors.getRow(1),
        uyxfactors.getRow(3), r_temp, n_eff, weights, cache, cplx, 1);
    double I_x_yu = res_temp.I;  // Before updating Y (and X).
    double Ik_x_yu = res_temp.I - res_temp.k;
    if ((is_continuous[var_idx[1]] == 1) && (r_old[1] > 1)) {
      int n_cuts_max = min(levels[var_idx[1]], maxbins);
      if (r_old[1] < n_cuts_max) {
        Ik_x_yu -= cache->getLogChoose(n_cuts_max - 1, r_old[1] - 1);
      }
    }

    if (is_continuous[var_idx[1]] == 1) {
      // I(x;yu), optimize on y
      int r0 = ruyx[2];  // ux
      int r1 = ruyx[0];  // u
      int sc_levels1 = r_old[0];
      int sc_levels2 = r_old[1];
      // Run optimization on Y. 2 factors ux and u
      optimizeCutPoints(data.getConstRow(var_idx[1]),
          data_idx.getConstRow(var_idx[1]), uyxfactors.getRow(2),
          uyxfactors.getRow(0), r0, r1, sc_levels1, sc_levels2,
          levels[var_idx[1]], weights, flag_sample_weights, maxbins, cache,
          cplx, cut.getRow(1));
    }

    // optimize I(x;u) over u
    // Either reset cutpoints or re-use I(y;u) cutpoints
    if(reuse_x_u_cuts){
      for(int l = 2; l < n_nodes; ++l){
        copy(cut_x_u.row_begin(l), cut_x_u.row_end(l), cut.row_begin(l));
      }
    } else {
      resetCutPoints(
          levels, is_continuous, var_idx, 2, n_nodes, u_initbins,
          n_samples, cut);
    }
    updateFactors (data, data_idx, cut, is_continuous,
                   var_idx, 2, n_nodes, datafactors, r);
    copy(begin(r) + 2, end(r), begin(r_old) + 2);
    for (int count = 0; (count < kMaxIterOnU) && !reuse_x_u_cuts; ++count) {
      for (int l = 2; l < n_nodes; ++l) {
        if (is_continuous[var_idx[l]] == 0) continue;

        setUyxJointFactors(datafactors, r_old, l, uyxfactors, ruyx);
        // init variables for the optimization run
        int r0 = ruyx[2];           // xu
        int r1 = ruyx[0];           // u
        int sc_levels1 = r_old[0];  // x
        int sc_levels2 = r_old[l];  // u
        // optimization run on var_idx[l], 2 factors xu and u
        optimizeCutPoints(data.getConstRow(var_idx[l]),
            data_idx.getConstRow(var_idx[l]), uyxfactors.getRow(2),
            uyxfactors.getRow(0), r0, r1, sc_levels1, sc_levels2,
            levels[var_idx[l]], weights, flag_sample_weights, maxbins, cache,
            cplx, cut.getRow(l));
      }  // for all Uis
      updateFactors (data, data_idx, cut, is_continuous,
                     var_idx, 2, n_nodes, datafactors, r);
      copy(begin(r) + 2, end(r), begin(r_old) + 2);

      if (n_ui == 1) break;
    }  // Iteration on ui

    //save u cutpoints
    for(int l = 2; l < n_nodes; ++l){
      copy(cut.row_begin(l), cut.row_end(l), cut_x_u.row_begin(l));
    }
    setUyxJointFactors(datafactors, r_old, -1, uyxfactors, ruyx);
    r_temp.assign({r_old[0], ruyx[0], ruyx[2]});
    res_temp = computeMI(datafactors.getRow(0), uyxfactors.getRow(0),
        uyxfactors.getRow(2), r_temp, n_eff, weights, cache, cplx, 1);
    double I_x_u = res_temp.I;  // After optimization on U.
    double Ik_x_u = res_temp.I - res_temp.k;

    // optimize I(y;u) over u
    // Either reset cutpoints or re-use I(y;u) cutpoints
    if(reuse_y_u_cuts){
      for(int l = 2; l < n_nodes; ++l){
        copy(cut_y_u.row_begin(l), cut_y_u.row_end(l), cut.row_begin(l));
      }
    } else {
      resetCutPoints(
          levels, is_continuous, var_idx, 2, n_nodes, u_initbins,
          n_samples, cut);
    }
    updateFactors (data, data_idx, cut, is_continuous,
                   var_idx, 2, n_nodes, datafactors, r);
    copy(begin(r) + 2, end(r), begin(r_old) + 2);
    for (int count = 0; (count < kMaxIterOnU) && !reuse_y_u_cuts; ++count) {
      for (int l = 2; l < n_nodes; ++l) {
        if (is_continuous[var_idx[l]] == 0) continue;

        setUyxJointFactors(datafactors, r_old, l, uyxfactors, ruyx);
        int r0 = ruyx[1];           // yu
        int r1 = ruyx[0];           // u
        int sc_levels1 = r_old[1];  // y
        int sc_levels2 = r_old[l];  // u
        // optimization run on var_idx[l], 2 factors yu and u
        optimizeCutPoints(data.getConstRow(var_idx[l]),
            data_idx.getConstRow(var_idx[l]), uyxfactors.getRow(1),
            uyxfactors.getRow(0), r0, r1, sc_levels1, sc_levels2,
            levels[var_idx[l]], weights, flag_sample_weights, maxbins, cache,
            cplx, cut.getRow(l));
      }  // for all Uis
      updateFactors (data, data_idx, cut, is_continuous,
                     var_idx, 2, n_nodes, datafactors, r);
      copy(begin(r) + 2, end(r), begin(r_old) + 2);

      if (n_ui == 1) break;
    }  // Iteration on ui

    //save u cutpoints
    for(int l = 2; l < n_nodes; ++l){
      copy(cut.row_begin(l), cut.row_end(l), cut_y_u.row_begin(l));
    }
    setUyxJointFactors(datafactors, r_old, -1, uyxfactors, ruyx);
    r_temp.assign({r_old[1], ruyx[0], ruyx[1]});
    res_temp = computeMI(datafactors.getRow(1), uyxfactors.getRow(0),
        uyxfactors.getRow(1), r_temp, n_eff, weights, cache, cplx, 1);
    double I_y_u = res_temp.I;  // After optimization on U.
    double Ik_y_u = res_temp.I - res_temp.k;

    // End of iteration: update X and Y cuts
    updateFactors (data, data_idx, cut, is_continuous,
                   var_idx, 0, 2, datafactors, r);
    copy(begin(r), end(r), begin(r_old));

    if (save_cuts) {
      auto& cp = cuts_info->cutpoints;
      copy(cut.row_begin(0), cut.row_end(0), cp.row_begin(step));
      copy(cut.row_begin(1), cut.row_end(1), cp.row_begin(step) + maxbins);
    }

    // Compute I(X;Y|U)
    double part_one = I_x_yu - I_x_u;
    double part_two = I_y_xu - I_y_u;
    // Either sum may be negative when we can't find U bins on I(X;YU) (or
    //I(Y;XY)) that are as good as those found on I(X;U) (I(Y;U)). If that is
    //the case, we reuse the better U bins for all terms for the next iteration.
    reuse_x_u_cuts = part_one < kPrecision || (reuse_x_u_cuts && (part_two < kPrecision));
    reuse_y_u_cuts = part_two < kPrecision || (reuse_y_u_cuts && (part_one < kPrecision));

    if(reuse_x_u_cuts || reuse_y_u_cuts) {
      bool one_bin = r[0] == 1 || r[1] == 1;
      if(one_bin && (step != (iter_max-1))) {
        // Next iteration will have 0 information with no further optimization possible.
        iter_max = step+2;
        if (save_cuts) cuts_info->n_iterations = step+1;
      }
      if(!one_bin) continue;
      // With one bin, we can keep the result and check for convergence
    }

    double cond_I = 0.5 * (part_one + part_two);
    double cond_Ik = 0.5 * (Ik_x_yu - Ik_x_u + Ik_y_xu - Ik_y_u);
    I_list[step] = cond_I;
    Ik_list[step] = cond_Ik;

    // Test stop condition on stop1
    bool converged{false};
    for (int i = step - 1; i >= 0; i--) {
      // If no real improvement over last information
      if (fabs(cond_Ik - Ik_list[i]) < kEpsilon) {
        converged = true;
        Ixy_ui = accumulate(begin(I_list) + i, begin(I_list) + step, 0.0);
        Ikxy_ui = accumulate(begin(Ik_list) + i, begin(Ik_list) + step, 0.0);
        Ikxy_ui /= (step - i);  // average over the periodic cycle
        Ixy_ui /= (step - i);
        break;
      }
    }
    if (save_cuts) cuts_info->n_iterations = step;
    if (converged) break;

    Ixy_ui = cond_I;
    Ikxy_ui = cond_Ik;
  }  // for step
  // I and Ik can always be 0 by choosing 1 bin on either X or Y.
  if (!negative_info && Ikxy_ui < 0) {
    Ixy_ui = 0;
    Ikxy_ui = 0;
  }
  if (save_cuts) {
    cuts_info->Ik = Ikxy_ui / n_eff;
    cuts_info->I = Ixy_ui / n_eff;
  }

  return InfoBlock(n_samples, Ixy_ui, Ixy_ui - Ikxy_ui);
}

}  // anonymous namespace

InfoBlock computeCondMutualInfo(const TempGrid2d<int>& data,
    const TempGrid2d<int>& data_idx, const TempVector<int>& levels,
    const TempVector<int>& is_continuous, const TempVector<int>& var_idx,
    const TempVector<double>& sample_weights, bool flag_sample_weights,
    int initbins, int maxbins, int cplx, bool negative_info,
    std::shared_ptr<CtermCache> cache,
    std::shared_ptr<CutPointsInfo> cuts_info) {
  if (data.n_rows() == 2) {
    return computeIxy(data, data_idx, is_continuous, var_idx, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx,
        negative_info, cache, cuts_info);
  } else {
    return computeIxyui(data, data_idx, is_continuous, var_idx, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx,
        negative_info, cache, cuts_info);
  }
}

// compute Rscore and three point mutual information I(x;y;z | u)
// return Info3PointBlock{score, N * Ixyz_ui, N * kxyz_ui}
Info3PointBlock computeInfo3PointAndScore(const TempGrid2d<int>& data,
    const TempGrid2d<int>& data_idx, const TempVector<int>& levels,
    const TempVector<int>& is_continuous, const TempVector<int>& var_idx,
    const TempVector<double>& sample_weights, bool flag_sample_weights,
    int initbins, int maxbins, int cplx, bool negative_info,
    std::shared_ptr<CtermCache> cache) {
  TempAllocatorScope scope;

  int n_ui = data.n_rows() - 3;
  // Optimize variables for each MI estimation for the R score
  // I(x,y|u,z)
  InfoBlock res_temp = computeIxyui(data, data_idx, is_continuous, var_idx,
      levels, sample_weights, flag_sample_weights, initbins, maxbins, cplx,
      negative_info, cache);
  double I_xy_zu = res_temp.I;
  double Ik_xy_zu = res_temp.I - res_temp.k;

  TempVector<int> var_idx_t(begin(var_idx), begin(var_idx) + n_ui + 2);
  // Do opt run on I(X;Y|U)
  if (n_ui > 0) {
    res_temp = computeIxyui(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx,
        negative_info, cache);
  } else {
    res_temp = computeIxy(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx,
        negative_info, cache);
  }
  double I_xy_u = res_temp.I;
  double Ik_xy_u = res_temp.I - res_temp.k;

  // I(x,z|u)
  var_idx_t[1] = var_idx.back();  // Z
  if (n_ui > 0) {
    res_temp = computeIxyui(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx,
        negative_info, cache);
  } else {
    res_temp = computeIxy(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx,
        negative_info, cache);
  }
  double Ik_xz_u = res_temp.I - res_temp.k;

  // I(y,z|u)
  var_idx_t[0] = var_idx[1];  // Y
  if (n_ui > 0) {
    res_temp = computeIxyui(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx,
        negative_info, cache);
  } else {
    res_temp = computeIxy(data, data_idx, is_continuous, var_idx_t, levels,
        sample_weights, flag_sample_weights, initbins, maxbins, cplx,
        negative_info, cache);
  }
  double Ik_yz_u = res_temp.I - res_temp.k;

  double xz = Ik_xz_u - Ik_xy_u;
  double yz = Ik_yz_u - Ik_xy_u;
  // Data processing inequality
  double dpi = std::fmin(xz, yz) - std::log1p(exp(-std::fabs(xz - yz)));

  // Conditional three point information
  double I_xyz_u = I_xy_u - I_xy_zu;
  // Ik_xyz_u can be seen as the probability of non-v-structure
  double Ik_xyz_u = Ik_xy_u - Ik_xy_zu;

  double Rscore = std::fmin(dpi, Ik_xyz_u);
  return Info3PointBlock{
      Rscore, I_xyz_u, I_xyz_u - Ik_xyz_u, I_xy_u, I_xy_u - Ik_xy_u};
}

}  // namespace computation
}  // namespace miic
