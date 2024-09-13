# v2.0.1

## Features

* Release to CRAN.

## Fixes and improvements

* Faster post-processing in R for datasets with large number of variables.

## Breaking changes

Consolidating long-pending breaking changes:

* Harmonization of exported function names using `camel case`.

* Harmonization of parameters and return values using `snake case`.

* Harmonization of abbreviations.

All the documentation has been updated accordingly, if you encounter any issue
upgrading to this version, please consult the help of the relevant function
for more information about its interface.

For the core `miic()` function, the main breaking changes in the interface
(when upgrading from the 1.5.3 release on CRAN) are:

in the parameters:

* `cplx`: renaming of the complexity term `"mdl"` &rarr; `"bic"`

* `ori_proba_ratio` &rarr; `ort_proba_ratio`  
  
in the miic object returned:
 
* `all.edges.summary` &rarr; `summary`
  * `Nxy_ai` &rarr; `n_xy_ai`
  * `log_confidence` &rarr; `info_shifted`
  * `infOrt` &rarr; `ort_inferred`
  * `trueOrt` &rarr; `ort_ground_truth`
  * `isOrtOk` &rarr; `is_inference_correct`
  * `isCausal` &rarr; `is_causal`
  * `proba` &rarr; `p_y2x`, `p_x2y`
  * `consensus` &rarr; `ort_consensus`

* `orientations.prob` &rarr; `triples`
  * `NI3` &rarr; `ni3`
  * `Error` &rarr; `conflict`
  
Still compared to 1.5.3, another important change in the behavior of `miic()`
is that, by default, `miic()` no longer propagates orientations 
and allows latent variables discovery during orientation step.

## Known issues

* Conditioning on a (very) large number of contributors can lead to a memory 
  fault.
  
# v2.0.0

## Features

* tMIIC version for temporal causal discovery on stationary time series:
  new mode of `miic()` to reconstruct networks from temporal stationary datasets.
  [Simon et al., eLife, reviewed preprint]
  (https://www.biorxiv.org/content/10.1101/2024.02.06.579177v1.abstract)  
  The temporal mode of `miic()` is not activated by default and can be enabled by
  setting the newly added parameter `mode` to `"TS"`(Temporal Stationary).
  A tuning of the temporal mode is possible through a set of new parameters:
  `max_nodes`, `n_layers`, `delta_t`, `mov_avg` and `keep_max_data`.

# v1.8.1

## Fixes and improvements

* The discretization of continuous variables has been improved when dealing 
  with variables having a large number of identical values.
  
* Fix for memory overflow on shared memory space.

# v1.8.0

## Features

* Addition of the 'is consequence' prior knowledge. Consequence variables are
  excluded from the possible contributors, edges between consequences are
  ignored and edges between a non consequence and a consequence are pre-oriented
  toward the consequence.  
  Information about consequence variables can be provided to `miic()`
  in the `state_order`, by supplying an `is_consequence` column.

# v1.7.0

## Features

* iMIIC version introducing contextual variables, genuine vs putative causes
  and multiple enhancements to deal with very large datasets.
  [Ribeiro-Dantas et al., iScience 2024]
  (https://arxiv.org/abs/2303.06423)  
  Information on contextual variables can be provided to `miic()`
  in the `state_order`, by supplying an `is_contextual`column and 
  genuine vs putative causes can be tuned by the newly added parameter
  `ort_consensus_ratio`.
  
# v1.6.0

## Features

* Enhancement of orientations using mutual information supremum principle for 
  finite datasets.
  [Cabeli et al., Why21 at NeurIPS 2021]
  (http://kinefold.curie.fr/isambertlab/PAPERS/cabeli_Why21-NeurIPS2021.pdf)  
  The use of enhanced orientations is controlled by the newly added parameter
  `negative_info` of `miic()` and is activated by default.

* By default, `miic()` does not propagate orientations anymore
  and allows latent variables discovery during orientation step. 

# v1.5.3

## Features

* Release to CRAN

# v1.5.2

## Fixes and improvements
* Further refactoring of the C++ code for the computation of information.

* Fix minor bugs in the continuous computation.

* Fix incompatibility with older versions of GCC (std::align).

# v1.5.1

## Fixes and improvements
* Fix various bugs in the computation of information in the presence of NA
  values in the dataset.

* An overhaul of the C++ code base, better memory management, computation time
  and code readability.

# v1.5.0

## Features
* Add a column `consensus` to the reconstructed graph's edges summary associated
  with the option `consistent`, and a new parameter `consensus_threshold`
  accordingly.

* Add a parameter `ori_proba_ratio` to have more control on the orientation of
  edges.

## Fixes and improvements
* Faster post processing in R.

* Rework plot functionality.

* Fix a bug in the orientation part about the log score.

* Refactor of the C++ code base (orientation).

# v1.4.2

## Fixes and improvements
* Various fixes of memory leaks and ambiguous function calls (at least for all
  that appear in CRAN check).

* Refactor of the C++ code base (confidence cut).

## Known issues
* Error when running the cosmicCancer example on CRAN's Solaris system.

## Miscellaneous
* Move from BitBucket to GitHub, the repo is now public.

# v1.4.1

Not released. A failed attempt to fix the check errors (though submitted to
CRAN).

# v1.4.0

## Incompatible changes
* Standardize the API naming convention: `snake_case` for parameters and
  `camelCase` for functions. This should have led to a major version increment to
  v2.0.0 given the previous version on CRAN is v1.0.3. But v1.0.3 and earlier
  versions were not properly maintained and versioned under a version control
  system (so we actually forgot to take them into consideration when releasing
  this version).

## Features
* The method now works with continuous variables (solely or mixed with discrete
  variables), thanks to the discretization method as described in
  [Cabeli et al., PLoS Comput. Biol. 2020](https://doi.org/10.1371/journal.pcbi.1007866).

* Add an option `consistent` to improve the reconstructed graph's
  interpretability based on schemes as described in
  [Li et al., NeurIPS 2019](https://papers.nips.cc/paper/9573-constraint-based-causal-structure-learning-with-consistent-separating-sets).

## Fixes and improvements
* Various fixes of memory leaks and typos.

* Major refactoring of the old C++ code base (still WIP) to improve readability
  and flexibility, and to enforce proper coding style and documentation.

* Enforce proper coding style for the R code base.

## Known issues
* Still have some memory leaks and CRAN check errors and notes on certain
  platforms.
