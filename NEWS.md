# Development version

# v1.5.1

## Fixes and improvements
- Fix various bugs in the computation of information in the presence of NA
  values in the dataset.

- An overhaul of the C++ code base, better memory management, computation time
  and code readability.

# v1.5.0

## Features
- Add a column `consensus` to the reconstructed graph's edges summary associated
  with the option `consistent`, and a new parameter `consensus_threshold`
  accordingly.

- Add a parameter `ori_proba_ratio` to have more control on the orientation of
  edges.

## Fixes and improvements
- Faster post processing in R

- Rework plot functionality

- Fix a bug in the orientation part about the log score

- Refactor of the C++ code base (orientation).

# v1.4.2

## Fixes and improvements
- Various fixes of memory leaks and ambiguous function calls (at least for all
  that appear in CRAN check).

- Refactor of the C++ code base (confidence cut).

## Known issues
- Error when running the cosmicCancer example on CRAN's Solaris system.

## Miscellaneous
- Move from BitBucket to GitHub, the repo is now public.

# v1.4.1

Not released. A failed attempt to fix the check errors (though submitted to
CRAN).

# v1.4.0

## Incompatible changes
- Standardize the API naming convention: `snake_case` for parameters and
  `camelCase` for functions. This should have led to a major version increment to
  v2.0.0 given the previous version on CRAN is v1.0.3. But v1.0.3 and earlier
  versions were not properly maintained and versioned under a version control
  system (so we actually forgot to take them into consideration when releasing
  this version).

## Features
- The method now works with continuous variables (solely or mixed with discrete
  variables), thanks to the discretization method as described in
  [Cabeli et al., PLoS Comp. Bio. 2020](https://doi.org/10.1371/journal.pcbi.1007866).

- Add an option `consistent` to improve the reconstructed graph's
  interpretability based on schemes as described in
  [Li et al., NeurIPS 2019](https://papers.nips.cc/paper/9573-constraint-based-causal-structure-learning-with-consistent-separating-sets).

## Fixes and improvements
- Various fixes of memory leaks and typos.

- Major refactoring of the old C++ code base (still WIP) to improve readability
  and flexibility, and to enforce proper coding style and documentation.

- Enforce proper coding style for the R code base.

## Known issues
- Still have some memory leaks and CRAN check errors and notes on certain
  platforms.
