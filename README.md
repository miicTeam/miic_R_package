# MIIC
  <!-- badges: start -->
  [![CRAN
  Status](https://www.r-pkg.org/badges/version/miic)](https://cran.r-project.org/package=miic)
  [![R build
  status](https://github.com/miicTeam/miic_R_package/workflows/R-CMD-check/badge.svg)](https://github.com/miicTeam/miic_R_package/actions)
  <!-- badges: end -->

This repository contains the source code for MIIC (**M**ultivariate **I**nformation based **I**nductive **C**ausation), a method based on constraint-based approaches that learns a large class of causal or non-causal graphical models from purely observational data while including the effects of unobserved latent variables. Starting from a complete graph, the method iteratively removes dispensable edges, by uncovering significant information contributions from indirect paths, and assesses edge-specific confidences from randomization of available data. The remaining edges are then oriented based on the signature of causality in observational data. This approach can be applied on a wide range of datasets and provide new biological insights on regulatory networks from single cell expression data, genomic alterations during tumor development and co-evolving residues in protein structures. For more information you can refer to Cabeli et al. PLoS Comp. Bio. 2020, Verny et al. PLoS Comp. Bio. 2017.

## References
Cabeli V., Verny L., Sella N., Uguzzoni G., Verny M., Isambert H.; Learning clinical networks from medical records based on information estimates in mixed-type data; PLoS computational biology., 2020. [doi:10.1371/journal.pcbi.1007866](https://doi.org/10.1371/journal.pcbi.1007866) | [code](https://github.com/vcabeli/miic_PLoS)

Li H., Cabeli V., Sella N., Isambert H.; Constraint-based causal structure learning with consistent separating sets; [In Advances in Neural Information Processing Systems 2019.](https://papers.nips.cc/paper/9573-constraint-based-causal-structure-learning-with-consistent-separating-sets) | [code](https://github.com/honghaoli42/consistent_pcalg)

Verny L., Sella N., Affeldt S., Singh PP., Isambert H.; Learning causal networks with latent variables from multivariate information in genomic data;  PLoS Comput. Biol., 2017. [doi:10.1371/journal.pcbi.1005662](https://doi.org/10.1371/journal.pcbi.1005662)

## Prerequisites
MIIC contains R and C++ sources.
- To compile from source, a compiler with support for c++14 language features is required.
- MIIC imports the following R packages: ppcor, scales, stats, Rcpp

## Installation

From CRAN (release):
```R
install.packages("miic")
```
Or from GitHub (development):
```R
# install.packages("remotes")
remotes::install_github("miicTeam/miic_R_package")
```

## Quick start

MIIC allows you to create a graph object from a dataset of observations of both discrete and continuous variables, potentially with missing values and taking into account unobserved latent variables.
You can find this example along others by calling the documentation of the main function `?miic` from R.
```R
library(miic)

# EXAMPLE HEMATOPOIESIS
data(hematoData)
# execute MIIC (reconstruct graph)
miic.res <- miic(
  input_data = hematoData, latent = "yes",
  n_shuffles = 10, conf_threshold = 0.001
)

# plot graph with igraph
if(require(igraph)) {
  plot(miic.res, method="igraph")
}
```

## Documentation
You can find the documentation pages in the "man" folder, in the auto generated [PDF](https://cran.r-project.org/web/packages/miic/miic.pdf), or use R functions `help()` and `?`.

## Authors
- Vincent Cabeli
- Honghao Li
- Marcel Ribeiro Dantas
- Verny Louis
- Sella Nadir
- Séverine Affeldt
- Hervé Isambert

## License
GPL-2 | GPL-3
