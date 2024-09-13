# MIIC
  <!-- badges: start -->
  [![CRAN
  Status](https://www.r-pkg.org/badges/version/miic)](https://cran.r-project.org/package=miic)
  [![R build
  status](https://github.com/miicTeam/miic_R_package/workflows/R-CMD-check/badge.svg)](https://github.com/miicTeam/miic_R_package/actions)
  <!-- badges: end -->

This repository contains the source code for MIIC (**M**ultivariate 
**I**nformation based **I**nductive **C**ausation), a method based on 
constraint-based approaches that learns a large class of causal 
or non-causal graphical models from purely observational data 
while including the effects of unobserved latent variables. 
Starting from a complete graph, the method iteratively removes dispensable 
edges, by uncovering significant information contributions from indirect paths, 
and assesses edge-specific confidences from randomization of available data. 
The remaining edges are then oriented based on the signature of causality in 
observational data. This approach can be applied on a wide range of datasets 
and provide new biological insights on regulatory networks from single cell 
expression data, genomic alterations during tumor development and co-evolving 
residues in protein structures. Since the version 2.0, MIIC can 
in addition process stationary time series to unveil temporal causal graphs.

## References
Simon F., Comes M. C., Tocci T., Dupuis L., Cabeli V., Lagrange N., 
Mencattini A., Parrini M. C., Martinelli E., Isambert H.,
[CausalXtract: a flexible pipeline to extract causal effects from live-cell 
time-lapse imaging data, eLife 2024](https://www.biorxiv.org/content/10.1101/2024.02.06.579177v1.abstract).

Ribeiro-Dantas M. D. C., Li H., Cabeli V., Dupuis L., Simon F., Hettal L., 
Hamy A. S., Isambert H.,
[Learning interpretable causal networks from very large datasets, application 
to 400,000 medical records of breast cancer patients, iScience, 2024](https://doi.org/10.1016/j.isci.2024.109736).

Cabeli V., Li H., Ribeiro-Dantas M., Simon F., Isambert H.,
[Reliable causal discovery based on mutual information supremum principle for 
finite dataset, Why21 at NeurIPS 2021](https://why21.causalai.net/papers/WHY21_24.pdf).

Cabeli V., Verny L., Sella N., Uguzzoni G., Verny M., Isambert H.,
[Learning clinical networks from medical records based on information estimates
in mixed-type data, PLoS Comput. Biol.  2020](https://doi.org/10.1371/journal.pcbi.1007866)
| [code](https://github.com/vcabeli/miic_PLoS)

Li H., Cabeli V., Sella N., Isambert H.,
[Constraint-based causal structure learning with consistent separating sets, 
In Advances in Neural Information Processing Systems 2019.](https://papers.nips.cc/paper/9573-constraint-based-causal-structure-learning-with-consistent-separating-sets)
| [code](https://github.com/honghaoli42/consistent_pcalg).

Verny L., Sella N., Affeldt S., Singh PP., Isambert H.,
[Learning causal networks with latent variables from multivariate information 
in genomic data, PLoS Comput. Biol. 2017](https://doi.org/10.1371/journal.pcbi.1005662).

Affeldt S., Isambert H.,
[Robust Reconstruction of Causal Graphical Models based on Conditional 2-point 
and 3-point Information, UAI 2015](https://auai.org/uai2015/proceedings.shtml)

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

MIIC allows you to create a graph object from a dataset of observations 
of both discrete and continuous variables, potentially with missing values 
and taking into account unobserved latent variables.
You can find this example along others by calling the documentation 
of the main function `?miic` from R.

```R
library(miic)

# EXAMPLE HEMATOPOIESIS
data(hematoData)
# execute MIIC (reconstruct graph)
miic_obj <- miic(
  input_data = hematoData, latent = "yes",
  n_shuffles = 10, conf_threshold = 0.001
)

# plot graph with igraph
if(require(igraph)) {
  plot(miic_obj, method="igraph")
}
```

## Documentation

You can find the documentation pages in the "man" folder, in the auto generated
[PDF](https://CRAN.R-project.org/package=miic/miic.pdf), 
or use R functions `help()` and `?`.
  
## Authors

- Tiziana Tocci
- Nikita Lagrange
- Orianne Debeaupuis
- Louise Dupuis
- Franck Simon
- Vincent Cabeli
- Honghao Li
- Marcel Ribeiro Dantas
- Verny Louis
- Sella Nadir
- Séverine Affeldt
- Hervé Isambert

## License

GPL-2 | GPL-3
