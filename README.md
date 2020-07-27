# MIIC
  <!-- badges: start -->
  [![R build status](https://github.com/miicTeam/miic_R_package/workflows/R-CMD-check/badge.svg)](https://github.com/miicTeam/miic_R_package/actions)
  <!-- badges: end -->

This repository contains the source code for MIIC, a method which learns a large class of causal or non-causal graphical models from purely observational data, while including the effects of unobserved latent variables, commonly found in many datasets. Starting from a complete graph, the method iteratively removes dispensable edges, by uncovering significant information contributions from indirect paths, and assesses edge-specific confidences from randomization of available data. The remaining edges are then oriented based on the signature of causality in observational data. This approach can be applied on a wide range of datasets and provide new biological insights on regulatory networks from single cell expression data, genomic alterations during tumor development and co-evolving residues in protein structures. For more information you can refer to: Cabeli et al. PLoS Comp. Bio. 2020 <doi:10.1371/journal.pcbi.1007866>, Verny et al. PLoS Comp. Bio. 2017 <doi:10.1371/journal.pcbi.1005662>.

## References
Verny L., Sella N., Affeldt S., Singh PP., Isambert H.; Learning causal networks with latent variables from multivariate information in genomic data;  PLoS Comput. Biol., 2017.

## Prerequisites
MIIC contains R and C++ sources. In order to compile MIIC you will need a c++ compiler with c++ 14.
The following R packages are dependencies : MASS, igraph, ppcor, scales, stats, Rcpp

## Installing

You may install the latest stable version that was submitted to CRAN directly from R using the package manager :
```{r}
install.packages("miic")
```
Or if you wish the latest version under development, you may clone this repository and install the package from the source using the console commmand :
```
cd miic_R_package
R CMD INSTALL .
```
## Documentation
You can find the documentation pages in the "man" folder.

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
