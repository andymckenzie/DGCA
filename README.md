[![Travis-CI Build Status](https://travis-ci.org/andymckenzie/DGCA.svg?branch=master)](https://travis-ci.org/andymckenzie/DGCA)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/ggplot2)](https://cran.r-project.org/package=DGCA)

# DGCA

The goal of DGCA is to calculate differential correlations across conditions.

It simplifies the process of seeing whether two correlations are different without having to rely solely on parametric assumptions by leveraging non-parametric permutation tests and adjusting the resulting empirical p-values for multiple corrections using the qvalue R package.

It also has several other options including calculating the average differential correlation between groups of genes, gene ontology enrichment analyses of the results, and differential correlation network identification via integration with MEGENA.  

## Installation

You can install DGCA from CRAN with:

```R
install.packages("DGCA")
```

You can install the development version of DGCA from github with:

```R
# install.packages("devtools")
devtools::install_github("andymckenzie/DGCA")
```

## Basic Example

```R
library(DGCA)
data(darmanis); data(design_mat)
ddcor_res = ddcorAll(inputMat = darmanis, design = design_mat, compare = c("oligodendrocyte", "neuron"))
head(ddcor_res, 3)
#   Gene1  Gene2 oligodendrocyte_cor oligodendrocyte_pVal neuron_cor neuron_pVal
# 1 CACYBP   NACA        -0.070261455           0.67509118  0.9567267           0
# 2 CACYBP    SSB        -0.055290516           0.74162636  0.9578999           0
# 3 NDUFB9    SSB        -0.009668455           0.95405875  0.9491904           0
#   zScoreDiff     pValDiff     empPVals pValDiff_adj Classes
# 1  10.256977 1.100991e-24 1.040991e-05    0.6404514     0/+
# 2  10.251847 1.161031e-24 1.040991e-05    0.6404514     0/+
# 3   9.515191 1.813802e-21 2.265685e-05    0.6404514     0/+
```

## Vignettes

There are three vignettes available in order to help you learn how to use the package:

- DGCA_basic: This will get you up-and-going quickly.
- DGCA: This is a more extended version that explains a bit about how the package works and shows several of the options available in the package.
- DGCA_modules: This will show you how to use the package to perform module-based and network-based analyses.

The second two vignettes can be found in inst/doc.

## Applications

You can view the manuscript describing DGCA in detail as well as several applications here:

- http://bmcsystbiol.biomedcentral.com/articles/10.1186/s12918-016-0349-1

Material for associated simulations and networks created from MEGENA can be found here:

- https://github.com/andymckenzie/dgca_manuscript
