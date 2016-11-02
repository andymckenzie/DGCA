# DGCA

The goal of DGCA is to calculate differential correlations across conditions.

It simplifies the process of seeing whether two correlations are different without having to rely solely on parametric assumptions by leveraging non-parametric permutation tests.

It also has several other options including calculating the average differential correlation between groups of genes, gene ontology enrichment analyses of the results, and differential correlation network identification via integration with MEGENA.  

## Installation

You can install DGCA from github with:

```R
# install.packages("devtools")
devtools::install_github("andymckenzie/DGCA")
```

## Basic Example

```R
library(DGCA)
data(darmanis); data(design_mat)
res = ddcorAll()
ddcor_res = ddcorAll(inputMat = darmanis, design = design_mat, compare = c("oligodendrocyte", "neuron"))
head(ddcor_res)
```

## More

There are three vignettes available in order to help you learn how to use the package:

- DGCA_basic: This will get you up-and-going quickly.
- DGCA: This is a more extended version that explains a bit about how the package works and shows several of the options available in the package.
- DGCA_modules: This will show you how to use the package to perform module-based and network-based analyses.

The second two vignettes can be found in inst/doc.
