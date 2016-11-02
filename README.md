## Why should I use it?

DGCA is an R package designed to calculate differential correlations across conditions.

It simplifies the process of seeing whether two correlations are different without having to rely solely on parametric assumptions by leveraging non-parametric permutation tests.

It also has several other options including calculating the average differential correlation between groups of genes, gene ontology enrichment analyses of the results, and differential correlation network identification via integration with MEGENA.  

## How do I use it?

There are three vignettes available in order to help you learn how to use the package, one with basic functions, one that goes more in-depth, and one about how to use the functions designed to help you with module-based analyses.

## How do I get it?

You can install via CRAN via install.packages("DGCA").

You can also install the development version via

install.packages("devtools")
library(devtools)
install_github("andymckenzie/DGCA")
