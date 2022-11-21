
# fastupgma

<!-- badges: start -->
<!-- badges: end -->

The goal of fastupgma is to ...

## Installation

You can install the development version of fastupgma from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("alrobles/fastupgma")
```

## Example

This is an r package that wraps the upgma function as hclust function 
with average method. We then wrap the hclust object as an ape phylo object
to create a phylogenetic tree.

We preload two distances matrices filtered from original consensus phylogenetic tree in [Upham et al 2019](https://doi.org/10.1371/journal.pbio.3000494). 
One with the Crocidura genus and another with Myotis genus.

``` r
library(fastupgma)

## basic example code

data(CrociduraDistance)
data(rodentiaDistance)
```

Then simply use the fastupgma function to create the tree
```r
treeOutput = fastupgma(CrociduraDistance)
```

This is a phylo object from ape and you can manipulated for further analysis.

