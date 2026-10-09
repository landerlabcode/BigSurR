# MergeSmallModules

MergeSmallModules

## Usage

``` r
MergeSmallModules(corr.matrix, modules, min.size = 15)
```

## Arguments

- corr.matrix:

  Sparse matrix containing the correlations of interest.

- modules:

  igraph communities object produced by FindCorrelationModules.

- min.size:

  Integer. Minimum size of output modules (if possible).

## Value

igraph communities object containing the merged gene modules.
