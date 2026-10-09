# FindCorrelationModules

FindCorrelationModules

## Usage

``` r
FindCorrelationModules(corr.matrix, num.steps = 4, weighted.clustering = FALSE)
```

## Arguments

- corr.matrix:

  Sparse matrix containing the correlations of interest.

- num.steps:

  Integer. Number of steps used for walktrap clustering. (Default = 4)

- weighted.clustering:

  Boolean. Determines whether or not the walktrap will be performed on
  the correlations matrix or the adjacency matrix. (Default = FALSE)

## Value

igraph communities object.
