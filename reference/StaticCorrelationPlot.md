# StaticCorrelationPlot

StaticCorrelationPlot

## Usage

``` r
StaticCorrelationPlot(corr.matrix, highlight = NULL, modules = NULL)
```

## Arguments

- corr.matrix:

  Sparse matrix containing the correlations of interest.

- highlight:

  Vector of gene symbols to be highlighted. If supplied, the relevant
  gene symbols will be written in red.

- modules:

  igraph communities object. If supplied, the nodes will change color
  and shape to indicate the module they belong to.

## Value

ggraph plot of correlations.
