# SubsetCorrelationMatrix

SubsetCorrelationMatrix

## Usage

``` r
SubsetCorrelationMatrix(
  seuratObject,
  pCutoff = F,
  minPosCorr = F,
  minNegCorr = F,
  genes = F
)
```

## Arguments

- seuratObject:

  Seurat object containing correlations matrix produced by BigSur.

- pCutoff:

  Boolean, Float. p-value threshold for correlations. Correlations with
  p-values larger than pCutoff will be set to 0.

- minPosCorr:

  Boolean, Float. Minimum positive correlation cutoff. Positive
  correlations with coefficients under minPosCorr will be set to 0.

- minNegCorr:

  Boolean, Float. Minimum negative correlation cutoff. Negative
  correlations with absolute values under abs(minNegCorr) will be set to
  0.

- genes:

  Boolean, Character Vector. The set of genes for the correlation matrix
  to be subset to.

## Value

Sparse matrix (dsCMatrix) containing correlations meeting subset
criteria.
