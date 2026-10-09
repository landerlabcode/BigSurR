# BigSur (Basic Informatics and Gene Statistics from Unnormalized Reads)

BigSur (Basic Informatics and Gene Statistics from Unnormalized Reads)

## Usage

``` r
BigSur(
  seurat.obj,
  assay = "RNA",
  counts.slot = "counts",
  cv.est.method = "MeanSpecific",
  null.distribution = "NB",
  variable.features = T,
  correlations = F,
  first.pass.cutoff = 2,
  inverse.fano.moments = T,
  fano.alpha = 0.05,
  min.fano = 1.5,
  cor.alpha = 0.05,
  block.size = 500,
  log.file = F,
  log.file.dir = paste0(getwd(), "/BigSurRun", Sys.Date(), ".txt")
)
```

## Arguments

- seurat.obj:

  Seurat object containing the raw transcript counts filtered for zero
  count genes.

- assay:

  Assay slot containing raw transcript counts (default "RNA").

- counts.slot:

  Slot within assay containing raw counts matrix (default "counts").

- cv.est.method:

  String. Sets the method used for determining the coefficient of
  variation used to define null distributions. There are three options
  here: 1) "Single": Determines a scalar value of c. 2) "MeanSpecific":
  Estimates the relationship between mean expression and the expected
  coefficient of variation and predicts a null value for each gene.
  3)"TwoComponent": Estimates a the relationship between mean expression
  and the expected coefficient of variation using a two parameter fit.

- null.distribution:

  String. Sets the null distribution from which to estimate p-values for
  Fano factors and correlations. "NB": Negative binomial, "PLN": Poisson
  log-normal. (Defaults to negative binomial.)

- variable.features:

  Boolean. If true, BigSur will identify select variable features based
  on the modified corrected Fano factor.

- correlations:

  Boolean. If true, BigSur will identify statistically significant
  gene-gene correlations.

- first.pass.cutoff:

  Integer. Removes roots before p-value calculations if the root is
  below Abs\[Sqrt(2)\*InverseErfc(2\*10^-first.pass.cutoff)\]. The
  higher the number, the more correlations are removed in initial
  screening.

- inverse.fano.moments:

  Boolean. If true, BigSur will calculate the moments for the inverse
  Fano factor pairs before performing Cornish Fisher expansion.

- fano.alpha:

  Double. Desired false discovery cutoff for labeling of variable
  features. (Default 0.05).

- min.fano:

  Double. Minimum mcFano value considered for variable genes.

- cor.alpha:

  Double. Desired false discovery cutoff for labeling of statistically
  significant correlations.

- block.size:

  Integer. Determines the block size for correlation cumulant
  calculation to relieve memory pressure at the cost of some speed.
  Number indicates the number of gene pairs evaluated at a time.
  (Default 500)

- log.file:

  Boolean. If true, a log file will be created.

- log.file.dir:

  String. Path of desired location for log file.

## Value

If both variable features and correlations are identified, a list
containing the updated Seurat object and the statistically significant
correlations is returned. If only one process is selected, their
respective output is returned alone.

## Examples

``` r
if (FALSE) { # \dontrun{
out <- BigSur(example.seurat, variable.features = TRUE, correlations = TRUE)
} # }



```
