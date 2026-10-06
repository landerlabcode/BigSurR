# BigSur 2.0 (Basic Informatics and Gene Statistics from Unnormalized Reads)

BigSur is a tool for single cell RNA sequencing analysis which interfaces with Seurat objects to both
- Select statistically significant features to use in clustering.
- Identify statistically significant correlations between genes.

The principles surrounding these analyses are described in ["Leveraging gene correlations in single cell transcriptomic data"][1] and ["Statistically principled feature selection for single cell transcriptomics"][2].

## Installation
BigSurR can be installed directly from github using the devtools package.
```{r}
devtools::install_github("landerlabcode/BigSurR")
```

Note, Seurat recently changed the structure of their assay objects. BigSur's feature selection process will no longer work on versions below 5.1.0 (Assay5 class required).

## Usage

### General
To use BigSur, first create a Seurat object containing the raw transcript counts for each gene. As a preprocessing step, ensure that there are no genes which have zero counts across all cells.
```{r}
library(BigSur)
example.seurat <- CreateSeuratObject(data, min.cells=1)
```
Pass this object into the BigSur function with the desired parameters. The parameters are as follows:
- **seurat.obj**: Seurat object containing the raw transcript counts filtered for zero count genes.
- **assay**: Assay slot containing raw transcript counts (default "RNA").
- **counts.slot**: Slot within assay containing raw counts matrix (default "counts").
- **cv.est.method**: String. Sets the method used for determining the coefficient of variation used to define null distributions. There are three options here: 1) "Single": Determines a scalar value of c. 2) "MeanSpecific": Estimates the relationship between mean expression and the expected coefficient of variation and predicts a null value for each gene. 3)"TwoComponent": Estimates a the relationship between mean expression and the expected coefficient of variation using a two parameter fit.
- **null.distribution**: String. Sets the null distribution from which to estimate p-values for Fano factors and correlations. "NB": Negative binomial, "PLN": Poisson log-normal. (Defaults to negative binomial.)
- **variable.features**: Boolean. If true, BigSur will identify select variable features based on the modified corrected Fano factor. (Default true)
- **correlations**:Boolean. If true, BigSur will identify statistically significant gene-gene correlations. (Default false)
- **first.pass.cutoff**: Integer. Removes roots before p-value calculations if the root is below Abs[Sqrt(2)\*InverseErfc(2*10^-first.pass.cutoff)]. The higher the number, the more correlations are removed in initial screening.
- **inverse.fano.moments**: Boolean. If true, BigSur will calculate the moments for the inverse Fano factor pairs before performing Cornish Fisher expansion. Setting false will significantly speed up calculation at the cost of some accuracy.
- **fano.alpha**: Double. Desired false discovery cutoff for labeling of variable features. (Default 0.05).
- **min.fano**: Double. Minimum mcFano value considered for variable genes.
- **cor.alpha**: Double. Desired false discovery cutoff for labeling of statistically significant correlations.
- **log.file**: Boolean. If true, a log file will be created.
- **log.file.dir**: String. Path of desired location for log file. *Note: Default string is set to work on Unix based file structures (i.e., manually set this on Windows).

For example, to calculate both the highly variable features and the significant correlations in your dataset with default false discovery cutoffs you could use:
```{r}
example.output <- BigSur(example.seurat, variable.feature=T, correlations=T)
```
This function will return a new Seurat object with the following attachments:
If variable.features == True:
  -Modified corrected Fano factors (mcfanos) and their respective adjusted (Benjamini Hochberg) p values are stored in the assay's @meta.data slot. The highly variable features boolean vector is stored in the typical slot used by FindVariableFeatures.
If correlations == True:
  -Modified corrected Pearson correlation coefficients and their respective adjusted (Benjamini Hochberg) p values (log10 format) are stored in the @misc slot.
Regardless of parameters:
  -The eta and theta parameter values used to estimate the per-gene value of the null coefficient of variation are stored in the @misc slot.
  -The modified corrected Pearson residuals matrix is stored in the assay's $data slot.


## 2.0 Changes
- Added support for the Negative Binomial null distribution.
- Empirically estimates a per-gene value of the null coefficient of variation rather than using a single value.
- Optimizes memory efficiency by minimizing the number of large matrices being stored in intermediate steps.
- Outputs are now completely integrated into the Seurat object.
- Added basic visualization and correlation matrix manipulation functions.

## Future updates
- Estimating the null inverse square root Fano factor moments is a stochastic process leading to small differences in the number of significant correlations in each run. We would like to minimize this.

[1]: https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05926-z "Leveraging gene correlations in single cell transcriptomic data"
[2]: https://link.springer.com/article/10.1186/s12859-025-06240-y "Statistically principled feature selection for single cell transcriptomics"
