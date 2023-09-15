# BigSurR (Basic Informatics and Gene Statistics from Unnormalized Reads R)

BigSur is a tool for single cell RNA sequencing analysis which interfaces with Seurat objects to both
- Select statistically significant features to use in clustering.
- Identify statistically significant correlations between genes.

The principles surrounding these analyses are described in ["Leveraging gene correlations in single cell transcriptomic data"][1].

## Installation
BigSurR can be installed directly from github using the devtools package.
```{r}
devtools::install_github("landerlabcode/BigSurR")
```

## Usage
To use BigSur, first create a Seurat object containing the raw transcript counts for each gene. As a preprocessing step, ensure that there are no genes which have zero counts across all cells. **Ensure that the active assay is the one containing the raw transcript counts.**
```{r}
library(BigSur)
example.seurat <- CreateSeuratObject(data, min.cells=1)
DefaultAssay(object = example.seurat) <- "RNA"
```
Pass this object into the BigSur function with the desired parameters. The parameters are as follows:
- **seurat.obj**: Seurat object containing the raw transcript counts filtered for zero count genes.
- **assay**: Assay slot containing raw transcript counts (default "RNA").
- **counts.slot**: Slot within assay containing raw counts matrix (default "counts").
- **c**: Boolean, double. The coefficient of variation of the dataset. If set to false, one will be estimated using the dataset.
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
If both variable features and correlations are identified, a list containing the updated Seurat object and the an adjacency matrix of the statistically significant correlations is returned. If only one process is selected, their respective output is returned alone.

To change:
- Add documentation for help function
- Add example usage
- Add more assertions to prevent improper use
- Change import structure to remove warnings from unused function conflicts
- Change minimum Fano factor from a set number to a quantile-based one.
- Add example correlation visualization using igraph
- Add custom igraph functions for visualization to package

In progress:
- Testing installation and code on Windows and Linux systems.

[1]: https://www.biorxiv.org/content/10.1101/2023.03.14.532643v1 "Leveraging gene correlations in single cell transcriptomic data"
