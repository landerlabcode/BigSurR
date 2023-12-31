% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BigSur.R
\name{BigSur}
\alias{BigSur}
\title{BigSur (Basic Informatics and Gene Statistics from Unnormalized Reads)}
\usage{
BigSur(
  seurat.obj,
  assay = "RNA",
  counts.slot = "counts",
  c = F,
  variable.features = T,
  correlations = F,
  first.pass.cutoff = 2,
  inverse.fano.moments = T,
  fano.alpha = 0.05,
  min.fano = 1.5,
  cor.alpha = 0.05,
  log.file = T,
  log.file.dir = paste0(getwd(), "/BigSurRun", Sys.Date(), ".txt")
)
}
\arguments{
\item{seurat.obj}{Seurat object containing the raw transcript counts filtered for zero count genes.}

\item{assay}{Assay slot containing raw transcript counts (default "RNA").}

\item{counts.slot}{Slot within assay containing raw counts matrix (default "counts").}

\item{c}{Boolean, double. The coefficient of variation of the dataset. If set to false, one will be estimated using the dataset.}

\item{variable.features}{Boolean. If true, BigSur will identify select variable features based on the modified corrected Fano factor.}

\item{correlations}{Boolean. If true, BigSur will identify statistically significant gene-gene correlations.}

\item{first.pass.cutoff}{Integer. Removes roots before p-value calculations if the root is below Abs[Sqrt(2)*InverseErfc(2*10^-first.pass.cutoff)]. The higher the number, the more correlations are removed in initial screening.}

\item{inverse.fano.moments}{Boolean. If true, BigSur will calculate the moments for the inverse Fano factor pairs before performing Cornish Fisher expansion.}

\item{fano.alpha}{Double. Desired false discovery cutoff for labeling of variable features. (Default 0.05).}

\item{min.fano}{Double. Minimum mcFano value considered for variable genes.}

\item{cor.alpha}{Double. Desired false discovery cutoff for labeling of statistically significant correlations.}

\item{log.file}{Boolean. If true, a log file will be created.}

\item{log.file.dir}{String. Path of desired location for log file.}
}
\value{
If both variable features and correlations are identified, a list containing the updated Seurat object and the
statistically significant correlations is returned. If only one process is selected, their respective output is returned alone.
}
\description{
BigSur (Basic Informatics and Gene Statistics from Unnormalized Reads)
}
\examples{
BigSur(example.seurat, variable.features=T, correlations=T)




}
