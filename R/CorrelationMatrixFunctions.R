SubsetCorrelationMatrix <- function(
    seuratObject,
    pCutoff=F,
    minPosCorr=F,
    minNegCorr=F,
    genes=F
){


  if(!(is.character(genes) || identical(genes, FALSE))){
    stop("Genes were not supplied as a vector of strings.")
  }
  if(!(is.numeric(pCutoff) || identical(pCutoff, FALSE))){
    stop("p-value cutoff not supplied as a number.")
  }
  if(!(is.numeric(minPosCorr) || identical(minPosCorr, FALSE))){
    stop("Minimum positive correlation cutoff not supplied as a number.")
  }
  if(!(is.numeric(minNegCorr) || identical(minNegCorr, FALSE))){
    stop("Minimum negative correlation cutoff not supplied as a number.")
  }

  corrs <- seuratObject@misc$BigSur.Correlations
  if (is.null(corrs)){
    stop("No BigSur correlations in this object.")
  }
  ps <- seuratObject@misc$BigSur.log.adj.pvalues
  if (is.numeric(pCutoff) && is.null(ps)){
    stop("No BigSur p-values in this object.")
  }

  if (is.numeric(pCutoff)) {
    stopifnot(identical(corrs@i, ps@i), identical(corrs@p, ps@p))
    corrs@x[ps@x > log(pCutoff)] <- 0
    corrs <- Matrix::drop0(corrs)
  }

  if (is.numeric(minPosCorr)) {
    corrs@x[corrs@x > 0 & corrs@x < minPosCorr] <- 0
    corrs <- Matrix::drop0(corrs)
  }

  if (is.numeric(minNegCorr)) {
    corrs@x[corrs@x < 0 & abs(corrs@x) < abs(minNegCorr)] <- 0
    corrs <- Matrix::drop0(corrs)
  }

  if(is.character(genes)){
    missing <-setdiff(genes, rownames(corrs))
    if(length(missing)){
      warning("Some requested genes not present in matrix: ", paste(missing, collapse=", "))
    }
    p.genes <- intersect(genes, rownames(corrs))
    corrs <- corrs[p.genes, p.genes, drop=F]
  }
  return(corrs)
}
