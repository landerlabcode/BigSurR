depth.scaling <- function(seuratob,  assay = "RNA", counts.slot="counts"){

  DefaultAssay(object = seuratob) <- assay
  raw.counts <- as.matrix(GetAssayData(object=seuratob, slot=counts.slot))

  if(!(all(raw.counts==floor(raw.counts)))){
    stop("Non-integer counts data detected. Please supply unnormalized data.")
  }
  gene.totals <- rowSums(raw.counts)

  if(0 %in% gene.totals){
    stop("Genes with zero counts were found, run quality control steps before attempting further analysis.")
  }

  cell.total.umis <- colSums(raw.counts)
  all.umis<- sum(cell.total.umis)

  num.cells <- ncol(raw.counts)
  num.genes <- nrow(raw.counts)

  depthlist <- cell.total.umis/all.umis

  return(depthlist)
}

