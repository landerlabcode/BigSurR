get.residuals <- function(seuratob, counts.slot, c = F){

  if((c!=F) & (c<=0)){
    stop("c must be greater than zero.")
  }

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
  ematrix <- gene.totals %o% depthlist

  if(c==F){
    c <- find.best.c(ematrix, depthlist, genetotals, num.cells, raw.counts)
  }

  residuals <- (raw.counts-ematrix)/sqrt(ematrix*(1+c^2*ematrix))

  data.list <- list(residuals, ematrix, depthlist, num.cells, gene.totals, c)
  data.list <- setNames(data.list, c("residuals","ematrix","depthlist","num.cells", "gene.totals", "c"))
  return(data.list)
}


find.best.c <- function(ematrix, depthlist, genetotals, num.cells, raw.counts){
  test.cs <- seq(from= 0, to= 1, by = 0.05)

  best.slope <- c("index"=0,"slope"=2000)

  for(i in 1:length(test.cs)){

    c <- test.cs[i]

    X <- log10(rowMeans(raw.counts))

    Y <- log10(1/(num.cells-1) * rowSums(((raw.counts-ematrix)/sqrt(ematrix*(1+c^2*ematrix)))^2))

    fit.indexes <- which((-1 < X) & (X < 2))

    X <- X[fit.indexes]
    Y <- Y[fit.indexes]

    test.fit <- lm(Y~X)

    coefficients <- summary(test.fit)$coefficients

    slope <- coefficients[2,1]

    if(abs(slope-0) < abs(best.slope[2]-0)){
      best.slope <- c("index"=i, "slope"=slope)
    }

    if(slope < 0){
      break
    }
  }


  return(test.cs[best.slope[1]])
}
