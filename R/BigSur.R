BigSur <- function(seurat.obj,
                   assay = "RNA",
                   counts.slot="counts",
                   c=F,
                   variable.features=T,
                   correlations=F,
                   first.pass.cutoff=2,
                   inverse.fano.moments = T,
                   fano.alpha = 0.05,
                   min.fano = 1.5,
                   cor.alpha = 0.05,
                   log.file = T,
                   log.file.dir = paste0(getwd(), "/BigSurRun", Sys.Date(),".txt")
                   )
  {
  if((variable.features!=T) & (correlations!=T)){
    stop("Both variable.features and correlations are set to false.")
  }

  if(log.file==T){
    lf <- log_open(log.file.dir)
    log_print("Pipeline started execution.")
  }

  residuals<-get.residuals(seurat.obj, counts.slot, c)

  c <- residuals$c

  num.genes <- length(residuals$gene.totals)

  if(log.file==T){
    log_print("Modified corrected Pearson residuals calculated.")
  }
  mcfanos <- get.mcFanos(residuals)


  if(log.file==T){
    log_print("Modified corrected Fano factors calculated.")
  }

  if(variable.features==T){
    if(log.file==T){
      log_print("Beginning identification of significant mcFanos.")
    }
    fano.cumulants <- Cumulants.Fano(residuals, c)
    fanocoeffs <- CF.Coefficients.Fano(fano.cumulants[,2], fano.cumulants[,3], fano.cumulants[,4], fano.cumulants[,5], mcfanos, rownames(residuals$ematrix))
    fanoroots <- CF.AllRoots(fanocoeffs)
    pval <- CF.pval(fanoroots)
    p.df <- data.frame(mcfanos, pval)
    p.df$roots <- fanoroots
    fanoBH <- Fano.BH(p.df, num.genes)
    fano.selected <- Fano.HighlyVariable(fanoBH, fano.alpha, min.fano)
    top.features <- row.names(fano.selected[fano.selected[,5]==T,])
    new.seurat.obj <- seurat.obj
    VariableFeatures(object = new.seurat.obj) <- top.features
    feat.metadata <- as.data.frame(fano.selected[,c(1,4,5)])
    new.seurat.obj[[assay]]@meta.features <- new.seurat.obj[[assay]]@meta.features[match(row.names(feat.metadata), row.names(new.seurat.obj[[assay]]@meta.features)),]
    new.seurat.obj[[assay]]@meta.features <- cbind(new.seurat.obj[[assay]]@meta.features, feat.metadata)
    if(log.file==T){
      log_print("Highly variable features identified.")
    }
  }

  if(correlations==T){
    if(log.file==T){
      log_print("Beginning correlation calculation.")
    }
    pcc <- get.mcPCC2(residuals, mcfanos)
    if(log.file==T){
      log_print("Modified-corrected Pearson Correlation Coefficients calculated.")
    }

    if(inverse.fano.moments==T){
      inv.correction <- inv.sqrt.correction2(residuals, residuals$c)
      a <- max(2, min(residuals$gene.totals))
      e <- residuals$num.cells/50
      h <- max(residuals$gene.totals)
      points <- as.integer(c(a, a*(e/a)^(1/4), a*(e/a)^(1/2), a*(e/a)^(3/4), e, e*(h/e)^(1/3), e*(h/e)^(2/3), h))
      moment.interp <- inv.sqrt.moment.interpolation2(inv.correction, residuals$gene.totals, points)
      if(log.file==T){
        log_print("Inverse sqrt moments calculated.")
      }
    }

    else{
      onesmat <- matrix(1, nrow=num.genes, ncol=num.genes)
      moment.interp <- list(onesmat, onesmat, onesmat, onesmat)
    }

    cor.cumulants <- Cumulants.PCC(residuals, moment.interp)
    if(log.file==T){
      log_print("PCC cumulants calculated.")
    }

    cor.coefficients <- CF.Coefficients.PCC(cor.cumulants, pcc)
    if(log.file==T){
      log_print("PCC Cornish Fisher coefficients calculated.")
    }

    cor.roots <- CF.PCC.Roots(cor.coefficients, 2, residuals$gene.totals)

    cor.p <- CF.PCC.pval(cor.roots)
    if(log.file==T){
      log_print("P-values calculated.")
    }

    cor.signmat <- get.signMat(pcc)
    if(log.file==T){
      log_print("Sign matrix calculated.")
    }

    equivalent.pccs <- get.inferred.PCCs.2(cor.p, cor.signmat, residuals$num.cells, num.genes)
    if(log.file==T){
      log_print("Equivalent PCCs calculated in.")
    }
    sig.equivalent.pccs <- get.significant.inferred.PCCs(cor.p, equivalent.pccs, num.genes, alpha)

    if(log.file==T){
      log_print("Equivalent PCCs filtered for significants.")
      log_print(paste0("Number of remaining correlations:", Matrix::nnzero(sig.equivalent.pccs)))
    }
  }


  if(variable.features==T & correlations==T){
    return(list(new.seurat.obj, sig.equivalent.pccs))}
  else if(variable.features==T){
    return(new.seurat.obj)
  }
  else{
    return(sig.equivalent.pccs)
  }
  if(log.file==T){
    log_close()
  }

}
