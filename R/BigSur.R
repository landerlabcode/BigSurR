#' BigSur (Basic Informatics and Gene Statistics from Unnormalized Reads)
#'
#' @param seurat.obj Seurat object containing the raw transcript counts filtered for zero count genes.
#' @param assay Assay slot containing raw transcript counts (default "RNA").
#' @param counts.slot Slot within assay containing raw counts matrix (default "counts").
#' @param c Boolean, double. The coefficient of variation of the dataset. If set to false, one will be estimated using the dataset.
#' @param variable.features Boolean. If true, BigSur will identify select variable features based on the modified corrected Fano factor.
#' @param correlations Boolean. If true, BigSur will identify statistically significant gene-gene correlations.
#' @param first.pass.cutoff Integer. Removes roots before p-value calculations if the root is below Abs[Sqrt(2)*InverseErfc(2*10^-first.pass.cutoff)]. The higher the number, the more correlations are removed in initial screening.
#' @param inverse.fano.moments Boolean. If true, BigSur will calculate the moments for the inverse Fano factor pairs before performing Cornish Fisher expansion.
#' @param fano.alpha Double. Desired false discovery cutoff for labeling of variable features. (Default 0.05).
#' @param min.fano Double. Minimum mcFano value considered for variable genes.
#' @param depthlist Boolean, vector. If a vector is supplied, that vector of values will be used to scale counts to account for unequal sequencing depth across cells. If left as False, this scaling will be calculated during the BigSur run.
#' @param cor.alpha Double. Desired false discovery cutoff for labeling of statistically significant correlations.
#' @param return.ps Boolean. If true, the Benjamini-Hochberg corrected p-values associated with each equivalent PCC will be returned in a list with the equivalent PCC sparse matrix. The first object in this list will be the equivalent PCCs, the second will be the p-value matrix.
#' @param log.file Boolean. If true, a log file will be created.
#' @param log.file.dir String. Path of desired location for log file.
#'
#' @return If both variable features and correlations are identified, a list containing the updated Seurat object and the
#' statistically significant correlations is returned. If only one process is selected, their respective output is returned alone.
#' @export
#'
#' @examples BigSur(example.seurat, variable.features=T, correlations=T)
#'
#'
#'
#'
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
                   depthlist = F,
                   return.ps = F,
                   log.file = T,
                   log.file.dir = paste0(getwd(), "/BigSurRun", Sys.Date(),".txt")
                   )
  {
  if((variable.features!=T) & (correlations!=T)){
    stop("Both variable.features and correlations are set to false.")
  }

  if(packageVersion("Seurat") < "5.0.1"){
    stop("Older versions of Seurat still utilize 'meta.features' which has been replaced by 'meta.data' in newer versions. Please upgrade to 5.0.1 at minimum.")
  }
  print("Pipeline started execution.")
  if(log.file==T){
    fileConn <- file(log.file.dir, open ="wt")
    write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": Pipeline started execution."), file=fileConn, append=T)
  }

  residuals<-get.residuals(seurat.obj, assay, counts.slot, c, depthlist)

  c <- residuals$c

  num.genes <- residuals$num.genes
  print("Modified corrected Pearson residuals calculated.")
  if(log.file==T){
    write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": Modified corrected Pearson residuals calculated."), file=fileConn, append=T)
    }
  mcfanos <- get.mcFanos(residuals)

  print("Modified corrected Fano factors calculated.")
  if(log.file==T){
    write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": Modified corrected Fano factors calculated."), file=fileConn, append=T)

    }


  if(variable.features==T){
    print("Beginning identification of significant mcFanos.")
    if(log.file==T){
      write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": Beginning identification of significant mcFanos."), file=fileConn, append=T)
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
    feat.metadata <- as.data.frame(fano.selected[,c(1,4,5)])
    new.seurat.obj[[assay]]@meta.data <- feat.metadata
    VariableFeatures(new.seurat.obj) <- top.features
    new.seurat.obj[[assay]]$data <- residuals$residuals
    print("Highly variable features identified.")
    if(log.file==T){
      write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": Highly variable features identified."), file=fileConn, append=T)

      }
  }

  if(correlations==T){
    if(log.file==T){
      print("Beginning correlation calculation.")
      write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": Beginning correlation calculation."), file=fileConn, append=T)

       }
    pcc <- get.mcPCC2(residuals, mcfanos)
    print("Modified-corrected Pearson Correlation Coefficients calculated.")
    if(log.file==T){
      write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": Modified-corrected Pearson Correlation Coefficients calculated."), file=fileConn, append=T)

      }

    if(inverse.fano.moments==T){
      inv.correction <- inv.sqrt.correction2(residuals, residuals$c)
      a <- max(2, min(residuals$gene.totals))
      e <- residuals$num.cells/50
      h <- max(residuals$gene.totals)
      points <- as.integer(c(a, a*(e/a)^(1/4), a*(e/a)^(1/2), a*(e/a)^(3/4), e, e*(h/e)^(1/3), e*(h/e)^(2/3), h))
      moment.interp <- inv.sqrt.moment.interpolation2(inv.correction, residuals$gene.totals, points)
      print("Inverse sqrt moments calculated.")
      if(log.file==T){
        write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": Inverse sqrt moments calculated."), file=fileConn, append=T)
        }
    }

    else{
      onesmat <- matrix(1, nrow=num.genes, ncol=num.genes)
      moment.interp <- list(onesmat, onesmat, onesmat, onesmat)
    }

    cor.cumulants <- Cumulants.PCC(residuals, moment.interp)
    print("PCC cumulants calculated.")
    if(log.file==T){
      write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": PCC cumulants calculated."), file=fileConn, append=T)

      }

    cor.coefficients <- CF.Coefficients.PCC(cor.cumulants, pcc)
    print("PCC Cornish Fisher coefficients calculated.")
    if(log.file==T){
      write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": PCC Cornish Fisher coefficients calculated."), file=fileConn, append=T)

       }

    cor.roots <- CF.PCC.Roots(cor.coefficients, 2, residuals$gene.totals)

    cor.p <- CF.PCC.pval(cor.roots)
    print("P-values calculated.")
    if(log.file==T){
      write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": P-values calculated."), file=fileConn, append=T)

      }

    cor.signmat <- get.signMat(pcc)
    print("Sign matrix calculated.")
    if(log.file==T){
      write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": Sign matrix calculated."), file=fileConn, append=T)

      }

    equivalent.pccs <- get.inferred.PCCs(cor.p, cor.signmat, residuals$num.cells, num.genes)
    print("Equivalent PCCs calculated.")
    if(log.file==T){
      write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": Equivalent PCCs calculated."), file=fileConn, append=T)

      }
    sig.equivalent.pccs <- get.significant.inferred.PCCs(cor.p, equivalent.pccs, num.genes, cor.alpha, return.ps)
    if(is.list(sig.equivalent.pccs)==T){print(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": ", paste0("Number of remaining correlations:", Matrix::nnzero(sig.equivalent.pccs[[1]]))), fileConn)}
    else{
      print(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": ", paste0("Number of remaining correlations:", Matrix::nnzero(sig.equivalent.pccs))),fileConn)
    }
    if(log.file==T){
      write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": Equivalent PCCs filtered for significance."), file=fileConn, append=T)

      if(is.list(sig.equivalent.pccs)==T){writeLines(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": ", paste0("Number of remaining correlations:", Matrix::nnzero(sig.equivalent.pccs[[1]]))), fileConn)}
      else{
        write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": ", paste0("Number of remaining correlations:", Matrix::nnzero(sig.equivalent.pccs))), file=fileConn, append=T)

        }
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
  print("Pipeline complete.")
  if(log.file==T){
    write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": Pipeline complete."), file=fileConn, append=T)

    close(fileConn)
  }

}

