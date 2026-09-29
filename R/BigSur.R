#' BigSur (Basic Informatics and Gene Statistics from Unnormalized Reads)
#'
#' @param seurat.obj Seurat object containing the raw transcript counts filtered for zero count genes.
#' @param assay Assay slot containing raw transcript counts (default "RNA").
#' @param counts.slot Slot within assay containing raw counts matrix (default "counts").
#' @param cv.est.method String. Sets the method used for determining the coefficient of variation used to define null distributions. There are three options here: \n1) "Single": Determines a scalar value of c. \n2) "MeanSpecific": Estimates the relationship between mean expression and the expected coefficient of variation and predicts a null value for each gene. \n3)"TwoComponent": Estimates a the relationship between mean expression and the expected coefficient of variation using a two parameter fit.
#' @param variable.features Boolean. If true, BigSur will identify select variable features based on the modified corrected Fano factor.
#' @param correlations Boolean. If true, BigSur will identify statistically significant gene-gene correlations.
#' @param first.pass.cutoff Integer. Removes roots before p-value calculations if the root is below Abs[Sqrt(2)*InverseErfc(2*10^-first.pass.cutoff)]. The higher the number, the more correlations are removed in initial screening.
#' @param inverse.fano.moments Boolean. If true, BigSur will calculate the moments for the inverse Fano factor pairs before performing Cornish Fisher expansion.
#' @param fano.alpha Double. Desired false discovery cutoff for labeling of variable features. (Default 0.05).
#' @param min.fano Double. Minimum mcFano value considered for variable genes.
#' @param cor.alpha Double. Desired false discovery cutoff for labeling of statistically significant correlations.
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
                   cv.est.method="MeanSpecific",
                   variable.features=T,
                   correlations=F,
                   first.pass.cutoff=2,
                   inverse.fano.moments = T,
                   fano.alpha = 0.05,
                   min.fano = 1.5,
                   cor.alpha = 0.05,
                   return.ps = T,
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

  if (log.file) {
    fileConn <- file(log.file.dir, open = "at")
    on.exit(close(fileConn), add = TRUE)
    write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": Pipeline started execution."),
          file = fileConn)
  }

  residuals<-get.residuals(seurat.obj, assay, counts.slot, cv.est.method)
  new.seurat.obj <- seurat.obj
  new.seurat.obj[[assay]]$data <- residuals$residuals
  Misc(new.seurat.obj, slot="BigSur.Eta")<-residuals$eta
  Misc(new.seurat.obj, slot="BigSur.Theta")<-residuals$theta
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
    feat.metadata <- as.data.frame(fano.selected[,c(1,4,5)])
    new.seurat.obj[[assay]]@meta.data <- feat.metadata
    VariableFeatures(new.seurat.obj) <- top.features

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
      inv.correction <- inv.sqrt.correction2(residuals, residuals$eta, residuals$theta)
      moment.interp <- inv.sqrt.moment.interpolation2(inv.correction, residuals$gene.totals)
      print("Inverse sqrt moments calculated.")
      if(log.file==T){
        write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": Inverse sqrt moments calculated."), file=fileConn, append=T)
        }
    }

    else{
      onesmat <- matrix(1, nrow=num.genes, ncol=num.genes)
      moment.interp <- list(onesmat, onesmat, onesmat, onesmat)
    }

    #cor.cumulants <- Cumulants.PCC(residuals, moment.interp)
    #print("PCC cumulants calculated.")
    #if(log.file==T){
    # write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": PCC cumulants calculated."), file=fileConn, append=T)

    # }

    #cor.coefficients <- CF.Coefficients.PCC(cor.cumulants, pcc)
    #print("PCC Cornish Fisher coefficients calculated.")
   # if(log.file==T){
    #  write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": PCC Cornish Fisher coefficients calculated."), file=fileConn, append=T)
    #
     #  }

    cor.coefficients <-  CF.PCC.blocked(residuals, moment.interp, pcc, first.pass.cutoff)
    print("Correlation cumulants calculated.")
    if(log.file==T){
    write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": Correlation cumulants calculated."), file=fileConn, append=T)
    }

    cor.roots <- CF.PCC.Roots2(cor.coefficients, first.pass.cutoff)

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

    sig.pccs <- get.significant.PCCs(pcc, cor.p, num.genes, cor.alpha)
    print(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": ", paste0("Number of remaining correlations:", Matrix::nnzero(sig.pccs[[1]])/2)))
    if(log.file==T){
      write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": modified-corrected PCCs filtered for significance."), file=fileConn, append=T)
      writeLines(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": ", paste0("Number of remaining correlations:", Matrix::nnzero(sig.pccs[[1]])/2)), fileConn)

    }
    Misc(new.seurat.obj, slot="BigSur.Correlations") <- sig.pccs$pccs
    Misc(new.seurat.obj, slot="BigSur.log.adj.pvalues") <- sig.pccs$logp
    Misc(new.seurat.obj, slot="BigSur.orig.alpha") <- sig.pccs$alpha
}

  print("Pipeline complete.")
  if(log.file==T){
    write(paste0(format(Sys.time(), "%a %b %d %X %Y"), ": Pipeline complete."), file=fileConn, append=T)
  }
  return(new.seurat.obj)
}

