library(Seurat)
library(BigSur)

data.dir <- file.path(tempdir(), "pbmc3k")
if (!dir.exists(data.dir)) {
  url <- "https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"
  tf <- tempfile(fileext = ".tar.gz")
  download.file(url, tf)
  untar(tf, exdir = data.dir)
}
pbmc_data <- Read10X(file.path(data.dir, "filtered_gene_bc_matrices/hg19"))
pbmc <- CreateSeuratObject(pbmc_data, features=names(keep_genes[keep_genes]))
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 1000 & percent.mt < 5)

#We need to also remove any genes with 0 UMI counts.
counts <- Seurat::GetAssayData(pbmc, assay = "RNA", layer = "counts")
keep_genes <- Matrix::rowSums(counts) > 10
pbmc <- subset(pbmc, features = names(keep_genes[keep_genes]))

set.seed(117)
pbmc_BS <- BigSur(pbmc, correlations=T, null.distribution = "NB")

corr.matrix <- pbmc_BS@misc$BigSur.Correlations
pvals <- pbmc_BS@misc$BigSur.log.adj.pvalues

s <- Matrix::summary(corr.matrix)
pairs <- data.frame(
  row   = rownames(corr.matrix)[s$i],
  col   = colnames(corr.matrix)[s$j],
  value = s$x
)

genes <- c("CD4", "IL2", "RORC", "FOXP3","CCR7", "CD8A", "CD8B", "LEF1", "TCF7", "IL7R", "GZMB", "PRF1", "NKG7", "GNLY")
pairs[pairs$row %in% genes | pairs$col %in% genes, ]
check <- c("GNLY", "GZMB", "PRF1", "CD8B", "LEF1", "CCR7", "SELL")
pairs[(pairs$row %in% c("NKG7", "CD8A", "TCF7") & pairs$col %in% check) |
        (pairs$col %in% c("NKG7", "CD8A", "TCF7") & pairs$row %in% check), ]
head(pairs[order(-abs(pairs$value)), ], 40)

pbmc_Corrs <- pbmc_BS@misc$BigSur.Correlations
set.seed(117)
pbmc_Modules <- FindCorrelationModules(pbmc_Corrs)

pbmc_Modules_Merged <- MergeSmallModules(pbmc_Corrs, pbmc_Modules, min.size=15)


sizes(pbmc_Modules_Merged)

mods <- igraph::communities(pbmc_Modules)
mods <- mods[order(-sizes(pbmc_Modules))]
sort(mods[[5]])

mod1_Genes <- sort(mods[[5]])
mod1_Genes <- mod1_Genes[!startsWith(mod1_Genes, "RP")]
mod1_Matrix <- SubsetCorrelationMatrix(pbmc_BS, pCutoff=0.001, minPosCorr = 0.2, minNegCorr = -0.1, genes = mod1_Genes)
mod1_Plot <- StaticCorrelationPlot(mod1_Matrix, highlight=c("CCR7", "SELL", "LEF1", "TCF7", "IL7R", "KLF2", "FOXP1"))
mod1_Plot
ggsave("~/Desktop/mod1_Plot.png", mod1_Plot, width = 24, height = 24, limitsize = FALSE)
InterCommunityCorrelations(pbmc_Corrs, mods, span=1:10, max.pos.fraction = 0.02, max.neg.fraction = 0.01)

sort(mods[[1]])

mod3_Genes <- sort(mods[[1]])
mod3_Genes <- mod3_Genes[!startsWith(mod3_Genes, "RP")]

mod3_mod5_Genes <- union(mod3_Genes, mod1_Genes)
mod5_mod3_Matrix <- SubsetCorrelationMatrix(pbmc_BS, pCutoff=0.001, minPosCorr = 0.3, minNegCorr = 0, genes = mod3_mod5_Genes)
mod35_Plot <- StaticCorrelationPlot(mod5_mod3_Matrix, highlight=c("CCR7", "SELL", "LEF1", "TCF7", "IL7R", "KLF2", "FOXP1", "NKG7","GNLY", "GZMA", "KLRF1", "NCAM1"))
ggsave("~/Desktop/mod35_Plot.png", mod35_Plot, width = 24, height = 24, limitsize = FALSE)
