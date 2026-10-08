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
counts <- Seurat::GetAssayData(pbmc, assay = "RNA", layer = "counts")
keep_genes <- Matrix::rowSums(counts) > 0
pbmc <- CreateSeuratObject(pbmc_data,features = names(keep_genes[keep_genes]))
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)

#We need to also remove any genes with 0 UMI counts.
counts <- Seurat::GetAssayData(pbmc, assay = "RNA", layer = "counts")
keep_genes <- Matrix::rowSums(counts) > 0
pbmc <- subset(pbmc, features = names(keep_genes[keep_genes]))

pbmc_BS <- BigSur(pbmc, min.fano=1, fano.alpha = 0.05)

mcfanos <- pbmc_BS@assays$RNA@meta.data$mcfanos
mcfanoPvals <- pbmc_BS@assays$RNA@meta.data$BH.Corrected.Pvalue
genes <- rownames(pbmc_BS)

df <- data.frame(
  Gene <- genes,
  mcFano <- mcfanos,
  Adj.Pval <- mcfanoPvals
)

df_sorted <- df[order(df$mcFano, decreasing=T), ]
head(df_sorted, n=20)

length(VariableFeatures(pbmc_BS))

pbmc_BS <- BigSur(pbmc, min.fano=1.5, fano.alpha=0.05)
length(VariableFeatures(pbmc_BS))

pbmc_BS <- ScaleData(pbmc_BS)
pbmc_BS <- RunPCA(pbmc_BS)
set.seed(117)
pbmc_BS <- FindNeighbors(object = pbmc_BS, dims = 1:20)
pbmc_BS <- FindClusters(object = pbmc_BS, resolution=0.9)
pbmc_BS <- RunUMAP(object = pbmc_BS, dims = 1:20)
DimPlot(object = pbmc_BS, reduction = "umap")

pbmc <- NormalizeData(object = pbmc)
pbmc <- FindVariableFeatures(object = pbmc)
pbmc <- ScaleData(object = pbmc)
pbmc <- RunPCA(object = pbmc)
set.seed(117)
pbmc <- FindNeighbors(object = pbmc, dims = 1:20)
pbmc <- FindClusters(object = pbmc, resolution=0.9)
pbmc <- RunUMAP(object = pbmc, dims = 1:20)
DimPlot(object = pbmc, reduction = "umap")

Matrix::writeMM(pbmc_BS@assays$RNA@layers$counts, "/Users/bigcomputer/Desktop/UCI/LanderLab/PBMCTest/counts.mtx")
write.csv(rownames(pbmc_BS), "/Users/bigcomputer/Desktop/UCI/LanderLab/PBMCTest/genes.csv", col.names=F)
write.csv(colnames(pbmc_BS), "/Users/bigcomputer/Desktop/UCI/LanderLab/PBMCTest/cells.csv", col.names=F)

dim(pbmc_BS)

counts <- Seurat::GetAssayData(pbmc, assay = "RNA", layer = "counts")
keep_genes <- Matrix::rowSums(counts) > 0
pbmc <- CreateSeuratObject(pbmc_data,features = names(keep_genes[keep_genes]))
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- subset(pbmc, subset = nFeature_RNA > 1000 & percent.mt < 5)

counts <- GetAssayData(pbmc, layer = "counts")
num_cells <- Matrix::rowSums(counts > 0)
n <- 10
keep_genes <- names(num_cells[num_cells >= n])

pbmc <- subset(pbmc, features = keep_genes)
dim(pbmc)
pbmc_BS <- BigSur(pbmc, correlations=T, null.distribution = "NB")

corr.matrix <- pbmc_BS@misc$BigSur.Correlations
pvals <- pbmc_BS@misc$BigSur.log.adj.pvalues

identical(rownames(corr.matrix), rownames(pbmc_BS))

s <- Matrix::summary(corr.matrix)

pairs <- data.frame(
  row   = rownames(corr.matrix)[s$i],
  col   = colnames(corr.matrix)[s$j],
  value = s$x
)

det <- Matrix::rowSums(counts > 0)
top <- head(pairs[order(-abs(pairs$value)), ], 40)
top$n.row  <- det[top$row]
top$n.col  <- det[top$col]
top$n.both <- mapply(function(a, b) sum(counts[a, ] > 0 & counts[b, ] > 0),
                     top$row, top$col)
top

counts <- pbmc_BS

identical(rownames(counts), rownames(r))
identical(colnames(counts), colnames(r))
mean((counts[a, ] > 0) == (r[a, ] > 0))   # should be 1: residual > 0 exactly where the gene is detected

top <- head(pairs[order(-abs(pairs$value)), ], 20)
top$p <- mapply(function(a, b) pvals[a, b], top$row, top$col)
top

genes <- c("CD4", "IL2", "RORC", "FOXP3","CCR7", "CD8A", "CD8B", "LEF1", "TCF7", "IL7R", "GZMB", "PRF1", "NKG7", "GNLY")
pairs[pairs$row %in% genes | pairs$col %in% genes, ]

check <- c("GNLY", "GZMB", "PRF1", "CD8B", "LEF1", "CCR7", "SELL")
pairs[(pairs$row %in% c("NKG7", "CD8A", "TCF7") & pairs$col %in% check) |
        (pairs$col %in% c("NKG7", "CD8A", "TCF7") & pairs$row %in% check), ]
head(pairs[order(-abs(pairs$value)), ], 40)
