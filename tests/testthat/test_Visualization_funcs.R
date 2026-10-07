#Test that subsetting works
testMat   <- readMM(test_path("fixtures","rajSubs.mtx"))[1:2000,]
testGenes <- read.csv(test_path("fixtures","rajSubGenes.csv"), header = FALSE)[[1]][1:2000]
testCells <- read.csv(test_path("fixtures","rajSubCells.csv"), header = FALSE)[[1]]
dimnames(testMat) <- list(testGenes, testCells)
testData <- as(testMat, "CsparseMatrix")
testSObj <- CreateSeuratObject(counts=testData)

outObj <- BigSur(testSObj, correlations=T, cor.alpha=0.001)
sub <- SubsetCorrelationMatrix(outObj, pCutoff=0.001, minPosCorr=0.2, minNegCorr=-0.2, genes = F)

x <- outObj@misc$BigSur.Correlations@x
y <- sub@x
c(total = length(x), negative = sum(x < 0), positive = sum(x > 0))
c(total = length(y), negative = sum(y < 0), positive = sum(y > 0))

y

cs <- outObj@misc$BigSur.Correlations
sum(cs@x > 0)
cs[cs > 0 & cs < 0.00002] <- 0
cs <- Matrix::drop0(cs)
sum(cs@x > 0)

ps <- outObj@misc$BigSur.log.adj.pvalues
c(stored = length(ps@x), under_0.001 = sum(ps@x < log(0.001)))


#Test module finding
corrMat <- outObj@misc$BigSur.Correlations
testModules <- FindCorrelationModules(corrMat)


#Test igraph functions separately
corr.matrix <- corrMat

lt.corr <- Matrix::tril(corr.matrix, k = -1)
lt.corr@x[lt.corr@x<0] <- 0
lt.corr <- Matrix::drop0(lt.corr)

adj.graph <- graph_from_adjacency_matrix(
  lt.corr,
  mode="lower",
  diag = FALSE,
  weighted = TRUE)

E(adj.graph)$color <- ifelse(E(adj.graph)$weight >0, "#005AB5", "#DC3220")

lay <- LayoutFromPos(adj.graph)

subModules <- FindCorrelationModules(sub)
mergedModules <- MergeSmallModules(sub, subModules)

modules <- FindCorrelationModules(corrMat)
mergedAll <- MergeSmallModules(corrMat, modules)

p<-StaticCorrelationPlot(corrMat, highlight=c("BRCA1","BRCA2"), modules = mergedAll)
ggsave("~/Desktop/network.pdf", p, width = 24, height = 24, limitsize = FALSE)
StaticCorrelationPlot(sub, highlight=c("BRCA1","BRCA2"), modules = mergedModules)

ggraph(adj.graph, layout = "manual", x = lay[, 1], y = lay[, 2]) +
  geom_edge_link(aes(color = weight > 0), width = 0.6, alpha = 0.7) +
  scale_edge_color_manual(values = c(`TRUE` = "#005AB5", `FALSE` = "#DC3220"),
                          guide = "none") +
  geom_node_point(aes(fill = fill, shape = shape), size = 3,
                  color = "grey20", stroke = 0.3) +
  scale_fill_identity() +
  scale_shape_identity() +
  geom_node_text(aes(label = name), repel = TRUE, size = 2.5,
                 max.overlaps = Inf) +
  theme_void()

testModules[[3]]
