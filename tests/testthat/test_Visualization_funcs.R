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
