#Test that pipeline function actually runs on 2000 genes.

testMat   <- readMM(test_path("fixtures","rajSubs.mtx"))[1:2000,]
testGenes <- read.csv(test_path("fixtures","rajSubGenes.csv"), header = FALSE)[[1]][1:2000]
testCells <- read.csv(test_path("fixtures","rajSubCells.csv"), header = FALSE)[[1]]

dimnames(testMat) <- list(testGenes, testCells)

testData <- as(testMat, "CsparseMatrix")
testSObj <- CreateSeuratObject(counts=testData)

outObj <- BigSur(testSObj, correlations=T, cor.alpha=0.02)

dim(outObj@misc$BigSur.Correlations)
dim(outObj@misc$BigSur.log.adj.pvalues)

a <- BigSur(testSObj, correlations = TRUE, cor.alpha = 0.02)
b <- BigSur(testSObj, correlations = TRUE, cor.alpha = 0.02)
c(length(a@misc$BigSur.Correlations@x), length(b@misc$BigSur.Correlations@x))

outObj@misc$BigSur.Theta


#Test that output is consistent with the Poisson log-normal version of the code.
x <- outObj@misc$BigSur.Correlations@x
c(total = length(x), negative = sum(x < 0), positive = sum(x > 0))
writeMM(outObj@misc$BigSur.Correlations, test_path("..", "mathematica", "rCorr2.mtx"))
writeMM(outObj@misc$BigSur.log.adj.pvalues, test_path("..", "mathematica", "rPs2.mtx"))
?writeMM
write.csv(outObj@assays$RNA@layers$data, "rRes.csv")

#The mcPCCs are not matching between the two. Let's try to isolate where it goes wrong.
residuals <- get.residuals(testSObj, "RNA", "counts", "MeanSpecific")
mcfanos <- get.mcFanos(residuals)
pcc <- get.mcPCC2(residuals, mcfanos)

write.csv(pcc, test_path("..", "mathematica", "rAllPCC.csv"))
write.csv(mcfanos, test_path("..", "mathematica", "rAllmcFanos.csv"))

identical(residuals$residuals, as.matrix(outObj[["RNA"]]$data))
identical(names(mcfanos), rownames(residuals$residuals))

P <- residuals$residuals / sqrt((residuals$num.cells - 1) * mcfanos)
c(hand = sum(P[18, ] * P[11, ]), func = pcc[18, 11])

inv.correction <- inv.sqrt.correction2(residuals, residuals$eta, residuals$theta)
moment.interp <- inv.sqrt.moment.interpolation2(inv.correction, residuals$gene.totals)
cor.coefficients <-  CF.PCC.blocked(residuals, moment.interp, pcc, 2)
cor.roots <- CF.PCC.Roots2(cor.coefficients, 2)
cor.p <- CF.PCC.pval(cor.roots)

sapply(c(0.05, 0.02, 0.01), function(al)
  sum(Matrix::summary(get.significant.PCCs(pcc, cor.p, residuals$num.genes, al)$pccs)$x != 0))

cmatrix<-cor.coefficients
np <- vapply(seq_len(nrow(cmatrix)), function(i) {
  r <- polyroot(cmatrix[i, 3:7]); sum(Re(r[abs(Im(r)) < 1e-5]) > 0)
}, numeric(1))
table(np)

roots.matrix <- cor.roots
two <- which(np == 2)
mma <- vapply(two, function(i) {
  r <- polyroot(cmatrix[i, 3:7]); r <- Re(r[abs(Im(r)) < 1e-5]); min(r[r > 0])
}, numeric(1))
mine <- roots.matrix[two, 3]
c(agree = sum(abs(mma - mine) < 1e-8), differ = sum(abs(mma - mine) >= 1e-8))
summary(abs(mine) - abs(mma))

c(n_tested = length(cor.p[, "logp"]), all_pairs = choose(residuals$num.genes, 2))


#P-values are different between the two
lp <- -cor.p[, "logp"]/log(10)        # convert to -log10(p), same units
c(n = length(lp))
quantile(lp, c(0, .25, .5, .75, 1))

logp <- cor.p[, "logp"]
n <- length(logp)
o <- order(logp)
adj <- logp[o] + log(choose(residuals$num.genes, 2)) - log(seq_len(n))
adj <- rev(cummin(rev(adj)))
adj.orig <- numeric(n); adj.orig[o] <- adj
keep <- adj.orig < log(0.02)

c(n = length(cor.p[,"logp"]),
  cutoff_log = max(adj.orig[keep]),
  threshold = log(0.02))

n <- length(cor.p[, "logp"]); o <- order(cor.p[, "logp"])
sorted <- cor.p[o, "logp"]
cond <- sorted <= log(0.02) + log(seq_len(n)) - log(choose(residuals$num.genes, 2))
c(raw_cutoff_R = if (any(cond)) sorted[max(which(cond))] else -Inf,
  raw_cutoff_M = -10.3757,
  n_passing_R  = if (any(cond)) max(which(cond)) else 0)
