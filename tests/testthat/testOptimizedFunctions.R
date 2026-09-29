testMat   <- readMM(test_path("fixtures","rajSubs.mtx"))[1:2000,]
testGenes <- read.csv(test_path("fixtures","rajSubGenes.csv"), header = FALSE)[[1]][1:2000]
testCells <- read.csv(test_path("fixtures","rajSubCells.csv"), header = FALSE)[[1]]

dimnames(testMat) <- list(testGenes, testCells)

testData <- as(testMat, "CsparseMatrix")
testSObj <- CreateSeuratObject(counts=testData)

residuals <- get.residuals(testSObj, "RNA", "counts", "MeanSpecific")

c <- residuals$c
mcfanos <- get.mcFanos(residuals)
pcc <- get.mcPCC2(residuals, mcfanos)

inv.correction <- inv.sqrt.correction2(residuals, residuals$eta, residuals$theta)
moment.interp <- inv.sqrt.moment.interpolation2(inv.correction, residuals$gene.totals)

first.pass.cutoff <- 2

dim(cor.cumulants[[1]])
dim(pcc)
length(residuals$c)

length(cor.cumulants)
sapply(cor.cumulants, function(k) if (is.null(k)) "NULL" else paste(dim(k), collapse = "x"))

cor.coefficients <- CF.Coefficients.PCC(cor.cumulants, pcc)

sapply(cor.cumulants, function(k) sum(!is.finite(k)))

#Old version
system.time({cor.cumulants <- Cumulants.PCC(residuals, moment.interp)
cor.coefficients <- CF.Coefficients.PCC(cor.cumulants, pcc)

cmatrix.pruned.1 <- QuickTest6CF(cor.coefficients, first.pass.cutoff)

to.test.bool <- apply(cmatrix.pruned.1, 1,
                      function(x) ifelse((residuals$gene.totals[x[1]]<=84)|(residuals$gene.totals[x[2]]<=84),
                                         T, F))

cmatrix.pruned.1 <- cbind(cmatrix.pruned.1, to.test.bool)

cmatrix.pruned.2 <- cmatrix.pruned.1[cmatrix.pruned.1[,8]==F, ]

cmatrix.more.testing <- cmatrix.pruned.1[cmatrix.pruned.1[,8]==T, ]

cmatrix.passed <- SecondTestCF(cmatrix.more.testing, first.pass.cutoff)

cmatrix.pruned.2 <- rbind(cmatrix.pruned.2, cmatrix.passed)})

gc()
#New version
system.time({cmatrix.new <- CF.PCC.blocked(residuals, moment.interp, pcc, first.pass.cutoff)})

#Compare values
key <- function(m) paste(m[,1], m[,2], sep = ":")
kOld <- key(cmatrix.pruned.2); kNew <- key(cmatrix.new)

old<-cmatrix.pruned.2
new<-cmatrix.new

c(old = length(kOld), new = length(kNew),
  old_not_in_new = length(setdiff(kOld, kNew)),   # must be 0
  new_only       = length(setdiff(kNew, kOld)))

head(old[, 1:2]); head(new[, 1:2])
c(old_row_range = range(old[,1]), old_col_range = range(old[,2]),
  new_row_range = range(new[,1]), new_col_range = range(new[,2]))
head(kOld); head(kNew)


#Test roots
system.time(old.roots <- CF.PCC.Roots(cmatrix.pruned.2, first.pass.cutoff, residuals$gene.totals))
system.time(new.roots <- CF.PCC.Roots2(cmatrix.pruned.2, first.pass.cutoff))
table(sign(new.roots[,3]))
sum(abs(new.roots[,3]) != abs(old.roots[,3]))   # magnitude changes, not just sign

hist(new.roots[,3])

hist(old.roots[,3])

d <- abs(abs(new.roots[,3]) - abs(old.roots[,3]))
quantile(d, c(.5, .9, .99, 1))
sum(d > 1e-6)

i <- order(d, decreasing = TRUE)[1:5]
for (k in i) {
  cf <- cmatrix.pruned.2[k, 3:7]
  r  <- polyroot(cf); r <- sort(Re(r[abs(Im(r)) < 1e-5]))
  cat(sprintf("c1=%9.4f  real roots: %-45s old=%8.4f new=%8.4f\n",
              cf[1], paste(round(r,3), collapse=" "), old.roots[k,3], new.roots[k,3]))
}

c(exactly_zero = sum(old.roots[,3] == 0),
  near_zero    = sum(abs(old.roots[,3]) < 1e-8),
  total        = nrow(old.roots))
quantile(abs(old.roots[,3]), c(0, .25, .5, .75, 1))

k <- rbind(c(0.295, 1.2, -0.05, 0.001, 1e-5),
           c(-0.340, 0.9,  0.02, 0.002, 2e-5))

s1 <- apply(k, 1, polyroot);  cat("stage1 dim:", dim(s1), "\n")
s2 <- apply(s1, 1, function(x) ifelse(abs(Im(x)) < 1e-5, Re(x), NA))
cat("stage2 dim:", dim(s2), "\n")
s3 <- apply(s2, 1, function(x) ifelse(!all(is.na(x)), min(abs(x), na.rm=TRUE), NA))
print(s3)
print(lapply(1:2, function(i) polyroot(k[i, ])))
