#testData <- readRDS(test_path("fixtures","test_counts.rds"))
testMat   <- readMM(test_path("fixtures","rajSubs.mtx"))
testGenes <- read.csv(test_path("fixtures","rajSubGenes.csv"), header = FALSE)[[1]]
testCells <- read.csv(test_path("fixtures","rajSubCells.csv"), header = FALSE)[[1]]

dimnames(testMat) <- list(testGenes, testCells)

testData <- as(testMat, "CsparseMatrix")
testSObj <- CreateSeuratObject(counts=testData)

residuals <- get.residuals(testSObj, "RNA", "counts", "MeanSpecific")
residualsTwoC <-  get.residuals(testSObj, "RNA", "counts", "TwoComponent")
readref <- function(f) unname(as.matrix(read.csv(f, header = FALSE)))

cumFunc <- function(res){
  c <- res$c
  mcfanos <- get.mcFanos(res)
  fano.cumulants <- Cumulants.Fano(res, c)
  fanocoeffs <- CF.Coefficients.Fano(fano.cumulants[,2], fano.cumulants[,3], fano.cumulants[,4], fano.cumulants[,5], mcfanos, rownames(res$ematrix))

  pcc <- get.mcPCC2(res, mcfanos)

  inv.correction <- inv.sqrt.correction2(res, res$eta, res$theta)
  moment.interp <- inv.sqrt.moment.interpolation2(inv.correction, res$gene.totals)
  cor.cumulants <- Cumulants.PCC(res, moment.interp)
  cor.coefficients <- CF.Coefficients.PCC(cor.cumulants, pcc)

  list(mcFano=mcfanos, FanoCum=fano.cumulants, FanoCoef=fanocoeffs, pccs = pcc,
       InvCorrection=inv.correction, CorrCum = cor.cumulants)
}

fileHead <- "../mathematica/RajThetaOnly/"
#Residuals
test_that("Residuals match between methods.",
          {
            #Residuals
            r <- readref(paste0(fileHead, "residuals.csv"))
            expect_equal(unname(residuals$residuals), unname(r), tolerance=0.0001)
          })

rOut <- cumFunc(residuals)

mcpcc <- readref("/Users/bigcomputer/Desktop/UCI/LanderLab/BigSur/tests/mathematica/RajThetaOnly/pccmatrix.csv")

is_lower_triangular <- function(mat) {
  if (!is.matrix(mat) || nrow(mat) != ncol(mat)) return(FALSE)
  # Check if all elements strictly above the diagonal are 0
  all(mat[!lower.tri(mat, diag = TRUE)] == 0)
}

is_lower_triangular(mcpcc)
is_lower_triangular(rOut$pccs)


x <- unname(rOut$pccs); y <- mcpcc
keep <- residuals$gene.totals > 1
r <- (x/y)[keep, keep]
gt <- residuals$gene.totals[keep]

par(mfrow = c(1, 3))

# 1. distribution — where is the disagreement centred?
hist(log10(abs(as.vector(r))), breaks = 100,
     main = "log10 ratio", xlab = "")
abline(v = 0, col = "red")

# 2. does it track expression? (points at the interpolation grid edges show here)
plot(gt, r[, 1], log = "x", pch = ".", main = "ratio vs gene total")
abline(h = 1, col = "red")

# 3. is it structured in the matrix? (blocks = a subset of genes; stripes = one gene)
image(log10(abs(r)), main = "ratio, spatial")

d <- abs(x - y)
c(max = max(d), median = median(d), q99 = quantile(d, 0.99), n_over_1e6 = sum(d > 1e-6))

ck1 <- readref("/Users/bigcomputer/Desktop/UCI/LanderLab/BigSur/tests/mathematica/RajThetaOnly/ck1.csv")
ck2 <- readref("/Users/bigcomputer/Desktop/UCI/LanderLab/BigSur/tests/mathematica/RajThetaOnly/ck2.csv")
ck3 <- readref("/Users/bigcomputer/Desktop/UCI/LanderLab/BigSur/tests/mathematica/RajThetaOnly/ck3.csv")
ck4 <- readref("/Users/bigcomputer/Desktop/UCI/LanderLab/BigSur/tests/mathematica/RajThetaOnly/ck4.csv")
keep <- residuals$gene.totals > 1

for (k in 1:4) {
  x <- as.matrix(rOut$CorrCum[[k]])[keep, keep]
  y <- as.matrix(get(paste0("ck", k)))[keep, keep]
  s <- max(abs(y))
  d <- abs(x - y)
  cat(sprintf("CK%d  scale=%.3e  max/scale=%.3e  median/scale=%.3e  rmse/scale=%.3e\n",
              k, s, max(d)/s, median(d)/s, sqrt(mean(d^2))/s))
}

k <- 1
x <- as.matrix(rOut$CorrCum[[k]])[keep, keep]
y <- as.matrix(ck1)[keep, keep]
i <- which(abs(x - y) == max(abs(x - y)), arr.ind = TRUE)
i
c(R = x[i[1,1], i[1,2]], M = y[i[1,1], i[1,2]])
sum(abs(x - y) > 0.01 * max(abs(y)))


g   <- rownames(residuals$ematrix)[keep]
gtk <- residuals$gene.totals[keep]

fR <- sqrt(abs(diag(x))); fM <- sqrt(abs(diag(y)))
r  <- fM / fR
bad <- which(abs(r - 1) > 0.05)

data.frame(gene = g[bad], total = gtk[bad], R = fR[bad], M = fM[bad], ratio = r[bad])

c(range_bad = range(gtk[bad]), grid = range(rOut$InvCorrection$points),
  above = sum(gtk[bad] > max(rOut$InvCorrection$points)))

c(n_bad = length(bad), n_total = length(gtk))
quantile(gtk, c(0, .05, .1, .25, .5))        # where does 295 sit in the distribution?
sum(gtk <= 295)
o <- order(gtk)
plot(gtk[o], (fM/fR)[o], log = "x", pch = ".", ylim = c(0, 5))
abline(h = 1, col = "red")
abline(v = rOut$InvCorrection$points, col = "blue", lty = 3)

rOut$InvCorrection$trials
round(4e7 / (residuals$num.cells * (log10(rOut$InvCorrection$points)^(1/5) +
                                      0.5*log10(rOut$InvCorrection$points)^3)))

ic_mma_c <- inv.sqrt.correction2(residuals, eta = 0, theta = 0.261967)
mi <- inv.sqrt.moment.interpolation2(ic_mma_c, residuals$gene.totals)

cc2  <- Cumulants.PCC(residuals, mi)

x2 <- as.matrix(cc2[[1]])[keep, keep]
y  <- as.matrix(ck1)[keep, keep]

fR2 <- sqrt(abs(diag(x2))); fM <- sqrt(abs(diag(y)))
o <- order(gtk)
plot(gtk[o], (fM/fR2)[o], log = "x", pch = ".", ylim = c(0, 5))
abline(h = 1, col = "red")

test_that("mcFanos are consistent between versions",
          {
            #mcFanos
            mf <-readref(paste0(fileHead, "mcfanos.csv"))
            expect_equal(unname(rOut$mcFano), as.vector(unname(mf)), tolerance=0.0001)
})
test_that("Fano cumulants are consistent between versions",{
            #Fano cumulants
            fk1<-readref(paste0(fileHead, "FK1.csv"))
            fk2<-readref(paste0(fileHead, "FK2.csv"))
            fk3<-readref(paste0(fileHead, "FK3.csv"))
            fk4<-readref(paste0(fileHead, "FK4.csv"))

            expect_equal(unname(rOut$FanoCum[,2]), as.vector(unname(fk1)), tolerance=0.0001)
            expect_equal(unname(rOut$FanoCum[,3]), as.vector(unname(fk2)), tolerance=0.0001)
            expect_equal(unname(rOut$FanoCum[,4]), as.vector(unname(fk3)), tolerance=0.0001)
            expect_equal(unname(rOut$FanoCum[,5]), as.vector(unname(fk4)), tolerance=0.0001)
})

test_that("Fano CF coefficients are consistent between versions",{
            #Fano Cornish fisher coefficients
           fcf1 <-readref(paste0(fileHead, "FCF1.csv"))
           fcf2 <-readref(paste0(fileHead, "FCF2.csv"))
           fcf3 <-readref(paste0(fileHead, "FCF3.csv"))
           fcf4 <-readref(paste0(fileHead, "FCF4.csv"))
           fcf5 <-readref(paste0(fileHead, "FCF5.csv"))
           expect_equal(unname(rOut$FanoCoef[,1]), as.vector(unname(fcf1)), tolerance = 0.0001)
           expect_equal(unname(rOut$FanoCoef[,2]), as.vector(unname(fcf2)), tolerance = 0.0001)
           expect_equal(unname(rOut$FanoCoef[,3]), as.vector(unname(fcf3)), tolerance = 0.0001)
           expect_equal(unname(rOut$FanoCoef[,4]), as.vector(unname(fcf4)), tolerance = 0.0001)
           expect_equal(unname(rOut$FanoCoef[,5]), as.vector(unname(fcf5)), tolerance = 0.0001)
})

test_that("mcPCCs are consistent",{
            #mcPCCs
           mcpcc <- readref(paste0(fileHead, "pccmatrix.csv"))
           expect_equal(unname(rOut$pccs), as.matrix(unname(mcpcc)), tolerance=0.0001)
})

test_that("mcPCC cumulants are consistent",{
            #mcPCC cumulants
           ck1 <- readref(paste0(fileHead, "ck1.csv"))
           ck2 <- readref(paste0(fileHead, "ck2.csv"))
           ck3 <- readref(paste0(fileHead, "ck3.csv"))
           ck4 <- readref(paste0(fileHead, "ck4.csv"))

           expect_equal(unname(as.matrix(rOut$CorrCum[[1]])),
                        unname(as.matrix(ck1)), tolerance = 0.05)
           expect_equal(unname(as.matrix(rOut$CorrCum[[2]])),
                        unname(as.matrix(ck2)), tolerance = 0.05)
           expect_equal(unname(as.matrix(rOut$CorrCum[[3]])),
                        unname(as.matrix(ck3)), tolerance = 0.05)
           expect_equal(unname(as.matrix(rOut$CorrCum[[4]])),
                        unname(as.matrix(ck4)), tolerance = 0.05)

          })
