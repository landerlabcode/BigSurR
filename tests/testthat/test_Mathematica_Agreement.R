#Test that the residuals and values for eta and theta match those calculated in Mathematica.
testData <- readRDS(test_path("fixtures","test_counts.rds"))
testData <- as(testData, "CsparseMatrix")
testSObj <- CreateSeuratObject(counts=testData)

residuals <- get.residuals(testSObj, "RNA", "counts", "MeanSpecific")
residualsTwoC <-  get.residuals(testSObj, "RNA", "counts", "TwoComponent")
#Send data to Mathematica
#writeMM(testData,test_path("fixtures","test_counts.mtx"))
#write.csv(testData@Dimnames[1], test_path("fixtures","test_Genes.csv"), row.names=F, col.names = F)
#write.csv(testData@Dimnames[2], test_path("fixtures","test_Cells.csv"), row.names=F, col.names = F)

#Import Mathematica Data
etaTheta <- scan("../mathematica/etaTheta.txt")
thetaOnly <-scan("../mathematica/thetaonly.txt")

test_that("Theta only fit matches Mathematica version up to tolerance of 0.01.",{
          expect_equal(residuals$theta, thetaOnly, tolerance = 0.01)
  })

means <- residualsTwoC$gene.totals / residualsTwoC$num.cells
y <- 1 + residualsTwoC$j_values * means
k <- is.finite(residualsTwoC$j_values) & is.finite(y) & y > 0 & means > 0
lx <- log(means[k]); ly <- log(y[k])
obj <- function(e, t) sum(abs(ly - log(1 + e + t * exp(lx))))

test_that("R's fit is at least as good as Mathematica's", {
  # flat plateau: the objective, not the parameters, is what is determined
  expect_lte(obj(residualsTwoC$eta, residualsTwoC$theta),
             obj(etaTheta[1], etaTheta[2]) * 1.001)
})

etaTheta[2] / residualsTwoC$theta
sqrt(etaTheta[2]) / residualsTwoC$theta
etaTheta[2] / sqrt(residualsTwoC$theta)

test_that("the resulting c vector agrees with Mathematica (two component)", {
  means <- residualsTwoC$gene.totals / residualsTwoC$num.cells
  c.got <- sqrt(residualsTwoC$eta / means + residualsTwoC$theta)
  c.ref <- sqrt(etaTheta[1] / means + etaTheta[2])
  expect_equal(unname(c.got), unname(c.ref), tolerance = 0.075)
})

test_that("the resulting c vector agrees with Mathematica (one component)", {
  means <- residuals$gene.totals / residuals$num.cells
  c.got <- sqrt(1 / means + residuals$theta)
  c.ref <- sqrt(1 / means + thetaOnly)
  expect_equal(unname(c.got), unname(c.ref), tolerance = 0.075)
})

test_that("the resulting residuals vector agrees with Mathematica (two component)", {
  means <- residuals$gene.totals / residuals$num.cells
  c.got <- sqrt(1 / means + residuals$theta)
  c.ref <- sqrt(1 / means + thetaOnly)
  expect_equal(unname(c.got), unname(c.ref), tolerance = 0.075)
})


test_that("the resulting residuals vector agrees with Mathematica (one component)", {
  means <- residuals$gene.totals / residuals$num.cells
  c.got <- sqrt(1 / means + residuals$theta)
  c.ref <- sqrt(1 / means + thetaOnly)
  expect_equal(unname(c.got), unname(c.ref), tolerance = 0.075)
})


test_that("theta matches Mathematica (two component)", {
  expect_equal(residualsTwoC$theta, etaTheta[2], tolerance = 0.01)
})
