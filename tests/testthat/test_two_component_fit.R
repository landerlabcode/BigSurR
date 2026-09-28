#Testing various parts of the two-component fit using a subset of data (x) from

#Test that residuals are correct dimensions
testData <- readRDS(test_path("fixtures","test_counts.rds"))
testData <- as(testData, "CsparseMatrix")
testSObj <- CreateSeuratObject(counts=testData)

residuals <- get.residuals(testSObj, "RNA", "counts", "MeanSpecific")

test_that("Residuals dimensions match expectations", {
  expect_equal(dim(testData), dim(residuals$residuals))
})

#Test that number of mcFanos is correct
c <- residuals$c
num.genes <- residuals$num.genes

test_that("Number of c's matches number of genes", {
  expect_equal(num.genes, length(c))
})
mcfanos <- get.mcFanos(residuals)

test_that("Number of mcFanos matches number of genes", {
  expect_equal(num.genes, length(mcfanos))
})

#Test that the dimensions propagate correctly to the p-values
fano.cumulants <- Cumulants.Fano(residuals, c)
fanocoeffs <- CF.Coefficients.Fano(fano.cumulants[,2], fano.cumulants[,3], fano.cumulants[,4], fano.cumulants[,5], mcfanos, rownames(residuals$ematrix))
fanoroots <- CF.AllRoots(fanocoeffs)

test_that("Number of p-values matches number of genes", {
  expect_equal(num.genes, length(fanoroots))
})

#Check the inverse fano correction applies correctly
inv.correction <- inv.sqrt.correction2(residuals, residuals$eta, residuals$theta)

test_that("inv corrections produce 4 moments per grid point", {
  np <- length(inv.correction$points)          # not the literal 8
  expect_length(inv.correction$moments, np)
  expect_true(all(vapply(inv.correction$moments, length, integer(1)) == 4L))
  expect_true(all(is.finite(unlist(inv.correction$moments))))
  expect_false(is.unsorted(inv.correction$points))
  expect_equal(anyDuplicated(inv.correction$points), 0)
})

#Test that cumulants come out with correct shape and values
points <- inv.sqrt.points(residuals$gene.totals, residuals$num.cells)
moment.interp <- inv.sqrt.moment.interpolation2(inv.correction, residuals$gene.totals)

cor.cumulants <- Cumulants.PCC(residuals, moment.interp)
test_that("Corr. cumulants have appropriate shape and values",
          {
            m <- residuals$num.genes
            dims <- vapply(cor.cumulants, dim, integer(2))   # 2 x 4 matrix
            expect_true(all(dims == m))
            expect_equal(dim(dims), c(2L, 4L))
            for (k in seq_along(cor.cumulants))
              expect_true(isSymmetric(cor.cumulants[[k]]))
            expect_true(all(cor.cumulants[[1]] > 0))
          })


#Test that coefficients have correct shape and values
pcc <- get.mcPCC2(residuals, mcfanos)
cor.coefficients <- CF.Coefficients.PCC(cor.cumulants, pcc)
test_that("Cornish Fisher coefficients have appropriate shape and values", {
  m <- residuals$num.genes
  n.pairs <- m * (m - 1) / 2

  expect_true(is.matrix(cor.coefficients))
  expect_equal(dim(cor.coefficients), c(n.pairs, 7L))
  expect_equal(colnames(cor.coefficients),
               c("row", "col", "c1", "c2", "c3", "c4", "c5"))
  expect_true(all(is.finite(cor.coefficients)))

  # the row/col construction is the fiddly part
  expect_true(all(cor.coefficients[, "row"] > cor.coefficients[, "col"]))
  expect_true(all(cor.coefficients[, "row"] <= m))
  expect_true(all(cor.coefficients[, "col"] >= 1))
  expect_equal(anyDuplicated(paste(cor.coefficients[, "row"],
                                   cor.coefficients[, "col"])), 0)

})
test_that("flattening keeps coefficients aligned with their gene pair", {
  k2 <- cor.cumulants[[1]]; k3 <- cor.cumulants[[2]]
  k4 <- cor.cumulants[[3]]; k5 <- cor.cumulants[[4]]
  c1.full <- -pcc - k3/(6*k2) + 17*k3^3/(324*k2^4) -
    k3*k4/(12*k2^3) + k5/(40*k2^2)
  for (idx in c(1L, 7L, 30L, nrow(cor.coefficients))) {
    i <- cor.coefficients[idx, "row"]; j <- cor.coefficients[idx, "col"]
    expect_equal(unname(cor.coefficients[idx, "c1"]), c1.full[i, j],
                 label = sprintf("pair (%d,%d) at row %d", i, j, idx))
  }
})
