get.significant.inferred.PCCs <- function(p.matrix, inferred.PCCs, num.genes, alpha){

  print("Calculating significance for inferred PCCs.")
  sorted.pvals <- p.matrix[order(p.matrix[,"p.mpfr"], decreasing = T),]

  BH.corrections <- (10^-sorted.pvals[,"p.mpfr"])*((num.genes*(num.genes-1))/2)/c(1:length(sorted.pvals[,"p.mpfr"]))
  BH.corrections.mat <-  Matrix::sparseMatrix(i = sorted.pvals[,"row"], j = sorted.pvals[,"col"], x = BH.corrections, dims= dim(inferred.PCCs))

  BH.bool.mat <- apply(BH.corrections.mat, 2, function(x) ifelse(x < alpha, T, F))

  sig.PCCs <- inferred.PCCs * BH.bool.mat

  print("Done.")

  return(sig.PCCs)
}

get.inferred.PCCs.2 <- function(p.matrix, signMat,n, num.genes){

  print("Calculating inferred Pearson Correlation Coefficients.")

  p.sparse.mat <- Matrix::sparseMatrix(i = p.matrix[,"row"],
                                       j = p.matrix[,"col"],
                                       x = p.matrix[,"p.mpfr"],
                                       dims = c(num.genes, num.genes))

  erfcinv.ps <- erfcinv(10^-p.matrix[,"p.mpfr"])

  erfcinv.matrix <- Matrix::sparseMatrix(i = p.matrix[,"row"],
                                         j = p.matrix[,"col"],
                                         x = erfcinv.ps,
                                         dims = c(num.genes, num.genes))

  erfcinv.PCCs <- signMat*tanh((sqrt(2)*erfcinv.matrix)/sqrt(n-3))

  boolmat <- apply(erfcinv.PCCs, 2, function(x) ifelse(abs(x)== 1, T, F))

  estimate.mat <- p.sparse.mat*boolmat

  estimate.mat <- signMat * tanh(sqrt(2*log(10))/sqrt(n-3)*sqrt(estimate.mat))

  return.mat <- erfcinv.PCCs + estimate.mat - erfcinv.PCCs*boolmat


  print("Done.")
  return(return.mat)


}

get.signMat <- function(mcPCCs){
  signMat <- mcPCCs
  signMat[signMat < 0] <- -1
  signMat[signMat > 0] <- 1

  signMat[upper.tri(signMat)] <- 0

  diag(signMat) <- 0

  signMat <- as(signMat, "sparseMatrix")

  return(signMat)
}
