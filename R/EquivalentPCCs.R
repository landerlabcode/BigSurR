get.significant.PCCs <- function(mcpccs, p.matrix, num.genes, alpha){

  print("Calculating significance for modified-corrected PCCs.")
  logp <- p.matrix[, "logp"]
  n <- length(logp)
  o <- order(logp)
  adj <- logp[o] + log(choose(num.genes, 2)) - log(seq_len(n))
  adj <- rev(cummin(rev(adj)))
  adj.orig <- numeric(n); adj.orig[o] <- adj

  keep <- adj.orig < log(alpha)
  i <- p.matrix[keep, "row"]; j <- p.matrix[keep, "col"]

  gene.names <- rownames(mcpccs)

  mk <- function(x) Matrix::sparseMatrix(i = i, j = j, x = x,
                                         dims = c(num.genes, num.genes), symmetric = TRUE,
                                         dimnames = list(gene.names, gene.names))

  pccs <- mk(mcpccs[cbind(i, j)])
  logp <- mk(adj.orig[keep])
  print("Done.")

  list(pccs  = pccs,
       logp  = logp,
       alpha = alpha)

}
