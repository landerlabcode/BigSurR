#' SubsetCorrelationMatrix
#'
#' @param seuratObject Seurat object containing correlations matrix produced by BigSur.
#' @param pCutoff Boolean, Float. p-value threshold for correlations. Correlations with p-values larger than pCutoff will be set to 0.
#' @param minPosCorr Boolean, Float. Minimum positive correlation cutoff. Positive correlations with coefficients under minPosCorr will be set to 0.
#' @param minNegCorr Boolean, Float. Minimum negative correlation cutoff. Negative correlations with absolute values under abs(minNegCorr) will be set to 0.
#' @param genes Boolean, Character Vector. The set of genes for the correlation matrix to be subset to.
#'
#' @returns Sparse matrix (dsCMatrix) containing correlations meeting subset criteria.
#' @export
#'
SubsetCorrelationMatrix <- function(
    seuratObject,
    pCutoff=F,
    minPosCorr=F,
    minNegCorr=F,
    genes=F
){
  if(!(is.character(genes) || identical(genes, FALSE))){
    stop("Genes were not supplied as a vector of strings.")
  }
  if(!(is.numeric(pCutoff) || identical(pCutoff, FALSE))){
    stop("p-value cutoff not supplied as a number.")
  }
  if(!(is.numeric(minPosCorr) || identical(minPosCorr, FALSE))){
    stop("Minimum positive correlation cutoff not supplied as a number.")
  }
  if(!(is.numeric(minNegCorr) || identical(minNegCorr, FALSE))){
    stop("Minimum negative correlation cutoff not supplied as a number.")
  }

  corrs <- seuratObject@misc$BigSur.Correlations
  if (is.null(corrs)){
    stop("No BigSur correlations in this object.")
  }
  ps <- seuratObject@misc$BigSur.log.adj.pvalues
  if (is.numeric(pCutoff) && is.null(ps)){
    stop("No BigSur p-values in this object.")
  }

  if (is.numeric(pCutoff)) {
    stopifnot(identical(corrs@i, ps@i), identical(corrs@p, ps@p))
    corrs@x[ps@x > log(pCutoff)] <- 0
    corrs <- Matrix::drop0(corrs)
  }

  if (is.numeric(minPosCorr)) {
    corrs@x[corrs@x > 0 & corrs@x < minPosCorr] <- 0
    corrs <- Matrix::drop0(corrs)
  }

  if (is.numeric(minNegCorr)) {
    corrs@x[corrs@x < 0 & abs(corrs@x) < abs(minNegCorr)] <- 0
    corrs <- Matrix::drop0(corrs)
  }

  if(is.character(genes)){
    missing <-setdiff(genes, rownames(corrs))
    if(length(missing)){
      warning("Some requested genes not present in matrix: ", paste(missing, collapse=", "))
    }
    p.genes <- intersect(genes, rownames(corrs))
    corrs <- corrs[p.genes, p.genes, drop=F]
  }
  return(corrs)
}

#' FindCorrelationModules
#'
#' @param corr.matrix Sparse matrix containing the correlations of interest.
#' @param num.steps Integer. Number of steps used for walktrap clustering. (Default = 4)
#' @param weighted.clustering Boolean. Determines whether or not the walktrap will be performed on the correlations matrix or the adjacency matrix. (Default = FALSE)
#'
#' @returns igraph communities object.
#' @export
#'
FindCorrelationModules <- function(
    corr.matrix,
    num.steps = 4,
    weighted.clustering = FALSE
){
  #Lower triangularize and remove negative correlations
  lt.corr <- Matrix::tril(corr.matrix, k = -1)
  lt.corr@x[lt.corr@x<0] <- 0
  lt.corr <- Matrix::drop0(lt.corr)

  #Handle weighted and unweighted specification
  if(!weighted.clustering){
    lt.corr@x[]<-1
  }
  weights <- if (weighted.clustering) TRUE else NULL

  #Create adjacency graph for igraph
  adj.graph <- graph_from_adjacency_matrix(
    lt.corr,
    mode="lower",
    diag = FALSE,
    weighted = weights)
  adj.graph <- delete_vertices(adj.graph, degree(adj.graph)==0)
  #Find gene modules
  modules <- cluster_walktrap(
    adj.graph,
    steps = num.steps
  )
  return(modules)
}

MakeAdjGraph <- function(corr.matrix, keep.negative=TRUE){
  lt.corr <- Matrix::tril(corr.matrix, k = -1)
  if (!keep.negative) {
    lt.corr@x[lt.corr@x < 0] <- 0
  }
  lt.corr <- Matrix::drop0(lt.corr)

  adj.graph <- graph_from_adjacency_matrix(
    lt.corr,
    mode="lower",
    diag = FALSE,
    weighted = TRUE)

  return(adj.graph)
}


#' MergeSmallModules
#'
#' @param corr.matrix Sparse matrix containing the correlations of interest.
#' @param modules igraph communities object produced by FindCorrelationModules.
#' @param min.size Integer. Minimum size of output modules (if possible).
#'
#' @returns igraph communities object containing the merged gene modules.
#' @export
#'
MergeSmallModules <- function(corr.matrix, modules, min.size = 15) {
  graph <- MakeAdjGraph(corr.matrix, keep.negative = FALSE)
  graph <- delete_vertices(graph, degree(graph) == 0)
  mem <- membership(modules)
  stopifnot(identical(V(graph)$name, names(mem)))

  # Weighted adjacency if available, otherwise binary
  w.attr <- if ("weight" %in% edge_attr_names(graph)) "weight" else NULL
  A <- as_adjacency_matrix(graph, attr = w.attr, sparse = TRUE)

  comm <- as.integer(factor(mem))   # labels 1..k
  frozen <- integer(0)

  repeat {
    sizes <- tabulate(comm)
    small <- setdiff(which(sizes < min.size), frozen)
    if (length(small) == 0) break

    c0 <- small[which.min(sizes[small])]   # smallest first

    # Total edge weight between c0 and every community
    M <- Matrix::sparseMatrix(i = seq_along(comm), j = comm, x = 1,
                              dims = c(length(comm), length(sizes)))
    strength <- as.numeric(Matrix::crossprod(M, A %*% M[, c0]))
    strength[c0] <- -Inf

    best <- which.max(strength)
    if (strength[best] <= 0) {            # isolated: leave alone
      frozen <- c(frozen, c0)
      next
    }

    comm[comm == c0] <- best                # merge c0 into best
    comm[comm > c0] <- comm[comm > c0] - 1  # close the label gap
    frozen[frozen > c0] <- frozen[frozen > c0] - 1
  }

  names(comm) <- names(mem)
  make_clusters(graph, membership = comm, algorithm = "walktrap (merged)")
}

EigenvectorCentralityScores<- function(
  corr.matrix,
  pos.cutoff=0,
  neg.cutoff=0
){


}

ModuleMembershipScore<- function(
  bigsur.residuals,
  correlation.modules){

}

#' Plot positive and negative correlation density within and between communities
#'
#' @param corr.matrix Sparse gene x gene matrix of significant correlations
#'   (nonzero = significant), with gene names as dimnames. Symmetric or one
#'   triangle only; both work.
#' @param communities An igraph communities object or a list of gene vectors.
#' @param span Which communities to show (indices). Defaults to all.
#' @param max.pos.fraction,max.neg.fraction Fraction that maps to full opacity.
#' @param pos.color,neg.color Disk colors.
#' @param max.size Disk size for the largest block.
#' @export
InterModuleCorrelations <- function(corr.matrix,
                                       communities,
                                       span = NULL,
                                       max.pos.fraction = 0.02,
                                       max.neg.fraction = 0.01,
                                       pos.color = "forestgreen",
                                       neg.color = "#DC3220",
                                       max.size = 10) {
  if (inherits(communities, "communities")) {
    communities <- igraph::communities(communities)
  }
  if (is.null(span)) span <- seq_along(communities)
  coms <- communities[span]
  k <- length(coms)

  # Full matrix with both triangles and an empty diagonal
  # (adding the transpose fills in a missing triangle; only signs are used below)
  M <- methods::as(methods::as(corr.matrix, "CsparseMatrix"), "generalMatrix")
  M <- M + Matrix::t(M)
  Matrix::diag(M) <- 0
  M <- Matrix::drop0(M)

  # Gene indices for each community, and a gene x community indicator matrix
  idx <- lapply(coms, function(g) {
    i <- match(g, rownames(M))
    i[!is.na(i)]
  })
  Z <- Matrix::sparseMatrix(i = unlist(idx),
                            j = rep(seq_len(k), lengths(idx)),
                            x = 1, dims = c(nrow(M), k))

  # Number of edges between each pair of communities (each pair counted once)
  count.blocks <- function(A) {
    C <- as.matrix(Matrix::crossprod(Z, A %*% Z))
    diag(C) <- diag(C) / 2
    C
  }
  pos <- count.blocks((M > 0) * 1)
  neg <- count.blocks((M < 0) * 1)

  # Number of possible gene pairs in each block
  n <- lengths(idx)
  maxlen <- outer(n, n)
  diag(maxlen) <- choose(n, 2)

  # Long data frame for the lower triangle
  grid <- expand.grid(i = seq_len(k), j = seq_len(k))
  grid <- grid[grid$i >= grid$j, ]
  ij <- cbind(grid$i, grid$j)
  grid$maxlen <- maxlen[ij]
  grid$row.lab <- factor(span[grid$i], levels = rev(span))
  grid$col.lab <- factor(span[grid$j], levels = span)

  make.panel <- function(counts, max.frac, label) {
    d <- grid
    d$frac <- ifelse(d$maxlen > 0, counts[ij] / d$maxlen, 0)
    d$alpha <- pmin(d$frac / max.frac, 1)
    d$sign <- label
    d
  }
  df <- rbind(make.panel(pos, max.pos.fraction, "Positive"),
              make.panel(neg, max.neg.fraction, "Negative"))
  df$sign <- factor(df$sign, levels = c("Positive", "Negative"))

  ggplot2::ggplot(df, ggplot2::aes(x = col.lab, y = row.lab)) +
    ggplot2::geom_tile(fill = "grey97", color = "white") +
    # faint outline showing the disk size even when nothing is observed
    ggplot2::geom_point(ggplot2::aes(size = maxlen), shape = 1,
                        color = "grey85", stroke = 0.3) +
    ggplot2::geom_point(ggplot2::aes(size = maxlen, alpha = alpha, color = sign)) +
    ggplot2::scale_size_area(max_size = max.size, guide = "none") +
    ggplot2::scale_alpha_identity() +
    ggplot2::scale_color_manual(values = c(Positive = pos.color,
                                           Negative = neg.color),
                                guide = "none") +
    ggplot2::facet_wrap(~ sign) +
    ggplot2::coord_fixed() +
    ggplot2::labs(x = NULL, y = NULL,
                  caption = sprintf("Full opacity = observed fraction of %g (positive), %g (negative)",
                                    max.pos.fraction, max.neg.fraction)) +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid = ggplot2::element_blank())
}
