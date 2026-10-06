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
