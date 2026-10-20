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

  new.mods <- make_clusters(graph, membership = comm, algorithm = "walktrap (merged)")
  new.mods$names <- V(graph)$name

  return(new.mods)
}

eigen_centrality_by_component <- function(g) {
  comp <- igraph::components(g)
  out  <- setNames(numeric(igraph::vcount(g)), igraph::V(g)$name)
  for (k in seq_len(comp$no)) {
    idx <- which(comp$membership == k)
    if (length(idx) == 1) next                      # isolated gene -> weight 0
    v <- igraph::eigen_centrality(igraph::induced_subgraph(g, idx))$vector
    out[idx] <- v / sum(v) * (length(idx) - 1)
  }
  if (sum(out) > 0) out / sum(out) else out
}

EVCentralityScore <- function(corr.mat, poscutoff = 0, negcutoff = 0) {
  ev <- function(adj) {
    g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", diag = FALSE)
    eigen_centrality_by_component(g)
  }
  sort(sqrt(ev((corr.mat > poscutoff) * 1) * ev((corr.mat < negcutoff) * 1)),decreasing=TRUE)
}

MembershipScore<- function(
    bigsur.residuals,
    corr.mat,
    genes,
    cutoff=0){

  #Get gene weights for present genes
  genes <- genes[genes %in% rownames(bigsur.residuals)]
  adj <- as.matrix(corr.mat[genes,genes,drop=FALSE]) > cutoff
  w <- rowSums(adj)

  #Standardize genes across cells
  z <- as.matrix(bigsur.residuals[genes, , drop=FALSE])
  z <- t(scale(t(z)))

  #Drop zero-variance genes
  pass <- rowSums(!is.finite(z))==0
  if(any(!pass)){
    z <- z[pass, , drop=FALSE]
    w <- w[pass]
  }
  #Return per-cell score
  colSums(z * w)/sum(w)
}

ModuleMembershipScores <- function(seurat.obj, modules, which=NULL, cutoff=0){
  member <- if(inherits(modules, "communities")) igraph::communities(modules) else modules
  if(is.null(which)) which <- seq_along(member)
  if (!is.character(member[[1]])) {
    stop("Community members aren't gene names; build the graph with vertex names set.")
  }
  bigsur.residuals <- seurat.obj@assays$RNA$data
  corr.mat <- seurat.obj@misc$BigSur.Correlations

  scores <- vapply(member[which],
                   function(g) MembershipScore(bigsur.residuals, corr.mat, g, cutoff),
                   numeric(ncol(bigsur.residuals)))
  scores <- matrix(scores, ncol=length(which),
                   dimnames=list(colnames(bigsur.residuals), paste0("Module_", which)))
  AddMetaData(seurat.obj, metadata = as.data.frame(scores))
}
