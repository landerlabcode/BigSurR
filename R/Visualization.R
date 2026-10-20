#' StaticCorrelationPlot
#'
#' Plots gene-gene correlations as a network graph. Positive correlations are shown in blue, negative in red. The length of the positive edges indicates the strength of that correlation.
#'
#' @param corr.matrix Sparse matrix containing the correlations of interest.
#' @param highlight Vector of gene symbols to be highlighted. If supplied, the relevant gene symbols will be written in red.
#' @param modules igraph communities object. If supplied, the nodes will change color and shape to indicate the module they belong to.
#'
#' @returns ggraph plot of correlations.
#' @export
#'
StaticCorrelationPlot <- function(
    corr.matrix,
    highlight = NULL,
    modules = NULL
    ){
  #Create underlying graph embedding
  adj.graph <- MakeAdjGraph(corr.matrix)
  adj.graph <- delete_vertices(adj.graph, degree(adj.graph) == 0)
  lay <- LayoutFromPos(adj.graph)

  #Default node appearance
  node.fill  <- rep("#E1A73B", vcount(adj.graph))
  node.shape <- rep(21, vcount(adj.graph))

  #Change node color and shape if modules are supplied
  if (!is.null(modules)) {
    mem <- membership(modules)
    size.rank <- rank(-sizes(modules), ties.method = "first")
    m <- size.rank[as.integer(mem[V(adj.graph)$name])]  # NA = not in a module

    n.col  <- 8
    pal    <- hcl.colors(n.col, "Dark 3")
    shapes <- c(21, 24, 22, 23, 25)

    node.fill  <- ifelse(is.na(m), "grey85", pal[(m - 1) %% n.col + 1])
    node.shape <- ifelse(is.na(m), 21, shapes[((m - 1) %/% n.col) %% length(shapes) + 1])
  }

  V(adj.graph)$fill <-node.fill
  V(adj.graph)$shape <- node.shape

  #Mark genes to highlight
  V(adj.graph)$hl <- V(adj.graph)$name %in% highlight

  ggraph(adj.graph, layout = "manual", x = lay[, 1], y = lay[, 2]) +
    #geom_edge_link(aes(color = weight > 0), width = 0.6, alpha = 0.7) +
    geom_edge_link(aes(color = weight > 0, alpha = abs(weight)), width = 0.6) +
    scale_edge_alpha(range = c(0.05, 0.6), guide = "none")+
    scale_edge_color_manual(values = c(`TRUE` = "#005AB5", `FALSE` = "#E8756A"),
                            guide = "none") +
    geom_node_point(aes(fill = fill, shape = shape), size = 3,
                    color = "grey20", stroke = 0.3) +
    scale_fill_identity() +
    scale_shape_identity() +
    geom_node_text(aes(label = name, color = hl), repel = TRUE, size = 2.5,
                   fontface = "bold", max.overlaps = Inf) +
    scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "black"),
                       guide = "none") +
    theme_void()

}


LayoutFromPos <- function(adj.graph){
  g.pos <- delete_edges(adj.graph, E(adj.graph)[weight < 0])
  layout <- layout_components(g.pos, layout = layout_with_fr)
  return(layout)
}

FeaturePlotMembershipScore<- function(seurat.obj,
                                      modules,
                                      color.scaling=TRUE,
                                      pt.size=NULL,
                                      reduction=NULL,
                                      ncol=NULL,
                                      mid.color = "white",
                                      title=NULL){
  plots <- lapply(modules, function(f){
    v <- FetchData(seurat.obj, vars=f)[[1]]
    lim <- if(color.scaling) max(abs(v), na.rm=TRUE) else 1

    p <- suppressMessages(
      FeaturePlot(seurat.obj, features = f, pt.size = pt.size, reduction = reduction) +
        ggplot2::scale_color_gradient2(low="blue", mid=mid.color, high="red",
                                       midpoint = 0, limits = c(-lim, lim),
                                       oob = scales::squish, na.value="grey90")
    )
    if(!is.null(p$data[[f]])) p$data <- p$data[order(abs(p$data[[f]])), ]
    p
  })
  out <- wrap_plots(plots, ncol=ncol)
  if(!is.null(title)) out <- out + plot_annotation(title=title)
  out
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



