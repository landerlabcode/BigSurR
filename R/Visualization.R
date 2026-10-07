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
