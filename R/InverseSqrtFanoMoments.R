inv.sqrt.moment.interpolation2 <- function(correction, gene.totals) {
  moments.mat <- matrix(unlist(correction$moments), ncol = 4, byrow = TRUE)
  points <- correction$points

  # rule = 2 clamps at the endpoints instead of returning NA. Genes below
  # points[1] (a total of 1, say) would otherwise poison an entire row AND
  # column of every gene x gene moment matrix.
  int.moments <- lapply(1:4, function(k) {
    10^approx(log10(points), log10(moments.mat[, k]),
              xout = log10(gene.totals), rule = 2)$y
  })

  e.moments <- lapply(int.moments, function(v) v %o% v)
  attr(e.moments, "vectors") <- int.moments
  e.moments
}


InverseSqrtFanoMoments2 <- function(elist, c, n, trials) {

  samples <- rep(0, trials)
  x <- elist
  mu <- log(x / sqrt(1 + c^2))
  sigma <- sqrt(log(1 + c^2))


  for (i in 1:trials) {
    rate <- rlnorm(n, meanlog = mu, sdlog = sigma)
    pois.samples <- rpois(n, rate)
    samples[i] <- 1/sqrt(sum((pois.samples - x)^2/(x + c^2*x^2))/(n - 1))
  }
    results <- all.moments(samples, order.max = 5)
    results <- results[3:6]
    return(results)
  }

inv.sqrt.correction2 <- function(residuals.list, eta, theta){

  n <- residuals.list$num.cells
  if (n < 100) {
    stop("inverse-sqrt moment interpolation requires at least 100 cells; ",
         sprintf("got %d.", n))
  }
  points <- inv.sqrt.points(residuals.list$gene.totals, n)

  simemat <- outer(points, residuals.list$depthlist)

  c <- sqrt(eta * n / points + theta)

  #trials <- as.integer(4E7/(n*(log10(points)^(1/5)+0.5*log10(points)^3)))
  trials <- round(4e7/(n*(log10(points)^(1/5) + 0.5*log10(points)^3)))

  moments <- vector("list", length(points))

  for (i in seq_along(points)) {
    moments[[i]] <- InverseSqrtFanoMoments2(simemat[i, ], c[i], n, trials[i])
  }

  return(list(moments = moments, points = points, trials = trials))
}

inv.sqrt.points <- function(gene.totals, num.cells, faithful = TRUE) {
  lo  <- max(2, min(gene.totals))
  hi  <- max(gene.totals)
  mid <- num.cells / 50

  if (!(hi > lo)) {
    stop("all genes have the same total; cannot build an interpolation grid.")
  }

  if (!is.finite(mid)) stop("num.cells must be finite.")

  if (mid <= lo || mid >= hi) {
    if (faithful) {
      warning(sprintf(
        "breakpoint num.cells/50 = %g lies outside (%g, %g); the grid is non-monotonic, matching Mathematica.",
        mid, lo, hi))
    } else {
      mid <- sqrt(lo * hi)
    }
  }

  # 5 log-spaced points lo..mid, then 4 log-spaced mid..hi with mid dropped:
  # exactly Mathematica's {a, a(e/a)^(1/4), a(e/a)^(1/2), a(e/a)^(3/4), e,
  #                        e(h/e)^(1/3), e(h/e)^(2/3), h}
  c(exp(seq(log(lo),  log(mid), length.out = 5)),
    exp(seq(log(mid), log(hi),  length.out = 4))[-1])
}
