get.residuals <- function(seuratob, assay, counts.slot, cv.est.method){

  DefaultAssay(object = seuratob) <- assay
  raw.counts <- as.matrix(GetAssayData(object=seuratob, layer=counts.slot))

  if(!(all(raw.counts==floor(raw.counts)))){
    stop("Non-integer counts data detected. Please supply unnormalized data.")
  }
  gene.totals <- rowSums(raw.counts)

  if(0 %in% gene.totals){
    stop("Genes with zero counts were found, run quality control steps before attempting further analysis.")
  }

  num.cells <- ncol(raw.counts)

  num.genes <- nrow(raw.counts)


  cell.total.umis <- colSums(raw.counts)

  if (any(cell.total.umis <= 0)) {
    stop(sprintf("%d cells have zero counts; filter them before fitting.",
                 sum(cell.total.umis <= 0)))
  }

  all.umis<- sum(cell.total.umis)
  depthlist <- cell.total.umis/all.umis

  ematrix <- gene.totals %o% depthlist

  eta_theta <- find.best.c(ematrix, depthlist, gene.totals, num.cells, raw.counts, cv.est.method)
  c <-sqrt((eta_theta$eta)/(gene.totals/num.cells)+eta_theta$theta)
  residuals <- (raw.counts-ematrix)/sqrt(ematrix*(1+c^2*ematrix))

  data.list <- list(residuals, ematrix, depthlist, num.cells, num.genes, gene.totals, c, eta_theta$eta, eta_theta$theta, eta_theta$j)
  data.list <- setNames(data.list, c("residuals","ematrix","depthlist","num.cells","num.genes", "gene.totals", "c","eta","theta", "j_values"))
  return(data.list)
}


find.best.c <- function(ematrix, depthlist, gene.totals, num.cells, raw.counts, cv.est.method){
  #Calculate the single value of c for all genes
  if(!(cv.est.method %in% c("Single","MeanSpecific","TwoComponent"))){
    stop("Invalid option for parameter cv.est.method supplied. Please choose one of 'Single','MeanSpecific', or 'TwoComponent'.")
  }
  if(cv.est.method=="Single"){
    test.cs <- seq(from= 0, to= 1, by = 0.05)

    best.slope <- c("index"=0,"slope"=2000)

    for(i in 1:length(test.cs)){

      c <- test.cs[i]

      X <- log10(rowMeans(raw.counts))

      Y <- log10(1/(num.cells-1) * rowSums(((raw.counts-ematrix)/sqrt(ematrix*(1+c^2*ematrix)))^2))

      fit.indexes <- which((-1 < X) & (X < 2))

      X <- X[fit.indexes]
      Y <- Y[fit.indexes]

      test.fit <- lm(Y~X)

      coefficients <- summary(test.fit)$coefficients

      slope <- coefficients[2,1]

      if(abs(slope-0) < abs(best.slope[2]-0)){
        best.slope <- c("index"=i, "slope"=slope)
      }

      if(slope < 0){
        break
      }
    }
    return(list(eta=0,theta=test.cs[best.slope[["index"]]]^2))
    }
  #Calculate the mean-specific c curve and fit genes
  else if(cv.est.method=="TwoComponent"){
    fit<-fit.eta.theta(ematrix, raw.counts, gene.totals, depthlist, num.cells)
    cs <- list(eta=fit$eta, theta=fit$theta, j=fit$j)
    return(cs)
    }
  else if(cv.est.method=="MeanSpecific"){
    fit<-fit.theta.only(ematrix, raw.counts, gene.totals, depthlist, num.cells)
    cs <- list(eta=fit$eta, theta=fit$theta, j=fit$j)
    return(cs)
    }
  }

fit.js <- function(D, E, target, tol = 1e-11, maxit = 60) {
  m <- nrow(D)
  g  <- function(j) rowSums(D / (E * (1 + j * E)))
  dg <- function(j) -rowSums(D / (1 + j * E)^2)

  emax <- apply(E, 1, max)
  lo <- -1 / emax * (1 - 1e-9)         # g -> +Inf as j approaches this
  hi <- rep(1, m)

  # push hi out until g(hi) < target
  for (k in 1:80) {
    too.high <- g(hi) > target
    if (!any(too.high)) break
    hi[too.high] <- hi[too.high] * 2
  }
  if (any(g(hi) > target)) warning("some genes never bracketed; returning NA")

  j <- pmin(pmax(0, lo + 1e-8), hi)    # start near 0, inside the bracket
  for (k in seq_len(maxit)) {
    gj <- g(j) - target
    ok <- is.finite(gj)
    lo[ok & gj > 0] <- j[ok & gj > 0]  # g too big -> root is to the right
    hi[ok & gj < 0] <- j[ok & gj < 0]

    step <- gj / dg(j)
    jnew <- j - step
    # fall back to bisection wherever Newton leaves the bracket
    bad <- !is.finite(jnew) | jnew <= lo | jnew >= hi
    jnew[bad] <- 0.5 * (lo[bad] + hi[bad])

    if (max(abs(jnew - j), na.rm = TRUE) < tol * max(1, max(abs(j)))) {
      j <- jnew
      break
    }
    j <- jnew
  }
  j[g(hi) > target] <- NA_real_
  j
}


fit.eta.theta<-function(ematrix, raw.counts, gene.totals, depthlist, num.cells, start = c(0.0, 0.1), bin.width=0.5){
  j.fit <- calc.lx.ly(ematrix, raw.counts, gene.totals, depthlist, num.cells)
  lx <- j.fit$lx
  ly <- j.fit$ly
  bin <- floor(lx / bin.width)
  w   <- 1 / sqrt(table(bin)[as.character(bin)])
  obj <- function(p) {
    if (any(p < 0)) return(1e10)
    sum(w*abs(ly - log1p(p[1] + p[2] * exp(lx))))
  }
  o <- optim(start, obj, method = "Nelder-Mead",
             control = list(maxit = 2000, reltol = 1e-12))
  list(eta = o$par[1], theta = o$par[2], objective = o$value,
       convergence = o$convergence, j = j.fit$j, means = j.fit$means, used = j.fit$used)
}

fit.theta.only<-function(ematrix, raw.counts, gene.totals, depthlist, num.cells, upper = 10, bin.width=0.5){
  j.fit <- calc.lx.ly(ematrix, raw.counts, gene.totals, depthlist, num.cells)
  lx <- j.fit$lx
  ly <- j.fit$ly
  bin <- floor(lx / bin.width)
  w   <- 1 / sqrt(table(bin)[as.character(bin)])
  obj <- function(theta) {
    sum(w*abs(ly - log1p(theta*exp(lx))))
  }
  o <- optimize(obj, interval=c(0, upper), tol=1e-12)
  list(eta=0, theta = o$minimum, objective = o$objective,
       convergence = 0L, j = j.fit$j, means = j.fit$means, used = j.fit$used)
}

calc.lx.ly<-function(ematrix, raw.counts, gene.totals, depthlist, num.cells){
  means <- gene.totals/num.cells
  N <- (raw.counts-ematrix)^2
  j <- fit.js(N, ematrix, target=num.cells-1)
  y  <- 1 + j * means
  keep <- is.finite(j) & is.finite(y) & y > 0 & means > 0
  lx <- log(means[keep])
  ly <- log(y[keep])
  out <- list(j=j, means=means, keep=keep, used=sum(keep),lx=lx, ly=ly)
  return(out)
}
