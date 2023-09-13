inv.sqrt.moment.interpolation2 <- function(sample.moments, gene.totals, points){

  moments.mat <- matrix(unlist(sample.moments), ncol=4, byrow=T)

  int.moments <- list(
    10^approx(log10(points), log10(moments.mat[,1]), xout= log10(gene.totals))$y,
    10^approx(log10(points), log10(moments.mat[,2]), xout= log10(gene.totals))$y,
    10^approx(log10(points), log10(moments.mat[,3]), xout= log10(gene.totals))$y,
    10^approx(log10(points), log10(moments.mat[,4]), xout= log10(gene.totals))$y
  )

  e.moments <- list(int.moments[[1]]%o%int.moments[[1]], int.moments[[2]]%o%int.moments[[2]],
                    int.moments[[3]]%o%int.moments[[3]], int.moments[[4]]%o%int.moments[[4]])
}

InverseSqrtFanoMoments2 <- function(elist, c, n, trials) {

  samples <- rep(0, trials)
  x <- elist
  mu <- log(x / sqrt(1 + c^2))
  sigma <- sqrt(log(1 + c^2))


  for (i in 1:trials) {
   pois.samples <- rep(0, n)
   for (y in 1:n){
     rate <- rlnorm(1, meanlog = mu[y], sdlog = sigma)
     pois.samples[y] <- rpois(1, rate)
   }

   sample <- 1/sqrt(sum((pois.samples-x)^2/(x+c^2*x^2))/(n-1))
   samples[i] <- sample
  }

  results <- all.moments(samples, order.max = 4)
  results <- results[2:5]
  return(results)
}

inv.sqrt.correction2 <- function(residuals.list, c){
  a <- max(2, min(residuals.list$gene.totals))
  e <- residuals.list$num.cells/50
  h <- max(residuals.list$gene.totals)
  n <- residuals.list$num.cells
  points <- as.integer(c(a, a*(e/a)^(1/4), a*(e/a)^(1/2), a*(e/a)^(3/4), e, e*(h/e)^(1/3), e*(h/e)^(2/3), h))

  simemat <- outer(points, residuals.list$depthlist)

  trials <- as.integer(4E7/(n*(log10(points)^(1/5)+0.5*log10(points)^3)))

  moments <- list()

  for(i in 1:8){
     moments[[i]] <- InverseSqrtFanoMoments2(simemat[i,], c, n, trials[i])
  }
  return(moments)
}
