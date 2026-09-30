CF.PCC.Roots2 <- function(cmatrix, first.pass.cutoff) {
  p.at <- function(k, x) k[,1] + k[,2]*x + k[,3]*x^2 + k[,4]*x^3 + k[,5]*x^4

  cf <- cmatrix[, 3:7, drop = FALSE]

  roots <- vapply(seq_len(nrow(cf)), function(i) cf.root.first.crossing(cf[i, ]), numeric(1))

  found <- !is.na(roots)

  # pairs with no crossing: fall back to the stationary point of the quartic,
  # i.e. the smallest real root of its derivative
  if (any(!found)) {
    d <- cbind(cf[!found, 2], 2*cf[!found, 3], 3*cf[!found, 4], 4*cf[!found, 5])
    roots[!found] <- vapply(seq_len(nrow(d)), function(i) {
      r <- polyroot(d[i, ])
      r <- Re(r[abs(Im(r)) < 1e-5])
      if (length(r) == 0) NA_real_ else r[which.min(abs(r))]
    }, numeric(1))
  }

  cbind(cmatrix[, 1:2, drop = FALSE], root = roots)
}

cf.root.first.crossing <- function(k) {
  r <- polyroot(k); r <- Re(r[abs(Im(r)) < 1e-5])
  r <- r[sign(r) == -sign(k[1])]
  if (length(r)) r[which.min(abs(r))] else NA_real_
}

CF.PCC.blocked <- function(residuals, inv.sqrt.moments, mcPCC,
                           first.pass.cutoff, block = 500, npts = 9, null.dist=c("PLN","NB")){
  null.dist<- match.arg(null.dist)
  n <- residuals$num.cells
  G <- nrow(residuals$ematrix)
  cut <- sqrt(2) * erfcinv(2 * 10^-first.pass.cutoff)

  xs  <- seq(-cut, cut, length.out = npts)
  tab <- t(outer(xs, 0:4, "^"))          # 5 x npts, built once

  m <- residuals$ematrix; cc <- residuals$c

  if(null.dist=="PLN"){
  K3  <- (1+cc^2*m*(3+cc^2*(3+cc^2)*m))/(sqrt(m)*(1+cc^2*m)^(3/2))
  K4  <- (1+m*(3+cc^2*(7+m*(6+3*cc^2*(6+m)+cc^4*(6+(16+15*cc^2+6*cc^4+cc^6)*m)))))/(m*(1+cc^2*m)^2)
  K52 <- 1/(m^(3/2)*(1+cc^2*m)^(5/2)) * (1 + 5*(2+3*cc^2)*m + 5*cc^2*(8+15*cc^2+5*cc^4)*m^2 +
                                           10*cc^4*(6+17*cc^2+15*cc^4+6*cc^6+cc^8)*m^3 +
                                           cc^6*(30+135*cc^2+222*cc^4+205*cc^6+120*cc^8+45*cc^10+10*cc^12+cc^14)*m^4)
  } else{
    xi <- cc^2
    K3   <- (1 + 2*xi*m) / sqrt(m*(1 + xi*m))
    K4   <- (1 + 3*(1 + 2*xi)*m*(1 + xi*m)) / (m*(1 + xi*m))
    K52 <- ((1 + 2*xi*m) * (1 + 2*(5 + 6*xi)*m*(1 + xi*m)))/(m*(1 + xi*m))^(3/2)
  }

  v <- attr(inv.sqrt.moments, "vectors")
  f2 <- v[[1]]; f3 <- v[[2]]
  f4 <- v[[3]]; f5 <- v[[4]]

  out <- vector("list", ceiling(G / block))

  for (b in seq_along(out)) {
    I <- ((b-1)*block + 1):min(b*block, G)

    g3  <- tcrossprod(K3[I, , drop = FALSE], K3)
    g4  <- tcrossprod(K4[I, , drop = FALSE], K4)
    g52 <- tcrossprod(K52[I, , drop = FALSE], K52)

    F2 <- outer(f2[I], f2); F3 <- outer(f3[I], f3)
    F4 <- outer(f4[I], f4); F5 <- outer(f5[I], f5)

    k2 <- F2 * n / (n-1)^2
    k3 <- F3 * g3 / (n-1)^3
    k4 <- (-3*n*F2^2 + F4*g4) / (n-1)^4
    k5 <- (-10*F2*F3*g3 + F5*g52) / (n-1)^5

    keep <- which(outer(I, 1:G, ">"))
    k2 <- k2[keep]; k3 <- k3[keep]; k4 <- k4[keep]; k5 <- k5[keep]
    pcc <- mcPCC[I, , drop = FALSE][keep]

    cf <- cbind(
      -pcc - k3/(6*k2) + 17*k3^3/(324*k2^4) - k3*k4/(12*k2^3) + k5/(40*k2^2),
      sqrt(k2) + 5*k3^2/(36*k2^(5/2)) - k4/(8*k2^(3/2)),
      k3/(6*k2) - 53*k3^3/(324*k2^4) + 5*k3*k4/(24*k2^3) - k5/(20*k2^2),
      -k3^2/(18*k2^(5/2)) + k4/(24*k2^(3/2)),
      k3^3/(27*k2^4) - k3*k4/(24*k2^3) + k5/(120*k2^2))

    surv <- which(abs(rowSums(sign(cf %*% tab))) == npts)

    if (length(surv)) {
      ij <- arrayInd(keep[surv], c(length(I), G))
      out[[b]] <- cbind(row = I[ij[,1]], col = ij[,2], cf[surv, , drop = FALSE])
    }
  }
  do.call(rbind, out)
}

CF.PCC.pval <- function(roots.matrix) {
  z <- roots.matrix[, 3]
  logp <- pnorm(-abs(z), log.p = TRUE) + log(2)
  cbind(roots.matrix, logp = logp, direction = sign(z))
}
