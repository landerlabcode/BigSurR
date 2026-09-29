Cumulants.PCC <- function(residuals, inv.sqrt.moments){

  n <- residuals$num.cells
  m <- residuals$ematrix
  c <- residuals$c

  options(matprod="default")

  f2 <- inv.sqrt.moments[[1]]
  f3 <- inv.sqrt.moments[[2]]
  f4 <- inv.sqrt.moments[[3]]
  f5 <- inv.sqrt.moments[[4]]

  k3.matrix <- (1+c^2*m*(3+c^2*(3+c^2)*m))/(sqrt(m)*(1+c^2*m)^(3/2))
  k4.matrix <- (1+m*(3+c^2*(7+m*(6+3*c^2*(6+m)+c^4*(6+(16+15*c^2+6*c^4+c^6)*m)))))/(m*(1+c^2*m)^2)
  #k5.matrix.1 <- (1+c^2 * m * (3+c^2*(3+c^2)*m))/(sqrt(m)*(1+c^2*m)^(3/2))
  k5.matrix.2 <- 1/(m^(3/2)*(1+c^2*m)^(5/2)) * (1 + 5*(2+3*c^2)*m + 5*c^2*(8+15*c^2+5*c^4)*m^2
                                                +10*c^4*(6+17*c^2+15*c^4+6*c^6+c^8)*m^3+
                                                  c^6*(30+135*c^2+222*c^4+205*c^6+120*c^8+45*c^10+10*c^12+c^14)*m^4)
  k3.crossprod <- tcrossprod(k3.matrix)
  k4.crossprod <-  tcrossprod(k4.matrix)
  #k5.crossprod.1 <- tcrossprod(k5.matrix.1)
  k5.crossprod.2  <- tcrossprod(k5.matrix.2)


  kappa2 <- 1/(n-1)^2 * f2 * n
  kappa3 <- 1/(n-1)^3 * f3 * k3.crossprod
  kappa4 <- 1/(n-1)^4 * (-3*n*f2^2 + f4 * k4.crossprod)
  kappa5 <- 1/(n-1)^5 * (-10 * f2 * f3 * k3.crossprod + f5 * k5.crossprod.2)


  k.list <- list(kappa2, kappa3, kappa4, kappa5)

  return(k.list)
}


CF.Coefficients.PCC <- function(k.list, mcPCCs) {
  G  <- nrow(k.list[[1]])
  lt <- lower.tri(k.list[[1]])

  k2 <- k.list[[1]][lt]; k3 <- k.list[[2]][lt]
  k4 <- k.list[[3]][lt]; k5 <- k.list[[4]][lt]

  c1 <- -mcPCCs[lt] - k3/(6*k2) + 17*k3^3/(324*k2^4) - k3*k4/(12*k2^3) + k5/(40*k2^2)
  c2 <- sqrt(k2) + 5*k3^2/(36*k2^(5/2)) - k4/(8*k2^(3/2))
  c3 <- k3/(6*k2) - 53*k3^3/(324*k2^4) + 5*k3*k4/(24*k2^3) - k5/(20*k2^2)
  c4 <- -k3^2/(18*k2^(5/2)) + k4/(24*k2^(3/2))
  c5 <- k3^3/(27*k2^4) - k3*k4/(24*k2^3) + k5/(120*k2^2)

  ij <- which(lt, arr.ind = TRUE)
  cbind(row = ij[, 1], col = ij[, 2], c1, c2, c3, c4, c5)
}


QuickTest6CF <- function(cmatrix, first.pass.cutoff){

  cut <- sqrt(2)*erfcinv(2*10^-first.pass.cutoff)

  testfunc.1 <- function(x, c1, c2, c3, c4, c5){c1+c2*x+c3*x^2+c4*x^3+c5*x^4}
  testfunc.2 <- function(x, c1, c2, c3, c4, c5){c1*(c1+c2*x+c3*x^2+c4*x^3+c5*x^4)}

  a <- cmatrix[,3]
  b <- cmatrix[,4]
  c <- cmatrix[,5]
  d <- cmatrix[,6]
  e <- cmatrix[,7]

  cut.vec <- cbind(pos = testfunc.1(cut, a, b, c, d, e),
                   neg = testfunc.1(-cut, a, b, c, d, e))

  cut.bool <- apply(cut.vec, 1, function(x) ifelse(x[1]*x[2]<0, F, T))

  cmatrix <- cmatrix[which(cut.bool==T),]

  a <- cmatrix[,3]
  b <- cmatrix[,4]
  c <- cmatrix[,5]
  d <- cmatrix[,6]
  e <- cmatrix[,7]

  cut.vec2 <-testfunc.2(cut, a, b, c, d, e)

  cut.bool2 <- unlist(lapply(cut.vec2, function(x) ifelse(x<0, F, T)))

  cmatrix <- cmatrix[which(cut.bool2==T),]

  return(cmatrix)
}

SecondTestCF <- function(cmatrix.more.testing, first.pass.cutoff){
  cut <- sqrt(2)*erfcinv(2*10^-first.pass.cutoff)
  dfunc <- function(x, c2, c3, c4, c5){c2 + 2*c3*x + 3*c4*x^2 +4*c5*x^3}

  test.conditions <- function(x){
      ifelse(
        dfunc(-cut, x[4], x[5], x[6], x[7])<0,
        T,
        ifelse(
          dfunc(cut, x[4], x[5], x[6], x[7])<0,
          F,
          ifelse(
            3*x[6]^2 < 8*x[5]*x[7],
            T,
            ifelse(
              (-cut<(3*x[6]-sqrt(9*x[6]^2-24*x[5]*x[7]))/(12*x[7]))
              &((3*x[6]-sqrt(9*x[6]^2-24*x[5]*x[7]))/(12*x[7])<cut)
              & ((45*x[6]^3-36*x[5]*x[6]*x[7]-15*x[6]^2*sqrt(9*x[6]^2-24*x[5]*x[7])+8*x[7]*(9*x[4]*x[7]-x[5]*sqrt(9*x[6]^2-24*x[5]*x[7])))<0),
              F,
              ifelse(
                (-cut<(3*x[6]+sqrt(9*x[6]^2-24*x[5]*x[7]))/(12*x[7]))
                &((3*x[6]+sqrt(9*x[6]^2-24*x[5]*x[7]))/(12*x[7])<cut)
                & ((45*x[6]^3-36*x[5]*x[6]*x[7]+15*x[6]^2*sqrt(9*x[6]^2-24*x[5]*x[7])+8*x[7]*(9*x[4]*x[7]+x[5]*sqrt(9*x[6]^2-24*x[5]*x[7])))<0),
                F,
                T
              )
            )
          )
        )
      )
  }

  test.results.bool <- apply(cmatrix.more.testing, 1, test.conditions)

  cmatrix.passed <- cmatrix.more.testing[which(test.results.bool==T),]

  return(cmatrix.passed)
}

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

CF.PCC.Roots <- function(cmatrix, first.pass.cutoff, gene.totals){
  print("Beginning root finding process for Cornish Fisher.")

  cmatrix.pruned.1 <- QuickTest6CF(cmatrix, first.pass.cutoff)

  print(sprintf("First pruning complete. Removed %s insignificant correlations.", dim(cmatrix)[1]-dim(cmatrix.pruned.1)[1]))

  to.test.bool <- apply(cmatrix.pruned.1, 1,
                             function(x) ifelse((gene.totals[x[1]]<=84)|(gene.totals[x[2]]<=84),
                                                T, F))

  cmatrix.pruned.1 <- cbind(cmatrix.pruned.1, to.test.bool)

  cmatrix.pruned.2 <- cmatrix.pruned.1[cmatrix.pruned.1[,8]==F, ]

  cmatrix.more.testing <- cmatrix.pruned.1[cmatrix.pruned.1[,8]==T, ]

  cmatrix.passed <- SecondTestCF(cmatrix.more.testing, first.pass.cutoff)

  cmatrix.pruned.2 <- rbind(cmatrix.pruned.2, cmatrix.passed)

  print(sprintf("Second pruning complete. %s correlations remain.", dim(cmatrix.pruned.2)[1]))

  print("Beginning root finding.")

  roots <- apply(apply( apply(cmatrix.pruned.2[,c(3,4,5,6,7)], 1, polyroot) , 1, function(x) ifelse( abs(Im(x)) < 0.00001, Re(x),NA)),1, function(x) ifelse( !all(is.na(x)), min(abs(x), na.rm=T), NA))

  cmatrix.pruned.2 <- cbind(cmatrix.pruned.2, roots)

  found.roots <- cmatrix.pruned.2[which(!is.na(cmatrix.pruned.2[,9])), ]

  unfound.roots <- cmatrix.pruned.2[which(is.na(cmatrix.pruned.2[,9])), ]

  if(length(unfound.roots[,1])==0){

    roots.matrix <- found.roots[, c(1,2,9)]

  }else{if(length(unfound.roots[,1])==1){

    single.d.root <- min(polyroot(c(unfound.roots[,4],
                                    2*unfound.roots[,5], 3*unfound.roots[,6], 4*unfound.roots[,7])), na.rm=T)

    unfound.roots <- cbind(unfound.roots, single.d.root)

    roots.matrix <- rbind(found.roots[, c(1,2,9)], unfound.roots[, c(1,2,10)])

  }else{

    d.coefficients <- cbind(unfound.roots[,4], 2*unfound.roots[,5], 3*unfound.roots[,6], 4*unfound.roots[,7])

    d.roots <- apply(apply( apply(d.coefficients, 1, polyroot) , 1, function(x) ifelse( abs(Im(x)) < 0.00001, Re(x),NA)),1, function(x) ifelse( !all(is.na(x)), min(x, na.rm=T), NA))

    unfound.roots <- cbind(unfound.roots, d.roots)

    roots.matrix <- rbind(found.roots[, c(1,2,9)], unfound.roots[, c(1,2,10)])

  }
  }

  print("Root finding complete.")
  return(roots.matrix)
}

CF.PCC.blocked <- function(residuals, inv.sqrt.moments, mcPCC,
                           first.pass.cutoff, block = 500, npts = 9) {
  n <- residuals$num.cells
  G <- nrow(residuals$ematrix)
  cut <- sqrt(2) * erfcinv(2 * 10^-first.pass.cutoff)

  xs  <- seq(-cut, cut, length.out = npts)
  tab <- t(outer(xs, 0:4, "^"))          # 5 x npts, built once

  m <- residuals$ematrix; cc <- residuals$c
  K3  <- (1+cc^2*m*(3+cc^2*(3+cc^2)*m))/(sqrt(m)*(1+cc^2*m)^(3/2))
  K4  <- (1+m*(3+cc^2*(7+m*(6+3*cc^2*(6+m)+cc^4*(6+(16+15*cc^2+6*cc^4+cc^6)*m)))))/(m*(1+cc^2*m)^2)
  K52 <- 1/(m^(3/2)*(1+cc^2*m)^(5/2)) * (1 + 5*(2+3*cc^2)*m + 5*cc^2*(8+15*cc^2+5*cc^4)*m^2 +
                                           10*cc^4*(6+17*cc^2+15*cc^4+6*cc^6+cc^8)*m^3 +
                                           cc^6*(30+135*cc^2+222*cc^4+205*cc^6+120*cc^8+45*cc^10+10*cc^12+cc^14)*m^4)
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
