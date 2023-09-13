Cumulants.PCC <- function(residuals.list, inv.sqrt.moments){

  n <- residuals.list$num.cells
  m <- residuals.list$ematrix
  c <- residuals.list$c

  options(matprod="default")

  f2 <- inv.sqrt.moments[[1]]
  f3 <- inv.sqrt.moments[[2]]
  f4 <- inv.sqrt.moments[[3]]
  f5 <- inv.sqrt.moments[[4]]

  k3.matrix <- (1+c^2*m*(3+c^2*(3+c^2)*m))/(sqrt(m)*(1+c^2*m)^(3/2))
  k4.matrix <- (1+m*(3+c^2*(7+m*(6+3*c^2*(6+m)+c^4*(6+(16+15*c^2+6*c^4+c^6)*m)))))/(m*(1+c^2*m)^2)
  k5.matrix.1 <- (1+c^2 * m * (3+c^2*(3+c^2)*m))/(sqrt(m)*(1+c^2*m)^(3/2))
  k5.matrix.2 <- 1/(m^(3/2)*(1+c^2*m)^(5/2)) * (1 + 5*(2+3*c^2)*m + 5*c^2*(8+15*c^2+5*c^4)*m^2
                                                +10*c^4*(6+17*c^2+15*c^4+6*c^6+c^8)*m^3+
                                                  c^6*(30+135*c^2+222*c^4+205*c^6+120*c^8+45*c^10+10*c^12+c^14)*m^4)
  k3.crossprod <- tcrossprod(k3.matrix)
  k4.crossprod <-  tcrossprod(k4.matrix)
  k5.crossprod.1 <- tcrossprod(k5.matrix.1)
  k5.crossprod.2  <- tcrossprod(k5.matrix.2)


  kappa2 <- 1/(n-1)^2 * f2 * n
  kappa3 <- 1/(n-1)^3 * f3 * k3.crossprod
  kappa4 <- 1/(n-1)^4 * (-3*n*f2^2 + f4 * k4.crossprod)
  kappa5 <- 1/(n-1)^5 * (-10 * f2 * f3 * k5.crossprod.1 + f5 * k5.crossprod.2)


  k.list <- list(kappa2, kappa3, kappa4, kappa5)

  return(k.list)
}


CF.Coefficients.PCC <- function(k.list, mcPCCs){

  k2 <- k.list[[1]]
  k3 <- k.list[[2]]
  k4 <- k.list[[3]]
  k5 <- k.list[[4]]

  c1 <- -mcPCCs-k3/(6*k2)+17*k3^3/(324*k2^4)-k3*k4/(12*k2^3)+k5/(40*k2^2)
  c2 <- sqrt(k2)+5*k3^2/(36*k2^(5/2))-k4/(8*k2^(3/2))
  c3 <- k3/(6*k2)-53*k3^3/(324*k2^4)+5*k3*k4/(24*k2^3)-k5/(20*k2^2)
  c4 <- -k3^2/(18*k2^(5/2))+k4/(24*k2^(3/2))
  c5 <- k3^3/(27*k2^4)-k3*k4/(24*k2^3)+k5/(120*k2^2)

  clist <- list(c1,c2,c3,c4,c5)

  mcPCCs.length <-nrow(clist[[1]])
  z <- sequence(mcPCCs.length)

  row <- unlist(lapply(2:mcPCCs.length, function(x) x:mcPCCs.length), use.names = FALSE)
  col <- rep(z[-length(z)], times = rev(tail(z, -1))-1)
  c1 <-  c(clist[[1]][lower.tri(clist[[1]])])
  c2 <- c(clist[[2]][lower.tri(clist[[2]])])
  c3 <- c(clist[[3]][lower.tri(clist[[3]])])
  c4 <- c(clist[[4]][lower.tri(clist[[4]])])
  c5 <- c(clist[[5]][lower.tri(clist[[5]])])

  cmatrix <-cbind(row, col, c1, c2, c3, c4, c5)

  return(cmatrix)
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

CF.PCC.pval <- function(roots.matrix){
  print("Estimating p-values.")

  p<- pnorm(abs(roots.matrix[,3]), log.p=T)
  p.matrix <- cbind(roots.matrix, p)

  p.mpfr <- apply(p.matrix,
                  1,
                  function(x){
                    ifelse(abs(x[3]) < 8.2,
                           -log10(1-exp(x[4])),
                            ifelse(abs(x[3]) >= 38.4,
                                   as.double(-log10((0.5 * exp(mpfr(-(x[3]^2) / 2, precBits=128))))),
                                   -log10(-x[4]/log(10)))
                            )
                  })

  p.matrix <- cbind(p.matrix, p.mpfr)

  print("P-value estimation complete.")
  return(p.matrix)
}
