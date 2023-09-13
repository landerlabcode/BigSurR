Cumulants.Fano <- function(residuals.list, c){

  chi <- 1+c^2
  n <- residuals.list$num.cells
  m <- residuals.list$ematrix

  kappa1 <- rep(n / (-1 + n), dim(m)[1])

  kappa2 <- rowSums((-1+((1+m*(-4 + 7*chi + 6*(1 - 2*chi + chi^3)*m + (-3 + 6*chi - 4*chi^3 + chi^6)*m^2))) / (m*(1 + (-1 + chi)*m)^2))) / (-1 + n)^2

  kappa3 <- rowSums(1 / (m^2 * (1 + (-1 + chi)*m)^3) * (1 + (-9 + 31*chi)*m + 2*(16 - 57*chi + 45*chi^3)*m^2 + (-56 + 180*chi - 21*chi^2 - 168*chi^3 + 65*chi^6)*m^3 + 3*(16 - 48*chi + 14*chi^2 + 40*chi^3 - 6*chi^4 - 21*chi^6 + 5*chi^10)*m^4 + (-16 + 48*chi - 24*chi^2 - 30*chi^3 + 12*chi^4 + 18*chi^6 - 3*chi^7 - 6*chi^10 + chi^15)*m^5)) / (-1 + n)^3

  kappa4 <- rowSums(1 / (m^3 * (1 + (-1 + chi)*m)^4) * (1 + (-15 + 127*chi)*m + (92 - 674*chi + 966*chi^3)*m^2 + (-302 + 1724*chi - 271*chi^2 - 2804*chi^3 + 1701*chi^6)*m^3 + 6*(96 - 452*chi + 174*chi^2 + 620*chi^3 - 102*chi^4 - 511*chi^6 + 175*chi^10)*m^4 + 2*(-320 + 1344*chi - 822*chi^2 - 1390*chi^3 + 672*chi^4 + 1124*chi^6 - 151*chi^7 - 590*chi^10 + 133*chi^15)*m^5 + 4*(96 - 384*chi + 312*chi^2 + 278*chi^3 - 276*chi^4 + 18*chi^5 - 194*chi^6 + 84*chi^7 - 9*chi^9 + 126*chi^10 - 15*chi^11 - 43*chi^15 + 7*chi^21)*m^6 + (-96 + 384*chi - 384*chi^2 - 160*chi^3 + 314*chi^4 - 48*chi^5 + 112*chi^6 - 120*chi^7 + 12*chi^8 + 24*chi^9 - 80*chi^10 + 24*chi^11 - 3*chi^12 + 32*chi^15 - 4*chi^16 - 8*chi^21 + chi^28)*m^7)) / (-1 + n)^4

  kappa5 <- rowSums(1 / (m^4 * (1 + (-1 + chi)*m)^5) * (1 + (-25 + 511*chi)*m + 30*(8 - 119*chi + 311*chi^3)*m^2 + 5*(-248 + 2540*chi - 561*chi^2 - 7208*chi^3 + 6821*chi^6)*m^3 + (3904 - 29880*chi + 15690*chi^2 + 68000*chi^3 - 12990*chi^4 - 86865*chi^6 + 42525*chi^10)*m^4 + (-7872 + 49360*chi - 39660*chi^2 - 81110*chi^3 + 46460*chi^4 + 98270*chi^6 - 13365*chi^7 - 74910*chi^10 + 22827*chi^15)*m^5 +  20*(512 - 2832*chi + 2898*chi^2 + 3168*chi^3 - 3654*chi^4 + 216*chi^5 - 3109*chi^6 + 1499*chi^7 - 240*chi^9 + 2953*chi^10 - 315*chi^11 - 1390*chi^15 + 294*chi^21)*m^6 + 10*(-832 + 4288*chi - 5136*chi^2 - 2940*chi^3 + 6222*chi^4 - 1008*chi^5 + 2276*chi^6 - 2778*chi^7 + 172*chi^8 + 806*chi^9 - 2636*chi^10 + 842*chi^11 - 65*chi^12 - 90*chi^13 + 1420*chi^15 - 140*chi^16 - 476*chi^21 + 75*chi^28)*m^7 + 5*(768 - 3840*chi + 5184*chi^2 + 1152*chi^3 - 5656*chi^4 + 1728*chi^5 - 936*chi^6 + 2432*chi^7 - 420*chi^8 - 960*chi^9 + 1448*chi^10 - 912*chi^11 + 186*chi^12 + 192*chi^13 - 720*chi^15 + 170*chi^16 - 12*chi^18 + 288*chi^21 - 28*chi^22 - 73*chi^28 + 9*chi^36)*m^8 + (-768 + 3840*chi - 5760*chi^2 + 320*chi^3 + 5280*chi^4 - 2536*chi^5 + 560*chi^6 - 2240*chi^7 + 840*chi^8 + 980*chi^9 - 1072*chi^10 + 880*chi^11 - 300*chi^12 - 210*chi^13 + 400*chi^15 -180*chi^16 + 20*chi^17 + 40*chi^18 - 170*chi^21 + 40*chi^22 + 50*chi^28 - 5*chi^29 - 10*chi^36 + chi^45)*m^9)) / (-1 + n)^5

  kappa.mat <- cbind(kappa1, kappa2, kappa3, kappa4, kappa5)

  rownames(kappa.mat) <- rownames(m)

  colnames(kappa.mat) <- c("K1", "K2", "K3", "K4", "K5")

  return(kappa.mat)
}

CF.Coefficients.Fano <- function(k2, k3, k4, k5, fano, gene.names){
  c1 <- 1-fano-k3/(6*k2)+17*k3^3/(324*k2^4)-k3*k4/(12*k2^3)+k5/(40*k2^2)
  c2 <- sqrt(k2)+5*k3^2/(36*k2^(5/2))-k4/(8*k2^(3/2))
  c3 <- k3/(6*k2)-53*k3^3/(324*k2^4)+5*k3*k4/(24*k2^3)-k5/(20*k2^2)
  c4 <- -k3^2/(18*k2^(5/2))+k4/(24*k2^(3/2))
  c5 <- k3^3/(27*k2^4)-k3*k4/(24*k2^3)+k5/(120*k2^2)

  cfc <- cbind(c1,c2,c3,c4,c5)
  rownames(cfc) <- gene.names
  colnames(cfc) <- c("CF1","CF2","CF3","CF4","CF5")

  return(cfc)
}


CF.FindRoot <- function(cfc){

  roots <- polyroot(c(cfc[1], cfc[2], cfc[3], cfc[4], cfc[5]))

  re.roots <- ifelse( abs(Im(roots)) < 0.00001, Re(roots),NA)

  root <- ifelse( !all(is.na(re.roots)), min(abs(re.roots), na.rm=T), NA)

  return(root)
}

CF.AllRoots <- function(cfc){
  num.genes <- dim(cfc)[1]
  roots <- rep(0, num.genes)

  for(i in 1:num.genes){
    cfc.row <- cfc[i,]
    roots[i] <- CF.FindRoot(cfc.row)
  }

  return(roots)
}

CF.pval <- function(root){

  p<- ifelse(root>8.2,  0.5 * exp(-(root^2) / 2), 1 - pnorm(root))

  return(p)
}

Fano.BH <- function(p.df, num.genes){

  ordered <- p.df[order(p.df$pval),]

  corrected.p <- ordered$pval*num.genes/c(1:length(p.df$pval))

  ordered$BH.Corrected.Pvalue <- corrected.p

  return(ordered)
}

Fano.HighlyVariable <- function(ordered.fanos, alpha, min.fano){
  highly.variable <- apply(ordered.fanos, 1,
                           function(x) ifelse((x[4]<=alpha)&(x[1]>=min.fano),
                                              T, F))

  feature.mat <- ordered.fanos
  feature.mat$highly.variable <- highly.variable

  return(feature.mat)
}
