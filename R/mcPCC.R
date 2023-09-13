get.mcPCC2 <- function(residuals.list, mcFanos){
  corrections <- as.matrix(correction.PCC(residuals.list, mcFanos))
  mcPCCs <- tcrossprod(corrections)

  return(mcPCCs)
}

correction.PCC <- function(residuals.list, mcFanos){
  out <- residuals.list$residuals/sqrt((residuals.list$num.cells-1)*mcFanos)
  return(out)
}
