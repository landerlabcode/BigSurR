get.mcFanos <- function(residuals.list){

  mcFanos <- 1/(residuals.list$num.cells-1)*rowSums((residuals.list$residuals)^2)

  return(mcFanos)
}
