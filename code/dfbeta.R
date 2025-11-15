my.dfbeta <- function(m){ #m is the model
  xx <- cbind(coef(m), coef(m)+t(apply(X=dfbeta(m), MARGIN=2, FUN=range)))
  colnames(xx) <- c("orig", "min", "max")
  return(xx)
}