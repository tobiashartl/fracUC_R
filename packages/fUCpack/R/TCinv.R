#' @export
#' @useDynLib fUCpack TC_inverse
TCinv <- function(Sd, B, nu, y, n, s, b) .Call(TC_inverse, S0t=Sd, S0c=B, nu=nu,
                                               y=y, n=as.integer(n), ma=s, ma_c=b)

