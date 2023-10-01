#' @export
#' @useDynLib fUCpack TC_inverse_corr
TCinv_cor <- function(Sd, B, sigma_eta, sigma_eps, cov_etaeps,
                  y, n, s, b) .Call(TC_inverse_corr, S0t=Sd, S0c=B, sigma_eta = sigma_eta,
                                    sigma_eps = sigma_eps, cov_etaeps = cov_etaeps,
                                    y=y, n=as.integer(n), ma=s, ma_c=b)

