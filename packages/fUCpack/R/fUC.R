#' @title fUC_comp
#' @description Kalman filter based calculation of fractional trend and cycle. The
#' function evaluates the analytical solution of the Kalman filter for fractional
#' trend and cycle given some parameters and a series to be decomposed
#' @param y Series to be decomposed into trend and cycle
#' @param d Memory parameter (integration order of trend)
#' @param nu Either a scalar giving the variance-ratio of cycle shock and trend shock
#' (if corr = FALSE), or the 2x2 covariance matrix of trend and cycle innovations Q.
#' @param corr Either TRUE or FALSE. If TRUE, trend and cycle are allowed to be correlated,
#' and nu must specify the 2x2 covariance matrix. If FALSE, trend and cycle are uncorrelated,
#' and nu can be a scalar.
#' @param ar Autoregressive coefficients of the cycle
#'
#'
#' @return A list with elements x, the estimated trend component, c, the estimated cycle component,
#' and v, the one-step ahead prediction error
#'
#' @export
fUC_comp <- function(y, d, nu, ar, corr = FALSE){


    n <- NROW(y)
    if(is.null(ar)) ar <- 0
    if(length(ar) < n-1){
        b <- c(ar, rep(0, n-length(ar)-1))
    }else{
        b <- ar
    }

    # Generate coefficient matrices
    S0t      <- embed0(ma_inf(d[1], n-1) , n)
    S0c      <- embed0(c(1, b), n)
    # Obtain coefficients for x, c
    ma       <- ma_inf(d[1],NROW(y))[-1]
    ma_c     <- c(ar, rep(0, n-length(ar)))

    # calculate Kalman-Filter
    ### Note: We use an algebraic identity that is computationally superior
    if(!corr){
        TC <- TCinv(S0t, S0c, nu, y, as.integer(n), ma, ma_c)[-(n+1),]
    }else{
        if(!is.matrix(nu)){
            warning('nu should be a 2x2 covariance matrix. It is now treated as a vector consisting of (sigma_eta, sigma_eps, sigma_eta_eps)')
            if(length(nu)!=3) error('nu is neither a covariance matrix, nor a 3-vector')
            TC <- TCinv_cor(S0t, S0c, nu[1], nu[2], nu[3], y, as.integer(n), ma, ma_c)[-(n+1),]
        }else{
            TC <- TCinv_cor(S0t, S0c, nu[1,1], nu[2,2], nu[1,2], y, as.integer(n), ma, ma_c)[-(n+1),]
        }
    }
    # return trend, cycle, prediction error
    return(list(x = TC[,1],
                c = TC[,2],
                v = y - TC[,1] - TC[,2]))
}



#' @title fUC_opt
#' @description Objective function for CSS estimator of fractional UC model
#' @param theta Parameter vector of form (d, mlogvech(Q), ar)
#' @param y Series to be decomposed into trend and cycle
#' @param nulim Interval for nu
#' @param quiet If FALSE, the value of the objective function is printed out
#' @param deterministics Either FALSE, or "const" for constant, "trend" for linear trend,
#' or "frac" for trend of order d
#' @param corr Either TRUE or FALSE. If TRUE, trend and cycle are allowed to be correlated,
#' and nu must specify the 2x2 covariance matrix. If FALSE, trend and cycle are uncorrelated,
#' and nu can be a scalar.
#' @param START First observable y that enters the objective function
#' @param return.det If true, estimates for the deterministic components are returned
#' @param d.int Interval for integration order#
#'
#'
#' @return Sum of squared prediction errors for a given parameter vector theta and an
#' observable series y
#'
#' @export
fUC_opt <- function(theta, y, nulim = c(0.05, 10), quiet = TRUE, deterministics = FALSE, corr = FALSE,
                    START=1, return.det = FALSE, d.int = c(0, 2), eta = NULL){

    d <- theta[1]
    n <- length(y)
    if(!corr){
        nu <- exp(theta[2])
        if(length(theta)>2){
            ar <- theta[-(1:2)]
        }else{
            ar <- NULL
        }
        if(nu[1] > nulim[2] | nu[1] < nulim[1])  return(.Machine$integer.max)

    }else{
        #Q <- mlogvech2mat(theta[2:4])
        #Q <- trans2mat(theta[2:4])
        # Note: Calculation of x and c is invariant to dividing Q by a constant.
        # Thus, Q[1,1] is fixed to unity
        Q  <- matrix(c(1, theta[3], theta[3], theta[2]), 2, 2)
        # check if sigma_eps is bounded
        if(theta[2] < nulim[1] | theta[2] > nulim[2]) return(.Machine$integer.max)
        # check if variance is pd
        corr <- theta[3] / sqrt(theta[2])
        if(corr > 1 | corr < -1) return(.Machine$integer.max)

        if(!is.null(eta)){
            ev <- tryCatch({
                eigen(cov2cor(Q))$values
            }, error = function(e) return(NA))
            if(any(is.na(ev)) | any(ev < 2*eta)) return(.Machine$integer.max)
        }

        nu <- c(Q[1,1], Q[2,2], Q[1,2])
        if(length(theta)>3){
            ar <- theta[-(1:3)]
        }else{
            ar <- NULL
        }
    }


    #checks
    if(d <= d.int[1] | d >= d.int[2]) return(.Machine$integer.max)
    if(!is.null(ar)){
        if(!toComp(-ar)$stable) return(.Machine$integer.max)
        if(!is.null(eta)){
            if(any(abs(toComp(-ar)$eigv)>(1-eta))) return(.Machine$integer.max)
        }
    }

    if(!corr){
        v <- fUC_comp(y, d, nu, ar, corr)$v
    }else{
        v <- fUC_comp(y, d, Q, ar, corr)$v
    }

    # account for deterministic terms
    # default is "no deterministics". (FALSE)
    # prio2 is "deterministic trend of order d" (TRUE)
    # prio3 is "constant" (constant)
    # prio4 is "constant + linear trend" (trend)

    if(is.logical(deterministics)){
        if(deterministics){
            det.type <- "frac"
        }else{
            det.type <- "none"
        }
    }else{
        det.type <- deterministics
        deterministics = TRUE
    }
    if(deterministics){
        if(det.type == "frac"){
            Z <- cbind(frac_diff(rep(1, length(y)), -d))
            if(!corr){
                v2 <- apply(Z, 2, function(x) fUC_comp(x, d, nu, ar)$v)
            }else{
                v2 <- apply(Z, 2, function(x) fUC_comp(x, d, Q, ar, corr)$v)
            }
            regmod <- lm(v[-(1:(START-1))] ~ -1+I(v2[-(1:(START-1)), ]))
            mu <- summary(regmod)$coefficients[,1]
            v  <- c(rep(NA, START-1), v[-(1:(START-1))] - v2[-(1:(START-1)), ] %*% matrix(mu, nrow = 1))

        }else{
            if(det.type == "const"){
                Z <- matrix(1, n, 1)
            }
            if(det.type == "trend"){
                Z <- cbind(1:n)
            }
            if(!corr){
                v2 <- apply(Z, 2, function(x) fUC_comp(x, d, nu, ar)$v)
            }else{
                v2 <- apply(Z, 2, function(x) fUC_comp(x, d, Q, ar, corr)$v)
            }
            regmod <- lm(v[-(1:(START-1))] ~ -1+I(v2[-(1:(START-1)), ]))
            mu <- summary(regmod)$coefficients[,1]
            v  <- c(rep(NA, START-1), v[-(1:(START-1))] - v2[-(1:(START-1)), ] %*% matrix(mu, nrow = 1))
        }
    }
    if(return.det) return(mu)


    if(!quiet) cat("ll = ", sum(v[START:n]^2), " , theta = ", theta, "\n")
    return(sum(v[START:n]^2))
}


#' @title fUC_opt_ML
#' @description Objective function for QML estimator of fractional UC model
#' @param theta Parameter vector of form (d, mlogvech(Q), ar)
#' @param y Series to be decomposed into trend and cycle
#' @param nulim Interval for nu
#' @param quiet If FALSE, the value of the objective function is printed out
#' @param approx If TRUE, the steady state Kalman filter will be used once the variance of
#' the prediction error has sufficiently converged
#' @param deterministics Either FALSE, or "const" for constant, "trend" for linear trend,
#' or "frac" for trend of order d
#' @param corr Either TRUE or FALSE. If TRUE, trend and cycle are allowed to be correlated,
#' and nu must specify the 2x2 covariance matrix. If FALSE, trend and cycle are uncorrelated,
#' and nu can be a scalar.
#' @param START First observable y that enters the objective function
#' @param return.det If true, estimates for the deterministic components are returned
#' @param d.int Interval for integration order#
#' @param diffuse If TRUE, the Kalman filter will be initialized diffusely
#' @param nu.opt If TRUE, optimization will be conducted setting the variance of
#' the long run innovations to unity. Highly recommended to set this to FALSE, unless
#' the QML estimator is compared to the CSS estimator in a simulation
#' @param Q.trans Type of transformation for covariance matrix in the estimation. Use
#' of "mlv" is highly recommended
#'
#'
#' @return Negative log likelihood for a given parameter vector theta and an
#' observable series y
#'
#' @export
fUC_opt_ML <- function(theta, y, nulim = c(0.05, 10), quiet = TRUE, approx = TRUE,
                       corr = FALSE, deterministics = FALSE, START=1, diffuse = FALSE,
                       return.det = FALSE, nu.opt = FALSE, d.int = c(0, 2),
                       eta = NULL, Q.trans = "mlv", flip=F, neg = F){
    d <- theta[1]
    n <- length(y)


    # =========================================================================
    # Parameter setup
    if(!corr){
        if(nu.opt){
            nu  <- exp(theta[2])
            Q   <- diag(c(1, nu))
            if(nu[1] > nulim[2] | nu[1] < nulim[1])  return(.Machine$integer.max)
        }else{
            sigma_par <- exp(theta[2:3])
            Q <- diag(sigma_par)
            theta <- theta[-2]
            if(sigma_par[1] > nulim[2] | sigma_par[1] < nulim[1])  return(.Machine$integer.max)
            if(sigma_par[2] > nulim[2] | sigma_par[2] < nulim[1])  return(.Machine$integer.max)
            nu <- sigma_par[2]/sigma_par[1]
        }

        Ht <- matrix(0, 1, 1)

        if(length(theta)>2){
            ar <- theta[-(1:2)]
        }else{
            ar <- NULL
        }
    }else{
        if(Q.trans == "beta"){
            if(theta[2]<=0 | theta[4]<0) return(.Machine$integer.max)
        }
        Q  <- pd_cov(theta[2:4], Q.trans,flip,neg)
        if(any(is.nan(Q))) return(.Machine$integer.max)
        if(!is.null(eta)){
            ev <- tryCatch({
                eigen(cov2cor(Q))$values
            }, error = function(e) return(NA))
            if(any(is.na(ev)) | any(ev < 2*eta)) return(.Machine$integer.max)
        }
        #Q <- trans2mat(theta[2:4])
        Ht <- matrix(0, 1, 1)
        nu <- c(Q[1,1], Q[2,2], Q[1,2])
        if(length(theta)>4){
            ar <- theta[-(1:4)]
        }else{
            ar <- NULL
        }
        if(nu[1] > nulim[2] | nu[1] < nulim[1] | nu[2]>nulim[2] | nu[2] < nulim[1])  return(.Machine$integer.max)

    }


    #checks
    if(d <= d.int[1] | d >= d.int[2]) return(.Machine$integer.max)
    if(!is.null(ar)){
        if(!toComp(-ar)$stable) return(.Machine$integer.max)
        if(!is.null(eta)){
            if(any(abs(toComp(-ar)$eigv)>(1-eta))) return(.Machine$integer.max)
        }
    }


    # =========================================================================
    # Kalman filter: prediction error
    if(!corr){
        v <- fUC_comp(y, d, nu, ar, corr=FALSE)$v
    }else{
        v <- fUC_comp(y, d, Q, ar, corr=TRUE)$v
    }



    # =========================================================================
    # Steady state KF: prediction error variance

    # iterations for P of the Kalman filter:
    Tt <- rbind(-frac_diff(c(1, rep(0, length(y)-1)), d)[-1],
                cbind(diag(length(y)-2), 0))
    if(!is.null(ar)){
        if(length(ar)>1){
            Tt_c <- rbind(-ar, cbind(diag(length(ar)-1), 0))
        }else{
            Tt_c <- matrix(-ar, 1, 1)
        }
        Tt <- bdiag(Tt, Tt_c)
    }else{
        Tt <- bdiag(Tt, 1)
    }

    Zt <- matrix(c(1, rep(0, length(y)-2)), ncol = length(y)-1)
    if(!is.null(ar)){
        Zt <- cbind(Zt, matrix(c(1, rep(0, length(ar)-1)), nrow = 1))
    }else{
        Zt <- cbind(Zt, 1)
    }

    Rt <- matrix(c(1, rep(0, NROW(Tt)-1)), ncol = 1)
    if(!is.null(ar)){
        Rt <- cbind(Rt, c(rep(0, length(y)-1), 1, rep(0, length(ar)-1)))
    }else{
        Rt <- cbind(Rt, c(rep(0, length(y)), 1))
    }


    P_v <- matrix(NA, length(y), 1)
    # initialization:
    if(diffuse){
        A <- toComp(-theta[-(1:4)])$CompMat
        p <- length(theta)-4
        sigma_sq <- Q[2,2]
        P1A <- matrix(solve(diag(p^2)-(A%x%A))%*%c(sigma_sq,rep(0,p^2-1)),p,p)
        P1Inf <- bdiag(diag(n-1)*0, P1A)
        Pt    <- P1Inf
    }else{
        Pt <- matrix(0, NCOL(Tt), NCOL(Tt))
        #     Pt[1,1] <- Q[1,1]
        #     Pt[length(y), length(y)] <- Q[2,2]

    }


    for(t in 1:n){
        P_tl <- tcrossprod(tcrossprod(Tt, Pt), Tt) + tcrossprod(tcrossprod(Rt, Q), Rt)
        F_t  <- tcrossprod(tcrossprod(Zt, P_tl), Zt)
        Pt   <- P_tl - tcrossprod(P_tl, Zt) %*% tcrossprod(Zt, P_tl) / c(F_t)
        P_v[t,] <- F_t
        if(approx & t>1){
            if(abs(P_v[t,1]/P_v[t-1,1]-1 )<0.001){
                P_v[((t+1):n), 1] <- P_v[t, 1]
                break
            }
        }
    }


    # =========================================================================
    # account for deterministic terms
    # default is "no deterministics". (FALSE)
    # prio2 is "deterministic trend of order d" (TRUE)
    # prio3 is "constant" (constant)
    # prio4 is "constant + linear trend" (trend)

    if(is.logical(deterministics)){
        if(deterministics){
            det.type <- "frac"
        }else{
            det.type <- "none"
        }
    }else{
        det.type <- deterministics
        deterministics = TRUE
    }
    if(deterministics){
        if(det.type == "frac"){
            Z <- cbind(frac_diff(rep(1, length(y)), -d))
            if(!corr){
                v2 <- apply(Z, 2, function(x) fUC_comp(x, d, nu, ar)$v)
            }else{
                v2 <- apply(Z, 2, function(x) fUC_comp(x, d, Q, ar, corr)$v)
            }
            regmod <- lm(v[-(1:(START-1))]/sqrt(P_v[-(1:(START-1))]) ~ -1+
                             I(v2[-(1:(START-1)), ]/sqrt(P_v[-(1:(START-1))])))
            mu <- summary(regmod)$coefficients[,1]
            v  <- c(rep(NA, START-1), v[-(1:(START-1))] - v2[-(1:(START-1)), ] %*% matrix(mu, nrow = 1))

        }else{
            if(det.type == "const"){
                Z <- matrix(1, n, 1)
            }
            if(det.type == "trend"){
                Z <- cbind(as.numeric(1:n))
            }
            if(!corr){
                v2 <- apply(Z, 2, function(x) fUC_comp(x, d, nu, ar)$v)
            }else{
                v2 <- apply(Z, 2, function(x) fUC_comp(x, d, Q, ar, corr)$v)
            }
            regmod <- lm(v[-(1:(START-1))]/sqrt(P_v[-(1:(START-1))]) ~ -1+
                             I(v2[-(1:(START-1)), ]/sqrt(P_v[-(1:(START-1))])))
            mu <- summary(regmod)$coefficients[,1]
            v  <- c(rep(NA, START-1), v[-(1:(START-1))] - v2[-(1:(START-1)), ] %*% matrix(mu, nrow = 1))

        }
    }
    if(return.det) return(mu)









    # checks: KFAS: ALL EQUAL!!!!
    # library(KFAS)
    # Pt <- matrix(0, NCOL(Tt), NCOL(Tt))
    # Pt[1,1] <-1
    # Pt[n+1, n+1] <- nu
    # mod <- SSModel(y ~ -1+SSMcustom(Z = Zt, T=Tt, R=Rt, Q=Qt, a1 = rep(0, NCOL(Tt)), P1=Pt), H = Ht)
    # smo <- KFS(mod, filtering = c("state"))
    # cbind(smo$v, v)
    # cbind(c(smo$F), P_v)
    # cbind(smo$P[1,1,-51], P_v)

    # calculate likelihood
    ll <- 1/2 * sum(log(P_v[START:n,])) + 1/2 * sum(v[START:n]^2/P_v[START:n,])

    if(!quiet){
        corr <- cov2cor(Q)[2,1]
        cat("log L = ", ll, ", theta = ", theta, ", corr = ", corr, "\n")
    }

    return(ll)
}








#' @title fUC_smooth
#' @description Kalman smoother based calculation of fractional trend and cycle. The
#' function evaluates the analytical solution of the Kalman smoother for fractional
#' trend and cycle given some parameters and a series to be decomposed
#' @param y Series to be decomposed into trend and cycle
#' @param d Memory parameter (integration order of trend)
#' @param nu Either a scalar giving the variance-ratio of cycle shock and trend shock
#' (if corr = FALSE), or the 2x2 covariance matrix of trend and cycle innovations Q.
#' @param corr Either TRUE or FALSE. If TRUE, trend and cycle are allowed to be correlated,
#' and nu must specify the 2x2 covariance matrix. If FALSE, trend and cycle are uncorrelated,
#' and nu can be a scalar.
#' @param ar Autoregressive coefficients of the cycle
#'
#'
#' @return A list with elements x, the estimated trend component, c, the estimated cycle component,
#' and v, the one-step ahead prediction error
#'
#' @export
fUC_smooth <- function(y, d, nu, ar, corr){


    n <- NROW(y)
    if(is.null(ar)) ar <- 0
    if(length(ar) < n-1){
        b <- c(ar, rep(0, n-length(ar)-1))
    }else{
        b <- ar
    }

    # Generate coefficient matrices
    S0t      <- embed0(ma_inf(d[1], n-1) , n)
    S0c      <- embed0(c(1, b), n)


    # calculate Kalman-Smoother
    if(!corr){
        inv <- solve(S0t%*%t(S0t)*nu + S0c %*% t(S0c))
        tau_t <- inv%*%S0c%*%t(S0c)%*%y[n:1]
        tau_c <- inv%*%S0t%*%t(S0t)%*%y[n:1]*nu
    }else{
        if(length(nu) == 3){
            a <- nu[1]
            b <- nu[2]
            c <- nu[3]
        }else{
            a <- nu[1,1]
            b <- nu[2,2]
            c <- nu[1,2]
        }
        inv <- solve(S0t%*%t(S0t)*b + S0c %*% t(S0c)*a + (S0c %*% t(S0t) + S0t %*% t(S0c))*c)
        tau_t <- inv %*% (a * S0c %*% t(S0c) + c * S0t %*% t(S0c)) %*% y[n:1]
        tau_c <- inv %*% (b * S0t %*% t(S0t) + c * S0c %*% t(S0t)) %*% y[n:1]
    }




    TC  <- cbind(tau_t[n:1], tau_c[n:1])

    # return trend, cycle, prediction error
    return(list(x = TC[,1],
                c = TC[,2],
                v = y - TC[,1] - TC[,2]))
}



UC_opt <- function(theta, y, nulim = c(0.05, 10), quiet = TRUE){
    d <- 1
    nu <- exp(theta[1])
    if(length(theta)>1){
        ar <- theta[-(1:1)]
    }else{
        ar <- NULL
    }

    #checks
    if(d <= 0 | d >= 10) return(.Machine$integer.max)
    if(!is.null(ar)){
        if(!toComp(-ar)$stable) return(.Machine$integer.max)
    }
    if(nu > nulim[2] | nu < nulim[1])  return(.Machine$integer.max)

    v <- fUC_comp(y, d, nu, ar)$v
    if(!quiet) cat("ll = ", sum(v^2), " , theta = ", theta, "\n")
    return(sum(v^2))
}
