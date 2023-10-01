UC_opt_ML <- function(theta, y, nulim = c(0.05, 10), quiet = TRUE, approx = TRUE,
                      corr = FALSE, deterministics = FALSE, START=1, diffuse = FALSE,
                      return.det = FALSE, nu.opt = TRUE, 
                      eta = NULL, Q.trans = "mlv", flip=F, neg = F, pq = c(2,0), m=30, penalty.corr=TRUE,
                      seas = F, period = NULL, irregular = F){
    n <- length(y)
    if(!corr){
        if(nu.opt){
            nu  <- exp(theta[1])
            Q   <- diag(c(1, nu))
            if(nu[1] > nulim[2] | nu[1] < nulim[1])  return(.Machine$integer.max)
        }else{
            sigma_par <- exp(theta[1:2])
            Q <- diag(sigma_par)
            theta <- theta[-1]
            if(sigma_par[1] > nulim[2] | sigma_par[1] < nulim[1])  return(.Machine$integer.max)
            if(sigma_par[2] > nulim[2] | sigma_par[2] < nulim[1])  return(.Machine$integer.max)
            nu <- sigma_par[2]/sigma_par[1]
        }
        
        Ht <- matrix(0, 1, 1)
        
        if(length(theta)>2){
            if(pq[1]>0){
                ar <- theta[2:(1+pq[1])]
            }else{
                ar <- NULL
            }
            if(pq[2]>0){
                ma <- theta[(2+pq[1]):(1+pq[1]+pq[2])]
            }else{
                ma <- NULL
            }
        }else{
            ar <- NULL
            ma <- NULL
        }
        penalty <- 0
        if(seas) sigma_gamma <- exp(theta[(2+pq[1]+pq[2])])
        if(irregular) sigma_u <- exp(theta[(2+pq[1]+pq[2] + seas)])
    }else{
        Q  <- pd_cov(theta[1:3], Q.trans,flip,neg)
        if(any(is.nan(Q))) return(.Machine$integer.max)
        #if(!is.null(eta)){
        #    if(any(eigen(cov2cor(Q))$values< 2*eta)) return(.Machine$integer.max)
        #}
        # penalize correlation if too close to +-1
        if(penalty.corr){
            correlation <- cov2cor(Q)[2,1]
            if(abs(correlation) > .99){
                suppressWarnings( penalty <- - log( (1 - abs(correlation))/.01 ) * 100)
            }else{
                penalty <- 0
            }
        }else{
            penalty <- 0
        }
        
        if(is.nan(penalty)) penalty <- .Machine$integer.max
        
        #Q <- trans2mat(theta[2:4])
        Ht <- matrix(0, 1, 1)
        nu <- c(Q[1,1], Q[2,2], Q[1,2])
        if(length(theta)>3){
            if(pq[1]>0){
                ar <- theta[4:(3+pq[1])]
            }else{
                ar <- NULL
            }
            if(pq[2]>0){
                ma <- theta[(4+pq[1]):(3+pq[1]+pq[2])]
            }else{
                ma <- NULL
            }
        }else{
            ar <- NULL
            ma <- NULL
        }
        if(nu[1] > nulim[2] | nu[1] < nulim[1] | nu[2]>nulim[2] | nu[2] < nulim[1])  return(.Machine$integer.max)
        if(seas) sigma_gamma <- exp(theta[(4+pq[1]+pq[2])])
        if(irregular) sigma_u <- exp(theta[(4+pq[1]+pq[2] + seas)])
    }
    
    #checks
    if(!is.null(ar)){
        if(!toComp(-ar)$stable) return(.Machine$integer.max)
        if (!is.null(eta)) {
            if (any(abs(toComp(-ar)$eigv) > (1 - eta)))
                return(.Machine$integer.max)
        }
    }
    
    
    if(!is.null(ma)){
        maInvert <- function(ma) {
            q <- length(ma)
            q0 <- max(which(c(1, ma) != 0)) - 1L
            if (!q0)
                return(ma)
            roots <- polyroot(c(1, ma[1L:q0]))
            ind <- Mod(roots) < 1
            if (all(!ind)){
                return(TRUE)
            }else{
                return(FALSE)
            }
        }
        if(!maInvert(ma)) return(.Machine$integer.max)
        b <- ARMAtoMA(ar = -ma, ma = ar, lag.max=n-1)
    }else{
        if(is.null(ar)){
            b <- NULL
            
        }else{
            b <- ARMAtoMA(ma = ar, lag.max=length(ar))
            
        }
    }
    
    # ====================================================================
    # Build SSM

    # fractional trend
    Tt <- t(toComp(c(1))$CompMat)
    Zt <- matrix(c(1), nrow=1)
    P1 <- diag(1)*0
    Rt <- matrix(c(1), ncol=1)
    a1 <- rep(0, 1)
    
    # cycle
    if(!is.null(ar)){
        Tt_c <- toComp(-b)$CompMat
    }else{
        Tt_c <- matrix(0, 1, 1)
    }
    Zt_c <- matrix(c(1, rep(0, ncol(Tt_c)-1)), nrow=1)
    Rt_c <- t(Zt_c)
    P1_c <- diag(ncol(Tt_c))*0
    a1_c <- rep(0, ncol(Tt_c))
    
    # Overall components
    Tt <- bdiag(Tt, Tt_c)
    Zt <- cbind(Zt, Zt_c)
    Rt <- bdiag(Rt, Rt_c)
    P1 <- bdiag(P1, P1_c)
    a1 <- c(a1, a1_c)
    
    # if seasonal is true
    if(seas){
        s <- period
        lambda <- 2 * pi * 1:s / s
        s_star <- floor((s-1)/2)
        C <- list()
        for(ss in 1:s_star){
            C[[ss]] <- matrix(c(cos(lambda[ss]), - sin(lambda[ss]), sin(lambda[ss]), 
                                cos(lambda[ss])), 2, 2
            )
        }
        Tt_gamma <- bdiag(bdiag(C), -1)
        Zt_gamma <- matrix(rep(c(1,0), length.out = ncol(Tt_gamma)), 1)
        Rt_gamma <- diag(s-1)
        Q_gamma  <- diag(s-1)*sigma_gamma
        a1_gamma <- rep(0, ncol(Tt_gamma))
        P1_gamma <- diag(ncol(Tt_gamma))*0
        
        
        Tt <- bdiag(Tt, Tt_gamma)
        Zt <- cbind(Zt, Zt_gamma)
        Rt <- bdiag(Rt, Rt_gamma)
        P1 <- bdiag(P1, P1_gamma)
        a1 <- c(a1, a1_gamma)
        Q <- bdiag(Q, Q_gamma)
    }
    
    if(irregular){
        Ht <- matrix(sigma_u, 1, 1)
    }else{
        Ht <- matrix(0, 1, 1)
    }
    
    
    if (diffuse) {
        A <- toComp(-theta[3:(2+pq[1])])$CompMat
        
        sigma_sq <- Q[2, 2]
        P1A <- matrix(solve(diag(pq[1]^2) - (A %x% A)) %*% c(sigma_sq, 
                                                             rep(0, pq[1]^2 - 1)), pq[1], pq[1])
        P1Inf <- bdiag(tcrossprod(rep(1, 1))*0, P1A)
        if(seas) P1Inf = bdiag(P1Inf, Q_gamma)
        P1 <- P1Inf
    }
    
    if(!corr){
        TCcomp <- fUC_KF_approx(y, Tt, Zt, Rt, a1, P1, Q, Ht)
        v      <- TCcomp$v
        Pt     <- TCcomp$P_t
    }else{
        TCcomp <- fUC_KF_approx(y, Tt, Zt, Rt, a1, P1, Q, Ht)
        v      <- TCcomp$v
        Pt     <- TCcomp$P_t
    }

    if(deterministics){
        Z <- as.double(matrix(1:n, n, 1))
        if(!corr){
            v2 <-matrix(fUC_KF_approx(Z, Tt, Zt, Rt, a1, P1, Q, Ht)$v, n, 1)
        }else{
            v2 <-matrix(fUC_KF_approx(Z, Tt, Zt, Rt, a1, P1, Q, Ht)$v, n, 1)
        }
                
        
        
        
        regmod <- lm(v[-(1:(START-1))]/sqrt(Pt[-(1:(START-1))]) ~ -1+I(v2[-(1:(START-1)), ]/sqrt(Pt[-(1:(START-1))])))
        mu <- summary(regmod)$coefficients[,1]
        v  <-  c(rep(NA, START-1), v[-(1:(START-1))] - v2[-(1:(START-1)), ] %*% matrix(mu, nrow = 1))
    }else{
        v2 <- rep(0, length(y))
        
        
        
        regmod <- lm(v[-(1:(START-1))] ~ -1+v2[-(1:(START-1))])
        mu <- summary(regmod)$coefficients[,1]
        v  <-  c(rep(NA, START-1), v[-(1:(START-1))] - v2[-(1:(START-1)), ] %*% matrix(mu, nrow = 1))
    }
    if(return.det) return(mu)
    
    Pt <- Pt + c(Ht)
    ll <- 1/2 * sum(log(Pt[START:n])) + 1/2 * sum(v[START:n]^2/Pt[START:n]) + penalty
    if (!quiet) {
        corr <- cov2cor(Q)[2, 1]
        cat("log L = ", ll, ", theta = ", theta, ", corr = ", 
            corr, "\n")
    }
    return(ll)
}




UC_i2_opt_ML <- function(theta, y, nulim = c(0.05, 10), quiet = TRUE, approx = TRUE,
                      corr = FALSE, deterministics = FALSE, START=1, diffuse = FALSE,
                      return.det = FALSE, nu.opt = TRUE, 
                      eta = NULL, Q.trans = "mlv", flip=F, neg = F, pq = c(2,0), m=30, penalty.corr=TRUE,
                      seas = F, period = NULL, irregular = F){
    n <- length(y)
    if(!corr){
        if(nu.opt){
            nu  <- exp(theta[1])
            Q   <- diag(c(1, nu))
            if(nu[1] > nulim[2] | nu[1] < nulim[1])  return(.Machine$integer.max)
        }else{
            sigma_par <- exp(theta[1:2])
            Q <- diag(sigma_par)
            theta <- theta[-1]
            if(sigma_par[1] > nulim[2] | sigma_par[1] < nulim[1])  return(.Machine$integer.max)
            if(sigma_par[2] > nulim[2] | sigma_par[2] < nulim[1])  return(.Machine$integer.max)
            nu <- sigma_par[2]/sigma_par[1]
        }
        
        Ht <- matrix(0, 1, 1)
        
        if(length(theta)>2){
            if(pq[1]>0){
                ar <- theta[2:(1+pq[1])]
            }else{
                ar <- NULL
            }
            if(pq[2]>0){
                ma <- theta[(2+pq[1]):(1+pq[1]+pq[2])]
            }else{
                ma <- NULL
            }
        }else{
            ar <- NULL
            ma <- NULL
        }
        penalty <- 0
        if(seas) sigma_gamma <- exp(theta[(2+pq[1]+pq[2])])
        if(irregular) sigma_u <- exp(theta[(2+pq[1]+pq[2] + seas)])
    }else{
        Q  <- pd_cov(theta[1:3], Q.trans,flip,neg)
        if(any(is.nan(Q))) return(.Machine$integer.max)
        #if(!is.null(eta)){
        #    if(any(eigen(cov2cor(Q))$values< 2*eta)) return(.Machine$integer.max)
        #}
        # penalize correlation if too close to +-1
        if(penalty.corr){
            correlation <- cov2cor(Q)[2,1]
            
            penalty <- tryCatch(if(abs(correlation) > .98){
                - log( (1 - abs(correlation))/.02 ) * 100
            }else{
                0
            }, error = function(e) .Machine$integer.max)
            if(is.nan(penalty)) penalty <- .Machine$integer.max
            if(penalty > .Machine$integer.max) penalty <- .Machine$integer.max
            
        }else{
            penalty <- 0
        }
        
        if(is.nan(penalty)) penalty <- .Machine$integer.max
        
        #Q <- trans2mat(theta[2:4])
        Ht <- matrix(0, 1, 1)
        nu <- c(Q[1,1], Q[2,2], Q[1,2])
        if(length(theta)>3){
            if(pq[1]>0){
                ar <- theta[4:(3+pq[1])]
            }else{
                ar <- NULL
            }
            if(pq[2]>0){
                ma <- theta[(4+pq[1]):(3+pq[1]+pq[2])]
            }else{
                ma <- NULL
            }
        }else{
            ar <- NULL
            ma <- NULL
        }
        if(nu[1] > nulim[2] | nu[1] < nulim[1] | nu[2]>nulim[2] | nu[2] < nulim[1])  return(.Machine$integer.max)
        if(seas) sigma_gamma <- exp(theta[(4+pq[1]+pq[2])])
        if(irregular) sigma_u <- exp(theta[(4+pq[1]+pq[2] + seas)])
    }
    
    #checks
    if(!is.null(ar)){
        if(!toComp(-ar)$stable) return(.Machine$integer.max)
        if (!is.null(eta)) {
            if (any(abs(toComp(-ar)$eigv) > (1 - eta)))
                return(.Machine$integer.max)
        }
    }
    
    
    if(!is.null(ma)){
        maInvert <- function(ma) {
            q <- length(ma)
            q0 <- max(which(c(1, ma) != 0)) - 1L
            if (!q0)
                return(ma)
            roots <- polyroot(c(1, ma[1L:q0]))
            ind <- Mod(roots) < 1
            if (all(!ind)){
                return(TRUE)
            }else{
                return(FALSE)
            }
        }
        if(!maInvert(ma)) return(.Machine$integer.max)
        b <- ARMAtoMA(ar = -ma, ma = ar, lag.max=n-1)
    }else{
        if(is.null(ar)){
            b <- NULL
            
        }else{
            b <- ARMAtoMA(ma = ar, lag.max=length(ar))
            
        }
    }
    
    # ====================================================================
    # Build SSM
    
    # fractional trend
    Tt <- t(toComp(c(2, -1))$CompMat)
    Zt <- matrix(c(1, 0), nrow=1)
    P1 <- diag(2)*0
    Rt <- matrix(c(1, 0), ncol=1)
    a1 <- rep(0, 2)
    
    # cycle
    if(!is.null(ar)){
        Tt_c <- toComp(-b)$CompMat
    }else{
        Tt_c <- matrix(0, 1, 1)
    }
    Zt_c <- matrix(c(1, rep(0, ncol(Tt_c)-1)), nrow=1)
    Rt_c <- t(Zt_c)
    P1_c <- diag(ncol(Tt_c))*0
    a1_c <- rep(0, ncol(Tt_c))
    
    # Overall components
    Tt <- bdiag(Tt, Tt_c)
    Zt <- cbind(Zt, Zt_c)
    Rt <- bdiag(Rt, Rt_c)
    P1 <- bdiag(P1, P1_c)
    a1 <- c(a1, a1_c)
    
    # if seasonal is true
    if(seas){
        s <- period
        lambda <- 2 * pi * 1:s / s
        s_star <- floor((s-1)/2)
        C <- list()
        for(ss in 1:s_star){
            C[[ss]] <- matrix(c(cos(lambda[ss]), - sin(lambda[ss]), sin(lambda[ss]), 
                                cos(lambda[ss])), 2, 2
            )
        }
        Tt_gamma <- bdiag(bdiag(C), -1)
        Zt_gamma <- matrix(rep(c(1,0), length.out = ncol(Tt_gamma)), 1)
        Rt_gamma <- diag(s-1)
        Q_gamma  <- diag(s-1)*sigma_gamma
        a1_gamma <- rep(0, ncol(Tt_gamma))
        P1_gamma <- diag(ncol(Tt_gamma))*0
        
        
        Tt <- bdiag(Tt, Tt_gamma)
        Zt <- cbind(Zt, Zt_gamma)
        Rt <- bdiag(Rt, Rt_gamma)
        P1 <- bdiag(P1, P1_gamma)
        a1 <- c(a1, a1_gamma)
        Q <- bdiag(Q, Q_gamma)
    }
    
    if(irregular){
        Ht <- matrix(sigma_u, 1, 1)
    }else{
        Ht <- matrix(0, 1, 1)
    }
    
    
    if (diffuse) {
        sigma_sq <- Q[2, 2]
        if(pq[1] > 0){ 
            A <- toComp(-theta[3:(2+pq[1])])$CompMat
            P1A <- matrix(solve(diag(pq[1]^2) - (A %x% A)) %*% c(sigma_sq, 
                                                                 rep(0, pq[1]^2 - 1)), pq[1], pq[1])
        }else{
            P1A <- matrix(sigma_sq, 1, 1)
        }
        
        
        
        P1Inf <- bdiag(diag(2)*0, P1A)
        if(seas) P1Inf = bdiag(P1Inf, Q_gamma)
        P1 <- P1Inf
    }
    
    if(!corr){
        TCcomp <- fUC_KF_approx(y, Tt, Zt, Rt, a1, P1, Q, Ht)
        v      <- TCcomp$v
        Pt     <- TCcomp$P_t
    }else{
        TCcomp <- fUC_KF_approx(y, Tt, Zt, Rt, a1, P1, Q, Ht)
        v      <- TCcomp$v
        Pt     <- TCcomp$P_t
    }
    
    if(is.logical(deterministics)){
        if(deterministics){
            Z <- as.double(matrix((1:n)^2, n, 1))
            if(!corr){
                v2 <-matrix(fUC_KF_approx(Z, Tt, Zt, Rt, a1, P1, Q, Ht)$v, n, 1)
            }else{
                v2 <-matrix(fUC_KF_approx(Z, Tt, Zt, Rt, a1, P1, Q, Ht)$v, n, 1)
            }
            
            
            
            
            regmod <- lm(v[-(1:(START-1))]/sqrt(Pt[-(1:(START-1))]) ~ -1+I(v2[-(1:(START-1)), ]/sqrt(Pt[-(1:(START-1))])))
            mu <- summary(regmod)$coefficients[,1]
            v  <-  c(rep(NA, START-1), v[-(1:(START-1))] - v2[-(1:(START-1)), ] %*% matrix(mu, nrow = 1))
        }else{
            
            v  <-  c(rep(NA, START-1), v[-(1:(START-1))] )
        }
    }else{
        if(deterministics == "double"){
            Z <- as.double(cbind(matrix(1:n, ncol=1), matrix((1:n)^2, n, 1)))
            if(!corr){
                v2 <-apply(Z, 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
            }else{
                v2 <-apply(Z, 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
            }
            regmod <- lm(v[-(1:(START-1))]/sqrt(Pt[-(1:(START-1))]) ~ -1+I(v2[-(1:(START-1)), ]/sqrt(Pt[-(1:(START-1))])))
            mu <- summary(regmod)$coefficients[,1]
            v  <-  c(rep(NA, START-1), v[-(1:(START-1))] - v2[-(1:(START-1)), ] %*% matrix(mu, nrow = 1))
        }
    }
    
    
    
    if(return.det) return(mu)
    
    Pt <- Pt + c(Ht)
    ll <- 1/2 * sum(log(Pt[START:n])) + 1/2 * sum(v[START:n]^2/Pt[START:n]) + penalty
    if (!quiet) {
        corr <- cov2cor(Q)[2, 1]
        cat("log L = ", ll, ", theta = ", theta, ", corr = ", 
            corr, "\n")
    }
    return(ll)
}
