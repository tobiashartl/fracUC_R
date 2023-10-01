fUC_opt_ML_ARMA_approx <- function(theta, y, nulim = c(0.05, 10), quiet = TRUE, approx = TRUE,
                                   corr = FALSE, deterministics = FALSE, START=1, diffuse = FALSE,
                                   return.det = FALSE, nu.opt = TRUE, d.int = c(0, 2),
                                   eta = NULL, Q.trans = "mlv", flip=F, neg = F, pq = c(2,0), m=30, penalty.corr=TRUE,
                                   seas = F, period = NULL, irregular = F){
    d <- (theta[1])
    n <- length(y)
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
            if(pq[1]>0){
                ar <- theta[3:(2+pq[1])]
            }else{
                ar <- NULL
            }
            if(pq[2]>0){
                ma <- theta[(3+pq[1]):(2+pq[1]+pq[2])]
            }else{
                ma <- NULL
            }
        }else{
            ar <- NULL
            ma <- NULL
        }
        penalty <- 0
        if(seas) sigma_gamma <- exp(theta[(3+pq[1]+pq[2])])
        if(irregular) sigma_u <- exp(theta[(3+pq[1]+pq[2] + seas)])
    }else{
        Q  <- pd_cov(theta[2:4], Q.trans,flip,neg)
        if(any(is.nan(Q))) return(.Machine$integer.max)
        #if(!is.null(eta)){
        #    if(any(eigen(cov2cor(Q))$values< 2*eta)) return(.Machine$integer.max)
        #}
        # penalize correlation if too close to +-1
        if(penalty.corr){
            corr <- cov2cor(Q)[2,1]
            if(abs(corr) > .99){
                suppressWarnings( penalty <- - log( (1 - abs(corr))/.01 ) * 100)
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
        if(length(theta)>4){
            if(pq[1]>0){
                ar <- theta[5:(4+pq[1])]
            }else{
                ar <- NULL
            }
            if(pq[2]>0){
                ma <- theta[(5+pq[1]):(4+pq[1]+pq[2])]
            }else{
                ma <- NULL
            }
        }else{
            ar <- NULL
            ma <- NULL
        }
        if(nu[1] > nulim[2] | nu[1] < nulim[1] | nu[2]>nulim[2] | nu[2] < nulim[1])  return(.Machine$integer.max)
        if(irregular) sigma_u <- exp(theta[(5+pq[1]+pq[2])])
    }
    
    #checks
    if(d <= d.int[1] | d >= d.int[2]) return(.Machine$integer.max)
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
    data(list=paste("d2arma",3,3,"_n",length(y),sep =""))
    arma_d <- get(paste("d2arma",3,3,sep=""))(d)
    
    # fractional trend
    Tt <- t(toComp(c(arma_d$ar,0))$CompMat)
    Zt <- matrix(c(1,0,0,0), nrow=1)
    P1 <- diag(4)*0
    Rt <- matrix(c(1, arma_d$ma))
    a1 <- rep(0, 4)
    
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
    
    if(irregular){
        Ht <- matrix(sigma_u, 1, 1)
    }else{
        Ht <- matrix(0, 1, 1)
    }
    
    
    if (diffuse) {
        A <- toComp(-theta[5:(4+pq[1])])$CompMat
        
        sigma_sq <- Q[2, 2]
        P1A <- matrix(solve(diag(pq[1]^2) - (A %x% A)) %*% c(sigma_sq, 
                                                             rep(0, pq[1]^2 - 1)), pq[1], pq[1])
        P1Inf <- bdiag(tcrossprod(rep(1, 4))*0, P1A)
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
    if(seas){
        Dummies <- embed0(rep(c(1, rep(0, period-1)), length.out = length(y)), period-1)
    }else{
        Dummies <- NULL
    }
    if(deterministics){
        if(det.type == "frac"){
            Z <-  cbind(frac_diff(rep(1, length(y)), -d), Dummies)
            
            if(!corr){
                v2 <- apply(matrix(Z, nrow = n), 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
            }else{
                v2 <- apply(matrix(Z, nrow = n), 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
            }
            
            regmod <- lm(v[-(1:(START-1))]/sqrt(Pt[-(1:(START-1))]) ~ -1+I(v2[-(1:(START-1)), ]/sqrt(Pt[-(1:(START-1))])))
            mu <- summary(regmod)$coefficients[,1]
            v  <- c(rep(NA, START-1), v[-(1:(START-1))] - v2[-(1:(START-1)), ] %*% matrix(mu, nrow = 1))
        }else{
            if(det.type == "const"){
                Z <- cbind(as.double(matrix(1, n, 1)), Dummies)
                if(!corr){
                    v2 <- apply(matrix(Z, nrow = n), 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
                }else{
                    v2 <- apply(matrix(Z, nrow = n), 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
                }
                
            }
            if(det.type == "trend"){
                Z <- cbind(as.double(matrix(1:n, n, 1)), Dummies)
                if(!corr){
                    v2 <- apply(matrix(Z, nrow = n), 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
                    
                }else{
                    v2 <- apply(matrix(Z, nrow = n), 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
                    
                }
                
            }
            
            if(det.type == "quadratic"){
                Z <- cbind(as.double(1:n), as.double((1:n)^2), Dummies)
                if(!corr){
                    v2 <- apply(Z, 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
                }else{
                    v2 <-apply(Z, 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
                }
                
            }
            
            if(det.type == "linfrac"){
                Z <- cbind(1:n, frac_diff(rep(1, n), -d), Dummies)
                if(!corr){
                    v2 <- apply(Z, 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
                }else{
                    v2 <-apply(Z, 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
                }
                
            }
            
            if (det.type == "trend_corona") {
                Z <- cbind(1:n, 0, Dummies)
                Z[294, 2] <- 1
                if(!corr){
                    v2 <- apply(Z, 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
                }else{
                    v2 <- apply(Z, 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
                }
                
            }
            
            if (det.type == "trend_corona_break") {
                Z <- cbind((1:n), 0, c(rep(0, 215), 1:(n-215)), Dummies)
                Z[294, 2] <- 1
                if(!corr){
                    v2 <- apply(Z, 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
                }else{
                    v2 <- apply(Z, 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
                }
                
            }
            
            
            
            regmod <- lm(v[-(1:(START-1))]/sqrt(Pt[-(1:(START-1))]) ~ -1+I(v2[-(1:(START-1)), ]/sqrt(Pt[-(1:(START-1))])))
            mu <- summary(regmod)$coefficients[,1]
            v  <- c(rep(NA, START-1), v[-(1:(START-1))] - v2[-(1:(START-1)), ] %*% matrix(mu, nrow = 1))
        }
    }else{
        if(seas){
            #Dummies <- embed0(rep(c(1, rep(0, period-1)), length.out = length(y)), period)
            Z <- Dummies
            if(!corr){
                v2 <- apply(matrix(Z, nrow = n), 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
            }else{
                v2 <- apply(matrix(Z, nrow = n), 2, function(x) fUC_KF_approx(x, Tt, Zt, Rt, a1, P1, Q, Ht)$v)
            }
            regmod <- lm(v[-(1:(START-1))]/sqrt(Pt[-(1:(START-1))]) ~ -1+I(v2[-(1:(START-1)), ]/sqrt(Pt[-(1:(START-1))])))
            mu <- summary(regmod)$coefficients[,1]
            v  <- c(rep(NA, START-1), v[-(1:(START-1))] - v2[-(1:(START-1)), ] %*% matrix(mu, nrow = 1))
        }else{
            mu <- NULL
        }
    }
    if(return.det) return(mu)
    
    Pt <- Pt + c(Ht)
    ll <- 1/2 * sum(log(Pt[START:n])) + 1/2 * sum(v[START:n]^2/Pt[START:n]) + penalty
    return(ll)
}







fUC_KF_approx <- function(y, Tt, Zt, Rt, at, Pt, Q, Ht=matrix(0, 1, 1)){
    # initialization:
    P_v <- matrix(NA, length(y), 1)
    TC <- matrix(NA, length(y), length(at))
    V  <- matrix(NA, length(y), 1)
    n <- length(y)
    #Pt[1,1] <-1
    #Ft <- Pt[1,1] + Ht
    #P_v[1,1] <- Pt[1,1]
    
    for(t in 1:n){
        P_tl <- tcrossprod(tcrossprod(Tt, Pt), Tt) + tcrossprod(tcrossprod(Rt, Q), Rt)
        a_tl <- Tt %*% at #+ tcrossprod(tcrossprod(Tt, Pt), Zt) / c(F_t)
        v_t  <- y[t] - Zt %*% a_tl
        
        
        F_t  <- tcrossprod(tcrossprod(Zt, P_tl), Zt) + Ht
        at   <- a_tl + tcrossprod(P_tl, Zt) / c(F_t) * c(v_t)
        Pt   <- P_tl - tcrossprod(P_tl, Zt) %*% tcrossprod(Zt, P_tl) / c(F_t)
        
        P_v[t,] <- c(F_t)
        V[t,] <- c(v_t)
        TC[t,] <- c(a_tl)
    }
    return(list(v = c(V),
                TC = TC,
                P_t = c(P_v)))
}


fUC_KS_ARMA_approx <- function(theta, y, nulim = c(0.05, 10), quiet = TRUE, approx = TRUE,
                               corr = FALSE, deterministics = FALSE, START=1, diffuse = FALSE,
                               return.det = FALSE, nu.opt = TRUE, d.int = c(0, 2),
                               eta = NULL, Q.trans = "mlv", flip=F, neg = F, pq = c(2,0), m=30, penalty.corr=TRUE,
                               seas = F, period = 4, irregular = F){
    d <- (theta[1])
    n <- length(y)
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
            if(pq[1]>0){
                ar <- theta[3:(2+pq[1])]
            }else{
                ar <- NULL
            }
            if(pq[2]>0){
                ma <- theta[(3+pq[1]):(2+pq[1]+pq[2])]
            }else{
                ma <- NULL
            }
        }else{
            ar <- NULL
            ma <- NULL
        }
        penalty <- 0
    }else{
        Q  <- pd_cov(theta[2:4], Q.trans,flip,neg)
        if(any(is.nan(Q))) return(.Machine$integer.max)
        #if(!is.null(eta)){
        #    if(any(eigen(cov2cor(Q))$values< 2*eta)) return(.Machine$integer.max)
        #}
        # penalize correlation if too close to +-1
        if(penalty.corr){
            corr <- cov2cor(Q)[2,1]
            if(abs(corr) > .99){
                suppressWarnings( penalty <- - log( (1 - abs(corr))/.01 ) * 100)
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
        if(length(theta)>4){
            if(pq[1]>0){
                ar <- theta[5:(4+pq[1])]
            }else{
                ar <- NULL
            }
            if(pq[2]>0){
                ma <- theta[(5+pq[1]):(4+pq[1]+pq[2])]
            }else{
                ma <- NULL
            }
        }else{
            ar <- NULL
            ma <- NULL
        }
        if(nu[1] > nulim[2] | nu[1] < nulim[1] | nu[2]>nulim[2] | nu[2] < nulim[1])  return(.Machine$integer.max)
        
    }
    
    #checks
    if(d <= d.int[1] | d >= d.int[2]) return(.Machine$integer.max)
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
    data(list=paste("d2arma",3,3,"_n",length(y),sep =""))
    arma_d <- get(paste("d2arma",3,3,sep=""))(d)
    
    # fractional trend
    Tt <- t(toComp(c(arma_d$ar,0))$CompMat)
    Zt <- matrix(c(1,0,0,0), nrow=1)
    P1 <- diag(4)*0
    Rt <- matrix(c(1, arma_d$ma))
    a1 <- rep(0, 4)
    
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
        sigma_gamma <- exp(theta[(3+pq[1]+pq[2])])
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
        Ht <- matrix(exp(theta[3 + pq[1] + pq[2] + seas]), 1, 1)
    }else{
        Ht <- matrix(0, 1, 1)
    }
    
    
    
    Mod    <- SSModel(y ~ -1 + SSMcustom(Zt, Tt, Rt, Q, a1, P1), H=Ht)
    TCcomp <- KFS(Mod, smoothing = "state")
    TC      <- TCcomp$alphahat
    x      <- TC[,1]
    c      <- TC[,5]
    if(seas){
        gamma <- TC[,5+pq[1]+1]
    }else{
        gamma <- NULL
    }
    return(list(x = x, 
                c = c,
                gamma=gamma))
}



fUC_KS_approx <- function(y, Tt, Zt, Rt, at, Pt, Q){
    # initialization:
    P_v <- matrix(NA, length(y), 1)
    TC <- matrix(NA, length(y), length(at))
    V  <- matrix(NA, length(y), 1)
    n <- length(y)
    #Pt[1,1] <-1
    #Ft <- Pt[1,1] + Ht
    #P_v[1,1] <- Pt[1,1]
    
    for(t in 1:n){
        P_tl <- tcrossprod(tcrossprod(Tt, Pt), Tt) + tcrossprod(tcrossprod(Rt, Q), Rt)
        a_tl <- Tt %*% at #+ tcrossprod(tcrossprod(Tt, Pt), Zt) / c(F_t)
        v_t  <- y[t] - Zt %*% a_tl
        
        
        F_t  <- tcrossprod(tcrossprod(Zt, P_tl), Zt)
        at   <- a_tl + tcrossprod(P_tl, Zt) / c(F_t) * c(v_t)
        Pt   <- P_tl - tcrossprod(P_tl, Zt) %*% tcrossprod(Zt, P_tl) / c(F_t)
        
        P_v[t,] <- c(F_t)
        V[t,] <- c(v_t)
        TC[t,] <- c(a_tl)
    }
    return(list(v = c(V),
                TC = TC,
                P_t = c(P_v)))
}

