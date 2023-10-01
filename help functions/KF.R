fUC_opt_KF <- function(theta, y, nulim = c(0.05, 10), quiet = TRUE){
    d <- theta[1]
    nu <- exp(theta[2])
    if(length(theta)>2){
        ar <- theta[-(1:2)]
    }else{
        ar <- NULL
    }
    
    #checks
    if(d <= 0 | d >= 10) return(.Machine$integer.max)
    if(!is.null(ar)){
        if(!toComp(ar)$stable) return(.Machine$integer.max)
    }
    if(nu > nulim[2] | nu < nulim[1])  return(.Machine$integer.max)
    
    v <- fUC_comp(y, d, nu, ar)$v
    
    # iterations for P of the Kalman filter: 
    Tt <- rbind(-frac_diff(c(1, rep(0, length(y))), d)[-1],
                cbind(diag(length(y)-1), 0))

    
    Zt <- matrix(c(1, rep(0, length(y)-1)), ncol = length(y))

    
    Qt <- matrix(1, 1, 1)
    Rt <- matrix(c(1, rep(0, NROW(Tt)-1)), ncol = 1)

    Ht <- matrix(nu, 1, 1)
    
    
    # initialization: 
    P_v <- matrix(NA, length(y), 1)
    Pt <- matrix(0, length(y), length(y))
    Pt[1,1] <-1
    Ft <- Pt[1,1] + Ht
    P_v[1,1] <- Pt[1,1]
    
    for (t in 2:length(y)){
        Pt <- Tt %*% (Pt - Pt[,1, drop = F] %*% Pt[1,,drop=F] / Ft[1,1]) %*% t(Tt) + Rt %*% t(Rt) * Qt[1,1]
        P_v[t, 1] <- Pt[1,1]
        Ft <- Pt[1,1] + Ht
    }

    ll <- 1/2 * sum(P_v) + 1/2 * sum(v^2/P_v)
    return(ll)
}







UC_opt_KF_i1 <- function(theta, y, nulim = c(0.05, 10), quiet = TRUE, 
                         return.mod = FALSE, ll = TRUE, corr = FALSE, START=1,
                         deterministics = FALSE, return.det = FALSE){
    n <- length(y)
    if(!corr){
        nu <- exp(theta[1])
        Qt <- matrix(1, 1, 1)
        Tt <- matrix(1, 1, 1)
        Zt <- matrix(1, 1, 1)
        Rt <- matrix(1, ncol = 1)
        
        if(nu[1] > nulim[2] | nu[1] < nulim[1])  return(.Machine$integer.max)
        
        if(length(theta)>1){
            ar <- theta[-(1)]
            if(length(ar)> 1){
                Tt_c <- rbind(ar, cbind(diag(length(ar)-1), 0))
            }else{
                Tt_c <- matrix(ar, 1, 1)
            }
            Tt <- bdiag(Tt, Tt_c)
            Zt <- cbind(Zt, matrix(c(1, rep(0, length(ar)-1)), nrow = 1))
            Ht <- matrix(0, 1, 1) 
            Rt <- bdiag(Rt, c(1, rep(0, length(ar)-1)))
            Qt <- diag(c(1, nu))
        }else{
            ar <- NULL
            Ht <- matrix(nu, 1, 1) 
        }
        
    }else{
        Qt <- matrix(c(1, theta[2], theta[2], theta[1]), 2, 2)
        if (theta[1] < nulim[1] | theta[1] > nulim[2]) 
            return(.Machine$integer.max)
        corr <- theta[2]/sqrt(theta[1])
        if (corr > 1 | corr < -1) 
            return(.Machine$integer.max)
        nu <- c(Qt[1, 1], Qt[2, 2], Qt[1, 2])
        Tt <- matrix(1, 1, 1)
        Zt <- matrix(1, 1, 1)
        Rt <- matrix(1, ncol = 1)
        # mod
        Ht <- matrix(0, 1, 1) 
        
        
        if(length(theta)>2){
            ar <- theta[-(1:2)]
            if(length(ar) > 1){
                Tt_c <- rbind(ar, cbind(diag(length(ar)-1), 0))
            }else{
                Tt_c <- matrix(ar, 1, 1)
            }
            Tt <- bdiag(Tt, Tt_c)
            Zt <- cbind(Zt, matrix(c(1, rep(0, length(ar)-1)), nrow = 1))
            Rt <- bdiag(Rt, c(1, rep(0, length(ar)-1)))
        }else{
            ar <- NULL
            Zt <- cbind(Zt, matrix(1, nrow =1, ncol=1))
            Tt <- bdiag(Tt, matrix(0, 1, 1))
            Rt <- bdiag(Rt, 1)
        }
    }
    
    #checks
    if(!is.null(ar)){
        if(!toComp(ar)$stable) return(.Machine$integer.max)
    }
    
    mod <- KFAS::SSModel(y ~ - 1 + SSMcustom(Z=Zt, T=Tt, R=Rt, Q=Qt, a1=rep(0, ncol(Zt)), P1=matrix(0, ncol(Tt), ncol(Tt))), H=Ht)
    
    
    
    if (is.logical(deterministics)) {
        if (deterministics) {
            det.type <- "trend"
        }else {
            det.type <- "none"
        }
    }else {
        det.type <- deterministics
        deterministics = TRUE
    }
    if (deterministics) {
        v <- KFS(mod, filtering = "state", smoothing = "none")$v
        if (det.type == "frac") {
            Z <- cbind(1, frac_diff(rep(1, length(y)), -1))
            if (!corr) {
                v2 <- apply(Z, 2, function(x) fUC_comp(x, 1, nu, ar)$v)
            }else{
                v2 <- apply(Z, 2, function(x) fUC_comp(x, 1, 
                                                       Q, ar, corr)$v)
            }
            regmod <- lm(v ~ -1 + v2)
            mu <- summary(regmod)$coefficients[, 1]
            v <- residuals(regmod)
        }else {
            if (det.type == "const") {
                Z <- matrix(1, n, 1)
            }
            if (det.type == "trend") {
                Z <- cbind(1, 1:n)
            }
            if (!corr) {
                v2 <- apply(Z, 2, function(x){
                    mod <- KFAS::SSModel(x ~ - 1 + SSMcustom(Z=Zt, T=Tt, R=Rt, Q=Qt, a1=rep(0, ncol(Zt)), P1=matrix(0, ncol(Tt), ncol(Tt))), H=Ht)
                    v <- KFS(mod, filtering = "state", smoothing = "none")$v
                    return(v)
                })
            }else {
                v2 <- apply(Z, 2, function(x){
                    mod <- KFAS::SSModel(x ~ - 1 + SSMcustom(Z=Zt, T=Tt, R=Rt, Q=Qt, a1=rep(0, ncol(Zt)), P1=matrix(0, ncol(Tt), ncol(Tt))), H=Ht)
                    v <- KFS(mod, filtering = "state", smoothing = "none")$v
                    return(v)
                })
            } 
        }
        regmod <- lm(v ~ -1 + v2)
        mu <- summary(regmod)$coefficients[, 1]
        v <- residuals(regmod)
        
        mod <- KFAS::SSModel(y - Z%*%matrix(mu, nrow=length(mu))~ - 1 + SSMcustom(Z=Zt, T=Tt, R=Rt, Q=Qt, a1=rep(0, ncol(Zt)), P1=matrix(0, ncol(Tt), ncol(Tt))), H=Ht)
    }
    
    if(return.det) return(mu)
    if(return.mod) return(mod)
    if(!ll){
        v <- KFS(mod, filtering = "state", smoothing = "none")$v
        return(sum(v[START:n]^2))
    }
    ll  <- -logLik(mod)
}






UC_opt_KF_i1_ML <- function(theta, y, nulim = c(0.05, 10), quiet = TRUE, 
                         return.mod = FALSE, ll = TRUE, corr = FALSE, START=1,
                         deterministics = FALSE){
    n <- length(y)
    if(!corr){
        nu <- c(exp(theta[1]), exp(theta[2]))
        Qt <- matrix(nu[1], 1, 1)
        Tt <- matrix(1, 1, 1)
        Zt <- matrix(1, 1, 1)
        Rt <- matrix(1, ncol = 1)
        
        if(nu[1] > nulim[2] | nu[1] < nulim[1])  return(.Machine$integer.max)
        if(nu[2] > nulim[2] | nu[2] < nulim[1])  return(.Machine$integer.max)
        
        
        if(length(theta)>2){
            ar <- theta[-(1:2)]
            if(length(ar)> 1){
                Tt_c <- rbind(ar, cbind(diag(length(ar)-1), 0))
            }else{
                Tt_c <- matrix(ar, 1, 1)
            }
            Tt <- bdiag(Tt, Tt_c)
            Zt <- cbind(Zt, matrix(c(1, rep(0, length(ar)-1)), nrow = 1))
            Ht <- matrix(0, 1, 1) 
            Rt <- bdiag(Rt, c(1, rep(0, length(ar)-1)))
            Qt <- bdiag(c(Qt, nu[2]))
        }else{
            ar <- NULL
            Ht <- matrix(nu[2], 1, 1) 
        }
        
    }else{
        Qt <- mlogvech2mat(theta[1:3])
        nu <- c(Qt[1,1], Qt[2,2], Qt[1,2])
        Tt <- matrix(1, 1, 1)
        Zt <- matrix(1, 1, 1)
        Rt <- matrix(1, ncol = 1)
        # mod
        Ht <- matrix(0, 1, 1) 
        
        
        if(length(theta)>3){
            ar <- theta[-(1:3)]
            if(length(ar) > 1){
                Tt_c <- rbind(ar, cbind(diag(length(ar)-1), 0))
            }else{
                Tt_c <- matrix(ar, 1, 1)
            }
            Tt <- bdiag(Tt, Tt_c)
            Zt <- cbind(Zt, matrix(c(1, rep(0, length(ar)-1)), nrow = 1))
            Rt <- bdiag(Rt, c(1, rep(0, length(ar)-1)))
        }else{
            ar <- NULL
            Zt <- cbind(Zt, matrix(1, nrow =1, ncol=1))
            Tt <- bdiag(Tt, matrix(0, 1, 1))
            Rt <- bdiag(Rt, 1)
        }
    }
    
    #checks
    if(!is.null(ar)){
        if(!toComp(ar)$stable) return(.Machine$integer.max)
    }
    
    mod <- KFAS::SSModel(y ~ - 1 + SSMcustom(Z=Zt, T=Tt, R=Rt, Q=Qt, a1=rep(0, ncol(Zt)), P1=matrix(0, ncol(Tt), ncol(Tt))), H=Ht)
    
    
    
    if (is.logical(deterministics)) {
        if (deterministics) {
            det.type <- "trend"
        }else {
            det.type <- "none"
        }
    }else {
        det.type <- deterministics
        deterministics = TRUE
    }
    if (deterministics) {
        trendlist <- list()
        trendlist[[1]] <- trendlist[[2]] <- 0
        mod <- KFAS::SSModel(y ~ - 1 + SSMtrend(degree = 2, Q=trendlist)
                             +SSMcustom(Z=Zt, T=Tt, R=Rt, Q=Qt, a1=rep(0, ncol(Zt)), P1=matrix(0, ncol(Tt), ncol(Tt))), H=Ht)
    }else{
        mod <- KFAS::SSModel(y ~ - 1 +SSMcustom(Z=Zt, T=Tt, R=Rt, Q=Qt, a1=rep(0, ncol(Zt)), P1=matrix(0, ncol(Tt), ncol(Tt))), H=Ht)
    }
    
    if(return.mod) return(mod)
    if(!ll){
        v <- KFS(mod, filtering = "state", smoothing = "none")$v
        return(sum(v[START:n]^2))
    }else{
        ll  <- -logLik(mod)
        return(ll)
    }
    
}





UC_opt_KF_i2 <- function(theta, y, nulim = c(0.05, 10), quiet = TRUE, 
                         return.mod = FALSE, ll = TRUE, corr = FALSE, START=1,
                         deterministics = FALSE, return.det = FALSE){
    n <- length(y)
    if(!corr){
        nu <- exp(theta[1])
        Qt <- (matrix(1, 1, 1))
        Tt <- matrix(c(2, -1, 1, 0), 2, 2, byrow = TRUE)
        Zt <- matrix(c(1, 0), 1, 2)
        Rt <- matrix(c(1, 0), ncol = 1)
        
        if(nu[1] > nulim[2] | nu[1] < nulim[1])  return(.Machine$integer.max)
        
        if(length(theta)>1){
            ar <- theta[-(1)]
            if(length(ar)> 1){
                Tt_c <- rbind(ar, cbind(diag(length(ar)-1), 0))
            }else{
                Tt_c <- matrix(ar, 1, 1)
            }
            Tt <- bdiag(Tt, Tt_c)
            Zt <- cbind(Zt, matrix(c(1, rep(0, length(ar)-1)), nrow = 1))
            Ht <- matrix(0, 1, 1) 
            Rt <- bdiag(Rt, c(1, rep(0, length(ar)-1)))
            Qt <- bdiag(Qt, nu)
        }else{
            ar <- NULL
            Ht <- matrix(nu, 1, 1) 
        }
        
    }else{
        Qt <- matrix(c(1, theta[2], theta[2], theta[1]), 2, 2)
        
        if (theta[1] < nulim[1] | theta[1] > nulim[2]) 
            return(.Machine$integer.max)
        corr <- theta[2]/sqrt(theta[1])
        if (corr > 1 | corr < -1) 
            return(.Machine$integer.max)
        nu <- c(Qt[1, 1], Qt[2, 2], Qt[1, 2])
        Tt <- matrix(c(2, -1, 1, 0), 2, 2, byrow = TRUE)
        Zt <- matrix(c(1, 0), 1, 2)
        Rt <- matrix(c(1, 0), ncol = 1)
        # mod
        Ht <- matrix(0, 1, 1) 
        
        
        if(length(theta)>2){
            ar <- theta[-(1:2)]
            if(length(ar) > 1){
                Tt_c <- rbind(ar, cbind(diag(length(ar)-1), 0))
            }else{
                Tt_c <- matrix(ar, 1, 1)
            }
            Tt <- bdiag(Tt, Tt_c)
            Zt <- cbind(Zt, matrix(c(1, rep(0, length(ar)-1)), nrow = 1))
            Rt <- bdiag(Rt, c(1, rep(0, length(ar)-1)))
        }else{
            ar <- NULL
            Zt <- cbind(Zt, matrix(1, nrow =1, ncol=1))
            Tt <- bdiag(Tt, matrix(0, 1, 1))
            Rt <- bdiag(Rt, 1)
        }
    }
    
    #checks
    if(!is.null(ar)){
        if(!toComp(ar)$stable) return(.Machine$integer.max)
    }
    
    mod <- KFAS::SSModel(y ~ - 1 + SSMcustom(Z=Zt, T=Tt, R=Rt, Q=Qt, a1=rep(0, ncol(Zt)), P1=matrix(0, ncol(Tt), ncol(Tt))), H=Ht)
    
    
    
    if (is.logical(deterministics)) {
        if (deterministics) {
            det.type <- "trend"
        }else {
            det.type <- "none"
        }
    }else {
        det.type <- deterministics
        deterministics = TRUE
    }
    if (deterministics) {
        v <- KFS(mod, filtering = "state", smoothing = "none")$v
        if (det.type == "frac") {
            Z <- cbind(1, frac_diff(rep(1, length(y)), -1))
            if (!corr) {
                v2 <- apply(Z, 2, function(x) fUC_comp(x, 1, nu, ar)$v)
            }else{
                v2 <- apply(Z, 2, function(x) fUC_comp(x, 1, 
                                                       Q, ar, corr)$v)
            }
            regmod <- lm(v ~ -1 + v2)
            mu <- summary(regmod)$coefficients[, 1]
            v <- residuals(regmod)
        }else {
            if (det.type == "const") {
                Z <- matrix(1, n, 1)
            }
            if (det.type == "trend") {
                Z <- cbind(1, 1:n)
            }
            if (!corr) {
                v2 <- apply(Z, 2, function(x){
                    mod <- KFAS::SSModel(x ~ - 1 + SSMcustom(Z=Zt, T=Tt, R=Rt, Q=Qt, a1=rep(0, ncol(Zt)), P1=matrix(0, ncol(Tt), ncol(Tt))), H=Ht)
                    v <- KFS(mod, filtering = "state", smoothing = "none")$v
                    return(v)
                })
            }else {
                v2 <- apply(Z, 2, function(x){
                    mod <- KFAS::SSModel(x ~ - 1 + SSMcustom(Z=Zt, T=Tt, R=Rt, Q=Qt, a1=rep(0, ncol(Zt)), P1=matrix(0, ncol(Tt), ncol(Tt))), H=Ht)
                    v <- KFS(mod, filtering = "state", smoothing = "none")$v
                    return(v)
                })
            } 
        }
        regmod <- lm(v ~ -1 + v2)
        mu <- summary(regmod)$coefficients[, 1]
        if(length(mu)!=ncol(Z)) return(.Machine$integer.max)
        v <- residuals(regmod)
        
        mod <- KFAS::SSModel(y - Z%*%matrix(mu, nrow=length(mu))~ - 1 + SSMcustom(Z=Zt, T=Tt, R=Rt, Q=Qt, a1=rep(0, ncol(Zt)), P1=matrix(0, ncol(Tt), ncol(Tt))), H=Ht)
    }
    
    if(return.det) return(mu)
    if(return.mod) return(mod)
    if(!ll){
        v <- KFS(mod, filtering = "state", smoothing = "none")$v
        return(sum(v[START:n]^2))
    }
    ll  <- logLik(mod)
}








UC_opt_KF_i2_ML <- function(theta, y, nulim = c(0.05, 10), quiet = TRUE, 
                            return.mod = FALSE, ll = TRUE, corr = FALSE, START=1,
                            deterministics = FALSE){
    n <- length(y)
    if(!corr){
        nu <- c(exp(theta[1], exp(theta[2])))
        Qt <- (matrix(nu[1], 1, 1))
        Tt <- matrix(c(2, -1, 1, 0), 2, 2, byrow = TRUE)
        Zt <- matrix(c(1, 0), 1, 2)
        Rt <- matrix(c(1, 0), ncol = 1)
        
        if(nu[1] > nulim[2] | nu[1] < nulim[1])  return(.Machine$integer.max)
        if(nu[2] > nulim[2] | nu[2] < nulim[1])  return(.Machine$integer.max)
        
        
        if(length(theta)>2){
            ar <- theta[-(1:2)]
            if(length(ar)> 1){
                Tt_c <- rbind(ar, cbind(diag(length(ar)-1), 0))
            }else{
                Tt_c <- matrix(ar, 1, 1)
            }
            Tt <- bdiag(Tt, Tt_c)
            Zt <- cbind(Zt, matrix(c(1, rep(0, length(ar)-1)), nrow = 1))
            Ht <- matrix(0, 1, 1) 
            Rt <- bdiag(Rt, c(1, rep(0, length(ar)-1)))
            Qt <- bdiag(c(Qt, nu[2]))
        }else{
            ar <- NULL
            Ht <- matrix(nu[2], 1, 1) 
        }
        
    }else{
        Qt <- mlogvech2mat(theta[1:3])
        nu <- c(Qt[1,1], Qt[2,2], Qt[1,2])
        
        Tt <- matrix(c(2, -1, 1, 0), 2, 2, byrow = TRUE)
        Zt <- matrix(c(1, 0), 1, 2)
        Rt <- matrix(c(1, 0), ncol = 1)
        # mod
        Ht <- matrix(0, 1, 1) 
        
        
        
        
        
        if(length(theta)>3){
            ar <- theta[-(1:3)]
            if(length(ar) > 1){
                Tt_c <- rbind(ar, cbind(diag(length(ar)-1), 0))
            }else{
                Tt_c <- matrix(ar, 1, 1)
            }
            Tt <- bdiag(Tt, Tt_c)
            Zt <- cbind(Zt, matrix(c(1, rep(0, length(ar)-1)), nrow = 1))
            Rt <- bdiag(Rt, c(1, rep(0, length(ar)-1)))
        }else{
            ar <- NULL
            Zt <- cbind(Zt, matrix(1, nrow =1, ncol=1))
            Tt <- bdiag(Tt, matrix(0, 1, 1))
            Rt <- bdiag(Rt, 1)
        }
    }
    
    #checks
    if(!is.null(ar)){
        if(!toComp(ar)$stable) return(.Machine$integer.max)
    }
    
    
    if (is.logical(deterministics)) {
        if (deterministics) {
            det.type <- "trend"
        }else {
            det.type <- "none"
        }
    }else {
        det.type <- deterministics
        deterministics = TRUE
    }
    if (deterministics) {
        # Build trend
        #Tt_det <- matrix(c(1, 1, 1, 0), 2, 2)
        #Zt_det <- matrix(c(1, 0), nrow=1)
        #Rt_det <- matrix(0, nrow = 2, ncol = ncol(Rt))
        
        #Tt <- bdiag(Tt, Tt_det) 
        #Zt <- cbind(Zt, Zt_det)
        #Rt <- rbind(Rt, Rt_det)
        trendlist <- list()
        trendlist[[1]] <- trendlist[[2]] <- 0
        mod <- KFAS::SSModel(y ~ - 1 + SSMtrend(degree=2, Q=trendlist)+
                                 SSMcustom(Z=Zt, T=Tt, R=Rt, Q=Qt, a1=rep(0, ncol(Zt)), P1=matrix(0, ncol(Tt), ncol(Tt))), H=Ht)
        
        
    }else{
        mod <- KFAS::SSModel(y ~ - 1 + SSMcustom(Z=Zt, T=Tt, R=Rt, Q=Qt, a1=rep(0, ncol(Zt)), P1=matrix(0, ncol(Tt), ncol(Tt))), H=Ht)
        
    }

    if(return.mod) return(mod)
    if(!ll){
        v <- KFS(mod, filtering = "state", smoothing = "none")$v
        return(sum(v[START:n]^2))
    }else{
        ll  <- -logLik(mod)
        return(ll)
    }
}

