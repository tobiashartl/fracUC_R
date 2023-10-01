fUC_comp_old <- function(y, d, nu, ar){


    n <- NROW(y)
    if(is.null(ar)) ar <- 0

    # Generate coefficient matrices
    S0t      <- embed0(ma_inf(-d[1], n-1) , n)
    S0c      <- embed0(c(1, ARMAtoMA(ar =- ar, lag.max=n-1)), n)
    # Adjust for signal-to-noise ratio
    Q        <- diag(2)
    Q[2,2]   <- nu
    # Obtain coefficients for x, c
    ma       <- ma_inf(-d[1],NROW(y))[-1]
    ma_c     <- ARMAtoMA(ar = -ar, lag.max=n)

    # calculate Kalman-Filter
    ### Note: We use an algebraic identity that is computationally superior
    TC <- .Call("TC_inv",
                S0t = S0t,
                S0c = S0c,
                Q = Q,
                y = y,
                n = as.integer(n),
                ma = ma,
                ma_c = ma_c,
                PACKAGE = "TCinv")

    # return trend, cycle, prediction error
    return(list(x = TC[,1],
                c = TC[,2],
                v = y - TC[,1] - TC[,2]))
}



fUC_comp2_old <- function(y, d, nu, ar){


    n <- NROW(y)
    if(is.null(ar)) ar <- 0

    if(length(ar) < n-1){
        b <- c(ar, rep(0, n-length(ar)-1))
    }else{
        b <- ar
    }
    # Build Sd
    Sd <- t(embed0(matrix( ma_inf(d, n-1), ncol = 1), n))
    B  <- t(embed0(matrix(c(1, b), ncol=1), n))
    s  <- ma_inf(d, n)[-1]

    SS <- crossprod(Sd)
    BB <- crossprod(B)

    X  <- SS + nu*BB
    # calculate Kalman-Filter
    ### Note: We use an algebraic identity that is computationally superior
    # TC <- .Call("TC_inv",
    #             S0t = S0t,
    #             S0c = S0c,
    #             Q = Q,
    #             y = y,
    #             n = as.integer(n),
    #             ma = ma,
    #             ma_c = ma_c,
    #             PACKAGE = "TCinv")
    inv.fn <- function(t){
        if(t == 1){
            return(c(0, 0))
        }else{
            xhat <- -s[1:(t-1)]%*%solve(X[1:(t-1), 1:(t-1)], BB[1:(t-1), 1:(t-1)]%*%y[(t-1):1])
            chat <- -b[1:(t-1)]%*%solve(X[1:(t-1), 1:(t-1)], BB[1:(t-1), 1:(t-1)]%*%y[(t-1):1])
            return(c(xhat, chat))
        }
    }
    TC <-  t(sapply(1:n, inv.fn))



    # return trend, cycle, prediction error
    return(list(x = TC[,1],
                c = TC[,2],
                v = y - TC[,1] - TC[,2]))
}


fUC_opt_old <- function(theta, y, nulim = c(0.05, 10), quiet = TRUE){
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
    if(!quiet) cat("ll = ", sum(v^2), " , theta = ", theta, "\n")
    return(sum(v^2))
}





fUC_smooth_old <- function(y, d, nu, ar){


    n <- NROW(y)
    if(is.null(ar)) ar <- 0

    # Generate coefficient matrices
    S0t      <- embed0(ma_inf(-d[1], n-1) , n)
    S0c      <- embed0(c(1, ARMAtoMA(ar =- ar, lag.max=n-1)), n)
    # Adjust for signal-to-noise ratio
    Q        <- diag(2)
    Q[2,2]   <- nu
    # Obtain coefficients for x, c
    ma       <- ma_inf(-d[1],NROW(y))[-1]
    ma_c     <- ARMAtoMA(ar = -ar, lag.max=n)

    # calculate Kalman-Filter
    ### Note: We use an algebraic identity that is computationally superior
    Sigma_y <- solve(tcrossprod(S0t[1:n, n:1], S0t[1:n, n:1]*Q[1,1]) +
                         tcrossprod(S0c[1:n, n:1], S0c[1:n, n:1]*Q[2,2]),
                     y[1:n])

    tau_t <- S0t%*%(crossprod(S0t[1:n, 1:n]*Q[1,1], Sigma_y))
    c_t <- S0c%*%crossprod(S0c[1:n, 1:n]*Q[2,2], Sigma_y)


    TC  <- cbind(tau_t, c_t)

    # return trend, cycle, prediction error
    return(list(x = TC[,1],
                c = TC[,2],
                v = y - TC[,1] - TC[,2]))
}




UC_opt_old <- function(theta, y, nulim = c(0.05, 10), quiet = TRUE){
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
        if(!toComp(ar)$stable) return(.Machine$integer.max)
    }
    if(nu > nulim[2] | nu < nulim[1])  return(.Machine$integer.max)

    v <- fUC_comp(y, d, nu, ar)$v
    if(!quiet) cat("ll = ", sum(v^2), " , theta = ", theta, "\n")
    return(sum(v^2))
}
