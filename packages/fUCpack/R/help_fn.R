# Coefficients of fractional difference operator
#' @export
ma_inf <- function(d, N) c(1,cumprod(((1:N)-1-d)/(1:N)))







#' @export
bdiag <- function (...)
{
    if (nargs() == 1)  x <- as.list(...)
    else x <- list(...)
    n <- length(x)
    if (n == 0) return(NULL)
    x <- lapply(x, function(y) if (length(y)) as.matrix(y)
                else stop("Zero-length component in x"))
    d <- array(unlist(lapply(x, dim)), c(2, n))
    rr <- d[1, ]
    cc <- d[2, ]
    rsum <- sum(rr)
    csum <- sum(cc)
    out <- array(0, c(rsum, csum))
    ind <- array(0, c(4, n))
    rcum <- cumsum(rr)
    ccum <- cumsum(cc)
    ind[1, -1] <- rcum[-n]
    ind[2, ] <- rcum
    ind[3, -1] <- ccum[-n]
    ind[4, ] <- ccum
    imat <- array(1:(rsum * csum), c(rsum, csum))
    iuse <- apply(ind, 2, function(y, imat) imat[(y[1] + 1):y[2],
                                                 (y[3] + 1):y[4]], imat = imat)
    iuse <- as.vector(unlist(iuse))
    out[iuse] <- unlist(x)
    return(out)
}


#' @export
EW <- function (x, m = NULL, alpha = 0.5, type = 2, interval = c(0,
                                                                 1))
{
    x <- na.omit(x)
    EW_LL <- function(d, m, x, lambda, n, type) {
        if (is.numeric(type))
            x <- demeanFracEW(x, d, type)
        I_d <- base::Re((stats::fft(frac_diff(x, d)) * base::Conj(stats::fft(frac_diff(x,
                                                                                       d))))[2:(m + 1)]/2/pi/n)
        return(log(mean(I_d)) - (2 * d) * mean(log(lambda)))
    }
    n <- length(x)
    if (is.null(m))
        m <- trunc(n^(alpha))
    d_hat <- optimize(EW_LL, m = m, n = length(x), lambda = (1:m) *
                          2 * pi/n, x = na.omit(x), type = type, interval = interval)$minimum
    return(d_hat)
}

#' @export
demeanFracEW <- function (z, d, type)
{
    x <- as.matrix(z)
    n <- NROW(x)
    k <- NCOL(x)
    switch(type, {
        dim(x) <- dim(z)
        y <- x
        return(y)
    }, {
        y <- x
        for (j in 1:k) {
            if (d[j] < 0.5) y[, j] <- x[, j] - mean(x[, j])
            if (d[j] > 0.75) y[, j] <- x[, j] - x[1, j]
            if (d[j] >= 0.5 & d[j] <= 0.75) {
                y[, j] <- x[, j] - 0.5 * (1 + cos(4 * pi * d[j])) *
                    mean(x[, j]) - (1 - 0.5 * (1 + cos(4 * pi *
                                                           d[j]))) * x[1, j]
            }
        }
        dim(y) <- dim(z)
        return(y)
    }, {
        Xhat <- as.matrix(lm(x ~ I(1:n))$resid)
        y <- x
        for (j in 1:k) {
            if (d[j] < 0.5) y[, j] <- Xhat[, j] - mean(Xhat[,
                                                            j])
            if (d[j] > 0.75) y[, j] <- Xhat[, j] - Xhat[1, j]
            if (d[j] >= 0.5 & d[j] <= 0.75) {
                y[, j] <- Xhat[, j] - 0.5 * (1 + cos(4 * pi *
                                                         d[j])) * mean(Xhat[, j]) - (1 - 0.5 * (1 +
                                                                                                    cos(4 * pi * d[j]))) * Xhat[1, j]
            }
        }
        dim(y) <- dim(z)
        return(y)
    })
}

#' @export
frac_demean <- function(y, d){
    xtilde <- frac_diff_multi(y, d)
    #seas_mat <- matrix((rep(diag(7), ceiling(n/7))), ncol=7, byrow = T)[1:n, ]
    #seas_mat_frac <- frac_diff_multi(seas_mat, d)
    c  <- rep(1, length(y))
    con  <- frac_diff(c, d)
    mu <- lm(xtilde ~ -1 + con)$coef
    y  <- y - c*mu
}

#' @title frac_diff_multi
#' @description Takes the fractional difference of a matrix of inputs (column-wise)
#' @param y A (T x N) matrix of N time series
#' @param d Order of differencing

#' @return The d-th difference of y
#'
#' @export
frac_diff_multi <- function (y, d)
{
    data2 <- as.data.frame(y)
    data2 <- mapply(frac_diff, y = data2, d = d)
    if (is.ts(y))
        data2 <- as.ts(data2, start = start(y), frequency = frequency(y))
    if (is.matrix(y))
        data2 <- as.matrix(data2)
    return(data2)
}

#' @title frac_diff
#' @description Takes the fractional difference of a vector time series
#' @param y A vector time series
#' @param d Order of differencing

#' @return The d-th difference of y
#'
#' @export
frac_diff <- function (y, d)
{
    orig_data <- y
    if (any(is.na(y))) {
        warning("There are NAs in frac_diff")
        y <- na.omit(y)
    }
    N <- length(y)
    coef <- c(1, cumprod(((1:N) - 1 - d)/(1:N)))
    orig_data[!is.na(orig_data)] <- stats::filter(c(rep(0, N), y), filter = coef,
                                                  sides = 1)[-(1:N)]
    return(orig_data)
}


#' @export
embed0 <- function (x, dimension = 1, fill = 0){
    if (is.matrix(x)){
        x0 <- rbind(matrix(fill,dimension-1,NCOL(x)))
        embed(rbind(x0,x), dimension = dimension)
    } else
    {
        x0 <- rep(fill,dimension-1)
        embed(c(x0,x), dimension = dimension)
    }
}

#' @title toComp
#' @description Constructs companion form of AR coefficients
#' @param x AR coefficients

#' @return Companion form of the AR
#'
#' @export
toComp <- function(x)
{
    if (is.null(dim(x))) x <- matrix(x,1)
    k <- nrow(x)
    p <- ncol(x)/k
    if (p>1) x <- rbind(x, cbind(diag((p-1)*k),matrix(0,(p-1)*k,k)))
    eigv <- eigen(x,only.values=TRUE)$values
    return(list(CompMat=x,eigv=eigv,stable=all(round(abs(eigv), digits = 2)<1)))
}

#' @export
vech <- function(x) x[lower.tri(x,diag=TRUE)]


#' @export
vecl <- function(x) x[lower.tri(x)]


#' @export
vech2mat <- function(x)
{
    k <- which(length(x)==((1:500)*(2:501)/2))
    X <- matrix(0,k,k)
    X[lower.tri(X,diag=TRUE)] <- x
    X <- X+t(X)-diag(diag(X))
    X
}


#' @export
mlog <- function(A)
{
    if (length(A)==1) return(log(A))
    L <- diag(eigen(A)$values)
    V <- eigen(A)$vectors
    V%*%diag(log(diag(L)))%*%t(V)
}

#' @export
mexp <- function(A, eta = NULL)
{
    if (length(A)==1) return(exp(A))
    l <- exp(eigen(A)$values)
    if(!is.null(eta)){
        if(any(l < eta)){
            l[l<eta] <- eta
        }
    }
    L <- diag(l)
    V <- eigen(A)$vectors
    V%*%L%*%t(V)
}

#' @title mlogvech
#' @description Matrix log-transformation
#' @param x Matrix to be log-transformed

#' @return log-transformed matrix
#'
#' @export
mlogvech <- function(x) vech(mlog(x))

#' @title mlogvech2mat
#' @description Matrix inverse log-transform
#' @param x Log-transformed matrix

#' @return Inverse log-transformed matrix
#'
#' @export
mlogvech2mat <- function(x, eta=NULL)
{
    if (length(x)==1) return(matrix(exp(x),1,1))
    k <- which(length(x)==((1:500)*(2:501)/2))
    if (any(is.na(x))) return(matrix(NA,k,k))
    X <- matrix(0,k,k)
    X[lower.tri(X,diag=TRUE)] <- x
    X <- X+t(X)-diag(diag(X))
    mexp(X, eta)
}

#' @export
mlogvech2mat_std <- function(x, eta=NULL)
{
    x <- c(0, x)
    #if (length(x)==1) return(matrix(exp(x),1,1))
    k <- which(length(x)==((1:500)*(2:501)/2))
    if (any(is.na(x))) return(matrix(NA,k,k))
    X <- matrix(0,k,k)
    X[lower.tri(X,diag=TRUE)] <- x
    X <- X+t(X)-diag(diag(X))
    mexp(X, eta)
    X <- X
    return(X)
}


#' @title pd_cov
#' @description Different matrix transformations that ensure p.d.
#' @param theta Transformed parameters to be optimized over
#' @param type Type of transformation. Possible types are "mlv" (matrix exponential),
#' "chol" (cholesky transformation), pchol (transformewd cholesky transformation),
#' and "beta" (beta transformation). "mlv" is suggested

#' @return p.d. matrix
#'
#' @export
pd_cov <- function(theta, type = "mlv", flip=F, neg=F){

    # cov(n_t, e_t) p.d.
    if(type == "mlv") Q <- mlogvech2mat(theta)
    if(type == "chol") Q <- chol_cov(theta,flip,neg)
    if(type == "pchol") Q <- pchol_cov(theta)
    if(type == "beta") Q <- beta_cov(theta)

    return(Q)
}

#' @export
beta_cov <- function(theta){
    # entries in theta: sigma, beta, gamma
    sigma <- theta[1]
    beta <- theta[2]
    gamma <- theta[3]
    Q <- sigma^2 * matrix(c(1, beta, beta, beta^2 + gamma^2), 2, 2)
    return(Q)
}


#' @export
beta_cov_inv <- function(Q){
    # entries in theta: sigma, beta, gamma
    sigma <- sqrt(Q[1,1])
    beta  <- Q[1,2]/sigma^2
    gamma <- sqrt(Q[2,2]/sigma^2 - beta^2 )
    theta <- c(sigma, beta, gamma)
    return(theta)
}

#' @export
chol_cov <- function(theta, switch=FALSE, neg=FALSE){

    psi       <- matrix(0,2,2)
    if(neg){
        psi[1,1]  <- exp(-theta[1])
        if(!switch){
            psi[1,2]  <- theta[2]
        }else{
            psi[2,1] <- theta[2]
        }
        psi[2,2]  <- exp(-theta[3])
    }else{
        psi[1,1]  <- exp(theta[1])
        if(!switch){
            psi[1,2]  <- theta[2]
        }else{
            psi[2,1] <- theta[2]
        }
        psi[2,2]  <- exp(theta[3])
    }

    Q         <- t(psi)%*%psi
    return(Q)
}

#' @export
cov_to_par <- function(Q, neg = FALSE, flip=FALSE){

    # returns inital parameters
    if(!flip){
        psi <- chol(Q)
        #if(flip) psi <- t(psi)
        if(neg){
            theta <- c(-log(psi[1,1]), psi[1,2], -log(psi[2,2]))
        }else{
            theta <- c(log(psi[1,1]), psi[1,2], log(psi[2,2]))

        }
    }else{
        if(neg){
            theta3 <- -log(Q[2,2])/2
            theta2 <- Q[1,2]/exp(-theta3)
            theta1 <- -log(Q[1,1] - theta2^2)/2
        }else{
            theta3 <- log(Q[2,2])/2
            theta2 <- Q[1,2]/exp(theta3)
            theta1 <- log(Q[1,1] - theta2^2)/2
        }
        theta <- c(theta1, theta2, theta3)
    }

    return(theta)
}



#' @export
pchol_cov <- function(theta) {
    a <- theta[1]
    b <- theta[2]
    c <- theta[3]

    V <- matrix(c(exp(a), b, 0, sqrt(exp(2*c) - b^2)), 2, 2)
    Q <- t(V) %*% V

    return(Q)
}


#' @export
inv_pchol_cov <- function(Q) {
    a <- log(Q[1,1])
    b <- Q[1,2]
    c <- log(Q[2,2] - b^2)

    theta <- c(a, b, c)

    return(theta)
}



#' @export
trans <- function(c0)
{
    c1 = c0
    c1[1:2] = exp(-c0[1:2])

    # Constrain cov(n_t, e_t) to be pd using Cholesky factorization
    p       = matrix(0,2,2)
    p[1,1]  = exp(-c0[1])
    p[2,1]  = c0[3]
    p[2,2]  = exp(-c0[2])
    varcov  = t(p)%*%p
    c1[1]   = (varcov[1,1])
    c1[2]   = (varcov[2,2])
    c1[3]   = varcov[2,1]

    return(c1)

}

#' @export
trans2mat <- function(c0){
    c1 <- trans(c0)
    Q  <- matrix(NA, 2, 2)
    Q[1,1] <- c1[1]
    Q[2,2] <- c1[2]
    Q[1,2] <- Q[2,1] <- c1[3]
    return(Q)
}




