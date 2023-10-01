# Simulation:
# Performance of fUC for parameter estimation
#   - High and normal signal to noise ratios
#   - AR plus non-AR
#   - compare d with other nonparametric estimators
gc()
rm(list = ls())
library(fUCpack)
library(dplyr)
library(parallel)
library(KFAS)
devtools::load_all()

# Wd, etc
#source("./help functions/help_functions.R")
#source("./help functions/KF.R")

# Settings
n <- 100
d <- 1.5
r <- 1

# Strong persistence
ar <- c(-1.6, 0.8)


set.seed(42)
# Check 1: non-correlated UC vs. correlated UC
Q <- diag(2)
Q[2,2] <- 6
Q[2,1] <- Q[1,2] <- 1
U <- sapply(1:1, function(x) mvtnorm::rmvnorm(n, mean=c(0, 0), sigma = Q))
x <- frac_diff_multi(U[1:n,], d=-d)
u <- U[-(1:n),]
c <- stats::filter(c(rep(0, length(ar)), u), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))]
y <- x + c + 1:n

plot(y, type="l")
lines(x, col="2")

theta <- c(d, mlogvech(Q), ar)

#Q[1,2] <- Q[2,1] <- 0.5

# correlated filter
x.hat1 <- fUC_comp(y, d, 1, ar, corr=FALSE)
x.hat2 <- fUC_comp(y, d, Q, ar, corr=TRUE)
fUC_opt(c(d, mlogvech(Q), ar), y, deterministics = TRUE, corr=TRUE, START=1)
S0tP <- embed0(ma_inf(d[1], n-1) , n)
(t(S0tP) + t(S0cP))[-n, -n] %*% x.hat2$x[-1]
cbind(x.hat1$x, x.hat2$x, x.hat1$c, x.hat2$c) #  OK


fUC_opt_ML(c(d, log(1/5), ar), y, nu.opt=TRUE)



if(length(ar) < n-1){
    b <- c(ar, rep(0, n-length(ar)-1))
}else{
    b <- ar
}
S0cP  <- embed0(c(1, b), n)
(t(S0tP)*Q[2,1] + t(S0cP)*Q[1,1])[1:(n-1), 1:(n-1)] %*% (S0cP[1:(n-1), 1:(n-1)])%*%y[1:(n-1)]

(t(S0tP)*Q[2,1] + t(S0cP)*Q[1,1])[1:(n-2), 1:(n-2)] %*% (S0cP[1:(n-2), 1:(n-2)])%*%y[1:(n-2)]


# correlated filter vs. state space model
Tt <- rbind(-frac_diff(c(1, rep(0, length(y))), d)[-1], cbind(diag(length(y) - 1), 0))
Tt_c <- rbind(-ar, cbind(diag(length(ar) - 1), 0))
Tt <- bdiag(Tt, Tt_c)
Zt <- matrix(c(1, rep(0, length(y) - 1)), ncol = length(y))
Zt <- cbind(Zt, matrix(c(1, rep(0, length(ar) - 1)), nrow = 1))
Rt <- matrix(c(1, rep(0, NROW(Tt) - 1)), ncol = 1)
Rt <- cbind(Rt, c(rep(0, length(y)), 1, rep(0, length(ar) - 1)))
P_v <- matrix(NA, length(y), 1)
Pt <- matrix(0, NCOL(Tt), NCOL(Tt))
Pt[1, 1] <- Q[1, 1]
Pt[length(y) + 1, length(y) + 1] <- Q[2, 2]
Pt[1, length(y) + 1] <- Pt[length(y) + 1, 1] <- Q[1,2]
P1 <- matrix(0, n+2, n+2)
SSMod <- SSModel(y ~ -1 + SSMcustom(Zt, Tt, Rt, Q, a1=rep(0, n+2), P1 = Pt), H = matrix(0, 1, 1))
xhat.kf <- KFS(SSMod, filtering = "state", smoothing = "none")
cbind(xhat.kf$a[-(n+1),c(1, n+1)], x.hat2$x, x.hat2$c)

plot(xhat.kf$a[-(n+1), 1], type="l")
lines(x.hat2$x, col="2")


# by hand:
TV_inv_hand <- function(y, d, ar, Q){
    if(length(ar) < n-1){
        b <- c(ar, rep(0, n-length(ar)-1))
    }else{
        b <- ar
    }
    n <- length(y)
    Sd <- embed0(ma_inf(d[1], n-1) , n)
    B  <- embed0(c(1, b), n)
    ma  <- ma_inf(d[1],NROW(y))[-1]
    ma_c  <- c(ar, rep(0, n-length(ar)))
    x_hat <- c_hat <- rep(0, n)
    for(t in 1:(n-1)){
        x_hat_t <- solve(Q[2,2]*crossprod(Sd[1:t, 1:t]) + Q[1,1]*crossprod(B[1:t, 1:t])
                         + Q[1,2]*crossprod(Sd[1:t, 1:t], B[1:t, 1:t]) + Q[1,2]*crossprod(B[1:t, 1:t], Sd[1:t, 1:t])) %*%
            (crossprod(B[1:t, 1:t])*Q[1,1] + crossprod(Sd[1:t, 1:t], B[1:t, 1:t])*Q[1,2]) %*% y[1:t]
        c_hat_t <- y[1:t] - x_hat_t
        x_hat[t+1] <- -ma[t:1]%*%x_hat_t
        c_hat[t+1] <- -ma_c[t:1]%*%c_hat_t
    }
    res <- cbind(x_hat, c_hat, y - x_hat - c_hat)
    colnames(res) <- c("x", "c", "v")
    return(res)
}
TV_inv_hand(y, d, ar, Q)

S0tP <- embed0(ma_inf(d[1], n-1) , n)
if(length(ar) < n-1){
    b <- c(ar, rep(0, n-length(ar)-1))
}else{
    b <- ar
}
S0cP  <- embed0(c(1, b), n)
Xacc <- matrix(0, n, n)
Bacc <- matrix(0, n, n)
a <- Q[1,1]
b <- Q[2,2]
c <- Q[1,2]

for (t in 1:n){
    k = t
    for(i in 1:k){
        for(j in i:k){
            tempA = 0.0
            tempB = 0.0
            tempD = 0.0
            tempE = 0.0

            if( i == 1){
                tempA = tempA + S0tP[i] * S0tP[j]
                tempB = tempB + S0cP[i] * S0cP[j]
                tempD = tempD + S0tP[i] * S0cP[j]
                tempE = tempE + S0cP[i] * S0tP[j]

                Xacc[(j-1)*n + i] = tempA * b + tempB * a + (tempD + tempE)*c
                Bacc[(j-1)*n + i] = tempB * a + tempD * c

            }else{
                Xacc[(j-1)*n + i] = Xacc[(j-2)*n + i - 1] + S0tP[i] * S0tP[j] * b +
                    S0cP[i] * S0cP[j] * a + (S0tP[i] * S0cP[j] + S0cP[i] * S0tP[j]) * c
            }
        }
    }
}

            for (i = 0; i < k; i++){
                for (j = i; j < k; j++){

                    tempA = 0.0;
                    tempB = 0.0;
                    tempD = 0.0;
                    tempE = 0.0;


                    if (i == 0){

                        /* in this case we have to compute the full row */
                            tempA += S0tP[ i ] * S0tP[ j ];
                            tempB += S0cP[ i ] * S0cP[ j ];
                            tempD += S0tP[ i ] * S0cP[ j ];
                            tempE += S0cP[ i ] * S0tP[ j ];


                            /* using triangular structure */
                                Xacc[ j * N + i] = tempA * b + tempB * a + (tempD + tempE) * c;
                                Bacc[ j * N + i] = tempB * a + tempD * c;
                                /*Bacc[ j * N + i] = tempB * a + tempE * c;*/


                                    /* if (t<=tCheck) { */
                                            /*   Rprintf("AA'[%i,%i]: %f\nIndex: %i\n", i,j, temp * a, */
                                                             /*           m * (N + 1) + j*N + i); */
                                            /* } */

                                    /* if (t==tCheck) { */
                                            /*   Rprintf("Xacc[%i,%i]: %f\n", i,j, Xacc[m * (N + 1) + j * N + i]); */
                                            /* } */

                    } else {

                        /* in this case, we only have to add the new stuff to the old one*/

                            //if (t==2) for (ii = 0; ii < k; ii++) Rprintf("A[%i,%i]: %f", i, ii, S0tP[m * (N + 1 + (k - ii - 1)) + i]);


                        /* if (t<=tCheck) { */
                                /*   Rprintf("AApre'[%i,%i]: %f\nIndex: %i\n", i,j, Xacc[m * (N + 1) + j*N + i], */
                                                 /*           m * (N + 1) + j*N + i); */
                                /* } */


                            Xacc[ j * N + i] = Xacc[ (j-1) * N + i - 1] + S0tP[i] * S0tP[j] * b + S0cP[i] * S0cP[j] * a + (S0cP[i] * S0tP[j] + S0tP[i] * S0cP[j])*c;
                            Bacc[ j * N + i] = Bacc[ (j-1) * N + i - 1] + S0cP[i] * S0cP[j] * a + S0tP[i] * S0cP[j] * c;
                            /*Bacc[ j * N + i] = Bacc[ (j-1) * N + i - 1] + S0cP[i] * S0cP[j] * a + S0cP[i] * S0tP[j] * c;*/



                                /* if (t==tCheck) { */
                                        /*   Rprintf("Xacc[%i,%i]: %f\n", i,j, Xacc[m * (N + 1) + j*N + i]); */
                                        /* } */

                    }
                }
            }
