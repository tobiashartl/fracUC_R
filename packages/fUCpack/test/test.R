devtools::load_all()
gc()

n <- 100
d <- 1.25
nu <- 3
R <- 1

# Strong persistence
ar <- c(-1.6, 0.8)

# ==============================================================================
# Check: TC_inv:
if(length(ar) < n-1){
    b <- c(ar, rep(0, n-length(ar)-1))
}else{
    b <- ar
}

# simulate data
x <- frac_diff(rnorm(n), -d)
c <- stats::filter(c(rep(0, length(ar)), rnorm(n)), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))]

y <- x + c
# Generate coefficient matrices
S0t      <- embed0(ma_inf(d[1], n-1) , n)
S0c      <- embed0(c(1, b), n)
# Obtain coefficients for x, c
ma       <- ma_inf(d[1],NROW(y))[-1]
ma_c     <- c(ar, rep(0, n-length(ar)))

Q <- matrix(c(1, 0, 0, nu), 2, 2)
theta <- c(d, mlogvech(Q), ar)

TV_inv_hand <- function(y, d, ar, nu){
    n <- length(y)
    Sd <- embed0(ma_inf(d[1], n-1) , n)
    B  <- embed0(c(1, b), n)
    ma  <- ma_inf(d[1],NROW(y))[-1]
    ma_c  <- c(ar, rep(0, n-length(ar)))
    x_hat <- c_hat <- rep(0, n)
    for(t in 1:(n-1)){
        x_hat_t <- solve(nu*crossprod(Sd[1:t, 1:t]) + crossprod(B[1:t, 1:t])) %*% crossprod(B[1:t, 1:t]) %*% y[1:t]
        c_hat_t <- nu*solve(nu*crossprod(Sd[1:t, 1:t]) + crossprod(B[1:t, 1:t])) %*% crossprod(Sd[1:t, 1:t]) %*% y[1:t]
        #c_hat_t <- y[1:t] - x_hat_t
        x_hat[t+1] <- -ma[t:1]%*%x_hat_t
        c_hat[t+1] <- -ma_c[t:1]%*%c_hat_t
    }
    res <- cbind(x_hat, c_hat, y - x_hat - c_hat)
    colnames(res) <- c("x", "c", "v")
    return(res)
}

system.time(TC1 <- TV_inv_hand(y, d, ar, nu))
system.time(TC2 <- TCinv(S0t, S0c, nu, y, n , ma, ma_c)[-(n+1),])
system.time(TC3 <- TCinv_cor(S0t, S0c, sigma_eta = 1,
                             sigma_eps = 3, cov_etaeps=0, y, n , ma, ma_c)[-(n+1),])




cbind(TC1[,1]-TC2[,1], TC1[,2]-TC2[,2])
# checked: works perfectly well
plot(x, type="l")
lines(TC2[,1], col="2")


# ==============================================================================
# Check: TC_inv_cor:
corr <- 0
Q <- matrix(c(1, corr, corr, nu), 2, 2)
err <- mvtnorm::rmvnorm(n, mean = rep(0, 2), sigma = Q)
if(length(ar) < n-1){
    b <- c(ar, rep(0, n-length(ar)-1))
}else{
    b <- ar
}

# simulate data
x <- frac_diff(err[,1], -d)
c <- stats::filter(c(rep(0, length(ar)), err[,2]), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))]

y <- x + c
# Generate coefficient matrices
S0t      <- embed0(ma_inf(d[1], n-1) , n)
S0c      <- embed0(c(1, b), n)
# Obtain coefficients for x, c
ma       <- ma_inf(d[1],NROW(y))[-1]
ma_c     <- c(ar, rep(0, n-length(ar)))




TV_inv_hand_cor <- function(y, d, ar, nu){
    n <- length(y)
    Sd <- embed0(ma_inf(d[1], n-1) , n)
    B  <- embed0(c(1, b), n)
    ma  <- ma_inf(d[1],NROW(y))[-1]
    ma_c  <- c(ar, rep(0, n-length(ar)))
    x_hat <- c_hat <- rep(0, n)
    for(t in 1:(n-1)){
        x_hat_t <- solve(nu*crossprod(Sd[1:t, 1:t]) + crossprod(B[1:t, 1:t]) -
                             ) %*% crossprod(B[1:t, 1:t]) %*% y[1:t]
        c_hat_t <- nu*solve(nu*crossprod(Sd[1:t, 1:t]) + crossprod(B[1:t, 1:t])) %*% crossprod(Sd[1:t, 1:t]) %*% y[1:t]
        x_hat[t+1] <- -ma[t:1]%*%x_hat_t
        c_hat[t+1] <- -ma_c[t:1]%*%c_hat_t
    }
    res <- cbind(x_hat, c_hat, y - x_hat - c_hat)
    colnames(res) <- c("x", "c", "v")
    return(res)
}

TC1 <- TV_inv_hand(y, d, ar, nu)
TC2 <- TCinv(S0t, S0c, nu, y, n , ma, ma_c)[-(n+1),]
cbind(TC1[,1]-TC2[,1], TC1[,2]-TC2[,2])
# checked: works perfectly well
plot(x, type="l")
lines(TC2[,1], col="2")








set.seed(42)

# Generate the data
x <- frac_diff_multi(matrix(rnorm(n*R, mean=0, sd=1), n, R), -d)
u <- matrix(rnorm(n*R, mean=0, sd=sqrt(nu)), n, R)
c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
y <- x + c

    optfn <- function(n, y, x){
        # Estimate
        tryCatch({
            (est <- optim(par=c(1, 0, -0.5, 0.5),
                          fn = fUC_opt,
                          y=y, method = "BFGS", nulim = c(1/100, 100),
                          deterministics = FALSE))


            par <- est$par
            nu  <- exp(par[2])
            ar <- est$par[-(1:2)]
            KF <- fUC_comp(y, par[1], nu, ar)
            #KF2 <- fUC_comp((1:length(y))^par[1], par[1], nu, ar)
            #mu  <- summary(lm(KF$v ~ -1 + KF2$v))$coefficients[1,1]


            KS <- fUC_smooth(y, par[1], nu, ar = ar)
            Rsq <- summary(lm(x ~ KS$x))$r.squared


            # Return results
            results <- c(par[1], exp(par[2]), ar, Rsq)
            names(results) <- c("d", "nu", "ar_1", "ar_2", "Rsq")
            return(results)
        }, error = function(e) return(rep(NA, 5))
        )

    }
    ### fUC part
    cl <- makeCluster(6)
    clusterExport(cl, c("optfn", "fUC_opt", "y", "fUC_comp", "ma_inf",
                        "embed0", "fUC_smooth", "x", "frac_diff", "lm", "n"))
    clusterEvalQ(cl, library(fUCpack))
    RESULTS <- parSapply(cl, 1:R, function(j) optfn(n, y[, j], x[, j]))
    stopCluster(cl)

    dEW_45 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2), alpha = .45))
    dEW_50 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2), alpha = .50))
    dEW_55 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2), alpha = .55))
    dEW_60 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2), alpha = .60))
    dEW_65 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2), alpha = .65))
    dEW_70 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2), alpha = .70))

    RESULTS <- rbind(RESULTS, c(dEW_45), c(dEW_50), c(dEW_55), c(dEW_60), c(dEW_65), c(dEW_70))
    rownames(RESULTS) <- c("d", "nu", "ar_1", "ar_2", "Rsq",
                           "d_45", "d_50", "d_55", "d_60", "d_65", "d_70")
    save(RESULTS, file = paste("./MC/MC_2/Sim_R", R, "_n", n, "_d", d, "_nu", nu, "_ar_high.RData", sep=""))
}

names <- list.files("./MC/MC_2/")[seq(from = 2, to = 72, by = 2)]
results.list <- list()
statistics <- matrix(NA, 36, 13)
for(i in 1:36){
    true.par          <- setups[i, ]
    load(paste("./MC/MC_2/Sim_R1000_n", true.par[1], "_d", true.par[2], "_nu", true.par[3],
               "_ar_high.RData", sep=""))
    results.list[[i]] <- RESULTS
    MSE_d <- rowMeans((results.list[[i]][c("d","d_45", "d_50", "d_55", "d_60", "d_65", "d_70"), ] - true.par$d)^2)
    SSR   <- mean(results.list[[i]][c("SSR"), ])
    Rsq   <- mean(results.list[[i]][c("Rsq"), ])
    MSE_nu <- mean((results.list[[i]][c("nu"), ] - true.par$r)^2)
    MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(-0.7, 0.4))^2)
    MSE_sigma <- mean((results.list[[i]][c("sigma_eta"), ] - 1/true.par$r)^2)
    statistics[i, ] <- c(MSE_d, MSE_nu, MSE_sigma, MSE_ar, SSR, Rsq)
}


hist(RESULTS["nu",], breaks = 50)




n <- 300
d<- 1
nu <- 1
ar <- -.5

eta <- rnorm(n, 0, 1)
eps <- rnorm(n, 0, nu)
if(!is.null(ar)){
    c <- stats::filter(c(rep(0, length(ar)), eps), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))]
}else{
    c <- eps
}


x <- frac_diff(eta, -d)
y <- x + c

if(length(ar) < n-1){
    b <- c(ar, rep(0, n-length(ar)-1))
}else{
    b <- ar
}
# Build Sd
Sd <- (embed0(matrix( ma_inf(d, n-1), ncol = 1), n))
B  <- (embed0(matrix(c(1, b), ncol=1), n))
s  <- ma_inf(d, n)[-1]
b  <- c(ar, rep(0, n-1))


#TC <- TC_old(y, Sd, B, nu, s, b)
gc()
(TCnew <- TCinv(Sd, B, nu, y, n, s, b))
system.time(fUC_new <- fUC_comp(y, d, nu, ar))

cbind(fUC_comp(y, d, nu, ar)$c, c)
plot()

# check: FKF
library(FKF)
# Build FKF matrices and vectors
Tx <- rbind(-ma_inf(d, n)[-1], cbind(diag(n-1), 0))
Tc <- -ar
Tt <- bdiag(Tx, Tc)
Zt <- matrix(c(1, rep(0, NCOL(Tt)-1)), nrow=1)
a0 <- (rep(0, NCOL(Tt)))
P0 <- diag(length(a0)) * 0
dt <- matrix(rep(0, NCOL(Tt)), ncol=1)
ct <- matrix(0, 1, 1)
Hht <- diag(n+1)*0
Hht[1,1] <- 1/nu
Hht[n+1, n+1] <- 1
Ggt <- matrix(0, 1, 1)
yt <- matrix(y, nrow = 1)

fkfest <- fkf(a0, P0, dt, ct, Tt, Zt, Hht, Ggt, yt)
fkfest$vt
(TCnew[-(n+1),] -  cbind(fUC_new$x, fUC_new$c))




library(TCinv)
#setwd("/Users/tobias/Dokumente/Projekte/filtering unknown persistence/R/")
system.time(fUC_old <- fUC_comp_old(y, d, nu, ar))
fUC_old$c-fUC_new$c

smooth_new <- fUC_smooth(y, d, nu, ar)
smooth_old <- fUC_smooth_old(y, d, nu, ar)

smooth_new$x - smooth_old$x
smooth_new$c - smooth_old$c


# time check

system.time(TC <- fUC_comp_old(y, d, nu, ar))
system.time(TCnew <- fUC_comp(y, d, nu, ar))


solve(Sd[1:10, 1:10]%*%t(Sd[1:10, 1:10]) + B[1:10, 1:10]%*%t(B[1:10, 1:10]))%*%Sd[1:10, 1:10]%*%t(Sd[1:10, 1:10])%*%y[10:1]


cbind(solve(Sd%*%t(Sd) + B%*%t(B))%*%Sd%*%t(Sd)%*%y[n:1], TCnew[,2],
      solve(Sd%*%t(Sd) + B%*%t(B))%*%B%*%t(B)%*%y[n:1], TCnew[,1])

cbind(solve(Sd%*%t(Sd) + B%*%t(B))%*%Sd%*%t(Sd)%*%y[n:1],
      solve(Sd%*%t(Sd) + B%*%t(B))%*%B%*%t(B)%*%y[n:1]) - TCnew[, c(2, 1)]


B%*%t(B)%*%y[n:1]
gc()
TCnew
plot(x, type="l")
lines(TC[,1], col="2")

Sd%*%t(Sd) + B%*%t(B)


