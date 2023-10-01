gc()
rm(list = ls())
library(fUCpack)
library(dplyr)

# Wd, etc
setwd("/Users/tobias/Dokumente/Projekte/filtering unknown persistence/R/code")
source("./help functions/KF.R")
# Settings
n <- c(100, 200, 300)
d <- c(0.75, 1, 1.25, 1.75)
r <- c(1, 5, 10)
R <- 1000


setups <- expand.grid(n, d, r)
colnames(setups) <- c("n", "d", "r")
setups <- setups[order(setups[, "n"]), ]


# Investigate the precision of x
# CSS
# Load CSS esimates
R <- 1000
ar <- c(-1.6, 0.8)
results.list <- list()
X.res <- list()
C.res <- list()
X.true <- list()
C.true <- list()
for(i in 1:nrow(setups)){
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    nu <- as.numeric(setup[3])
    d <- as.numeric(setup[2])
    
    set.seed(42)
    
    rho <- -0.8
    Q <- matrix(c(1, rho*sqrt(1)*sqrt(nu), rho*sqrt(1)*sqrt(nu), nu), 2, 2)
    
    # Generate the data
    U <- sapply(1:R, function(x) mvtnorm::rmvnorm(n, mean=c(0, 0), sigma = Q))
    x <- frac_diff_multi(U[1:n,], d=-d)
    u <- U[-(1:n),]
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- x + c
    
    try({
        load(paste("./MC/Sim_R1000_n", n, "_d", d, "_nu", nu,
                   "_ar_high_corr.RData", sep=""))
        results.list[[i]] <- RESULTS
        
        filterfun <- function(theta, y){
            nu  <- matrix(c(1, theta[3], theta[3], theta[2]), 2, 2)
            KS <- fUC_smooth(y, theta[1], nu, theta[-(1:3)], corr=TRUE)
            return(c(KS$x, KS$c))
        }
        
        XC <- sapply(1:R, function(j) filterfun(RESULTS[1:5,j], y[,j]))
        X.mat <- XC[1:n, ]
        C.mat <- XC[-(1:n), ]
    })
    X.res[[i]] <- X.mat
    C.res[[i]] <- C.mat
    cat("Iteration ", i, "\n")
}

save(X.res, C.res, results.list, X.true, C.true, 
     file = "./MC/results_cor.RData")


# ML
# Load esimates
R <- 1000
ar <- c(-1.6, 0.8)
results.list <- list()
X.res <- list()
C.res <- list()
X.true <- list()
C.true <- list()
for(i in 1:nrow(setups)){
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    nu <- as.numeric(setup[3])
    d <- as.numeric(setup[2])
    
    set.seed(42)
    
    rho <- -0.8
    Q <- matrix(c(1, rho*sqrt(1)*sqrt(nu), rho*sqrt(1)*sqrt(nu), nu), 2, 2)
    
    # Generate the data
    U <- sapply(1:R, function(x) mvtnorm::rmvnorm(n, mean=c(0, 0), sigma = Q))
    x <- frac_diff_multi(U[1:n,], d=-d)
    u <- U[-(1:n),]
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- x + c
    try({
        load(paste("./MC/ML/Sim_R1000_n", n, "_d", d, "_nu", nu,
                   "_ar_high_corr.RData", sep=""))
        results.list[[i]] <- RESULTS[, colSums(is.na(RESULTS))==0 & RESULTS[1, ] < 1e+06]
        
        
        filterfun <- function(theta, y){
            nu  <- mlogvech2mat(theta[2:4])
            KS <- fUC_smooth(y, theta[1], nu, theta[-(1:4)], corr=TRUE)
            return(c(KS$x, KS$c))
        }
        
        XC <- sapply(1:R, function(j) filterfun(RESULTS[2:7,j], y[,j]))
        X.mat <- XC[1:n, ]
        C.mat <- XC[-(1:n), ]
    })
    X.res[[i]] <- X.mat
    C.res[[i]] <- C.mat
    cat("Iteration ", i, "\n")
}

save(X.res, C.res, results.list, X.true, C.true, 
     file = "./MC/results_ML_corr.RData")


# I(1) CSS
# Load esimates
R <- 1000
ar <- c(-1.6, 0.8)
results.list <- list()
X.res <- list()
C.res <- list()
X.true <- list()
C.true <- list()
for(i in 1:nrow(setups)){
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    nu <- as.numeric(setup[3])
    d <- as.numeric(setup[2])
    
    #if(file.exists(file = paste("./MC/MC_3/Integer/Sim_R", R, "_n", n, "_d", d, "_nu", nu, "_ar_high.RData", sep=""))) next
    set.seed(42)
    
    rho <- -0.8
    Q <- matrix(c(1, rho*sqrt(1)*sqrt(nu), rho*sqrt(1)*sqrt(nu), nu), 2, 2)
    
    # Generate the data
    U <- sapply(1:R, function(x) mvtnorm::rmvnorm(n, mean=c(0, 0), sigma = Q))
    x <- frac_diff_multi(U[1:n,], d=-d)
    u <- U[-(1:n),]
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- x + c
    try({
        load(paste("./MC/INTEGER/Sim_R1000_n", n, "_d", d, "_nu", nu,
                   "_ar_high_corr.RData", sep=""))
        results.list[[i]] <- RESULTS
        
        
        filterfun <- function(theta, y){
            Q <- matrix(c(1, theta[2], theta[2], theta[1]), 2, 2)
            KS <- fUC_smooth(y, 1, Q, -theta[-(1:2)], corr=TRUE)
            return(c(KS$x, KS$c))
        }
        
        XC <- sapply(1:R, function(j) filterfun(RESULTS[1:4,j], y[,j]))
        X.mat <- XC[1:n, ]
        C.mat <- XC[-(1:n), ]
    })
    X.res[[i]] <- X.mat
    C.res[[i]] <- C.mat
    cat("Iteration ", i, "\n")
}

save(X.res, C.res, results.list, X.true, C.true, 
     file = "./MC/results_i1_corr.RData")


# I(1) ML
# Load esimates
R <- 1000
ar <- c(-1.6, 0.8)
results.list <- list()
X.res <- list()
C.res <- list()
X.true <- list()
C.true <- list()
for(i in 1:nrow(setups)){
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    nu <- as.numeric(setup[3])
    d <- as.numeric(setup[2])
    
    set.seed(42)
    
    rho <- -0.8
    Q <- matrix(c(1, rho*sqrt(1)*sqrt(nu), rho*sqrt(1)*sqrt(nu), nu), 2, 2)
    
    # Generate the data
    U <- sapply(1:R, function(x) mvtnorm::rmvnorm(n, mean=c(0, 0), sigma = Q))
    x <- frac_diff_multi(U[1:n,], d=-d)
    u <- U[-(1:n),]
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- x + c
    try({
        load(paste("./MC/INTEGER/ML/Sim_R1000_n", n, "_d", d, "_nu", nu,
                   "_ar_high_corr.RData", sep=""))
        results.list[[i]] <- RESULTS[, colSums(is.na(RESULTS))==0]
        
        
        filterfun <- function(theta, y){
            nu  <- matrix(c(theta[1], theta[3], theta[3], theta[2]), 2, 2)
            KS <- fUC_smooth(y, 1, nu, -theta[-(1:3)], corr=TRUE)
            
            
            
            return(c(KS$x, KS$c))
        }
        
        XC <- sapply(1:R, function(j) tryCatch(filterfun(RESULTS[1:5,j], y[,j]),
                                               error = function(e) return(rep(NA, 2*n))
        )
        )
        
        X.mat <- XC[1:n, ]
        C.mat <- XC[-(1:n), ]
    })
    
    
    calcRSQ <- function(x.true, x.est){
        Rsq <- mean(sapply(1:ncol(x.true), function(j) tryCatch(summary(lm(x.true[, j] ~ x.est[,j]))$r.squared,
                                                                error = function(e) return(NA))), na.rm = TRUE)
    }
    calcRSQ(x, X.mat)-> rsq
    
    X.res[[i]] <- X.mat
    C.res[[i]] <- C.mat
    cat("Iteration ", i, "\n")
}

save(X.res, C.res, results.list, X.true, C.true, 
     file = "./MC/results_i1_ML_corr.RData")









