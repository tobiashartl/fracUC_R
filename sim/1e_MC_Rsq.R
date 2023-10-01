gc()
rm(list = ls())
library(fUCpack)
library(dplyr)

# Wd, etc
setwd("/Users/tobias/Dokumente/Projekte/Filtering unknown persistence/R/code")
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
    true.par          <- setups[i, ]
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    nu <- as.numeric(setup[3])
    d <- as.numeric(setup[2])
    set.seed(42)
    x <- frac_diff_multi(matrix(rnorm(n*R, mean=0, sd=1), n, R), -d)
    u <- matrix(rnorm(n*R, mean=0, sd=sqrt(nu)), n, R)
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- x + c
    X.mat <- x
    C.mat <- c
    X.true[[i]] <- X.mat
    C.true[[i]] <- C.mat
    try({
        load(paste("./MC/Sim_R1000_n", true.par[1], "_d", true.par[2], "_nu", true.par[3],
                   "_ar_high_corr0.RData", sep=""))
        results.list[[i]] <- RESULTS
        
        
        filterfun <- function(theta, y){
            KS <- fUC_smooth(y, theta[1], theta[2], theta[-(1:2)], corr=FALSE)
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
     file = "./MC/results.RData")


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
    true.par          <- setups[i, ]
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    nu <- as.numeric(setup[3])
    d <- as.numeric(setup[2])
    set.seed(42)
    x <- frac_diff_multi(matrix(rnorm(n*R, mean=0, sd=1), n, R), -d)
    u <- matrix(rnorm(n*R, mean=0, sd=sqrt(nu)), n, R)
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- x + c
    X.mat <- x
    C.mat <- c
    X.true[[i]] <- X.mat
    C.true[[i]] <- C.mat
    try({
        load(paste("./MC/ML/Sim_R1000_n", true.par[1], "_d", true.par[2], "_nu", true.par[3],
                   "_ar_high_corr0.RData", sep=""))
        results.list[[i]] <- RESULTS
        
        
        filterfun <- function(theta, y){
            KS <- fUC_smooth(y, theta[1], theta[2], theta[-(1:2)], corr=FALSE)
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
     file = "./MC/results_ML.RData")


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
    true.par          <- setups[i, ]
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    nu <- as.numeric(setup[3])
    d <- as.numeric(setup[2])
    set.seed(42)
    x <- frac_diff_multi(matrix(rnorm(n*R, mean=0, sd=1), n, R), -d)
    u <- matrix(rnorm(n*R, mean=0, sd=sqrt(nu)), n, R)
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- x + c
    X.mat <- x
    C.mat <- c
    X.true[[i]] <- X.mat
    C.true[[i]] <- C.mat
    try({
        load(paste("./MC/INTEGER/Sim_R1000_n", true.par[1], "_d", true.par[2], "_nu", true.par[3],
                   "_ar_high.RData", sep=""))
        results.list[[i]] <- RESULTS
        
        
        filterfun <- function(theta, y){
            KS <- fUC_smooth(y, 1, theta[1], -theta[-(1:1)], corr=FALSE)
            return(c(KS$x, KS$c))
        }
        
        XC <- sapply(1:R, function(j) filterfun(RESULTS[1:3,j], y[,j]))
        X.mat <- XC[1:n, ]
        C.mat <- XC[-(1:n), ]
    })
    X.res[[i]] <- X.mat
    C.res[[i]] <- C.mat
    cat("Iteration ", i, "\n")
}

save(X.res, C.res, results.list, X.true, C.true, 
     file = "./MC/results_i1.RData")


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
    true.par          <- setups[i, ]
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    nu <- as.numeric(setup[3])
    d <- as.numeric(setup[2])
    set.seed(42)
    x <- frac_diff_multi(matrix(rnorm(n*R, mean=0, sd=1), n, R), -d)
    u <- matrix(rnorm(n*R, mean=0, sd=sqrt(nu)), n, R)
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- x + c
    X.mat <- x
    C.mat <- c
    X.true[[i]] <- X.mat
    C.true[[i]] <- C.mat
    try({
        load(paste("./MC/INTEGER/ML/Sim_R1000_n", true.par[1], "_d", true.par[2], "_nu", true.par[3],
                   "_ar_high_ML.RData", sep=""))
        results.list[[i]] <- RESULTS
        
        
        filterfun <- function(theta, y){
            KS <- fUC_smooth(y, 1, theta[1], -theta[-(1:1)], corr=FALSE)
            return(c(KS$x, KS$c))
        }
        
        XC <- sapply(1:R, function(j) filterfun(RESULTS[1:3,j], y[,j]))
        X.mat <- XC[1:n, ]
        C.mat <- XC[-(1:n), ]
    })
    X.res[[i]] <- X.mat
    C.res[[i]] <- C.mat
    cat("Iteration ", i, "\n")
}

save(X.res, C.res, results.list, X.true, C.true, 
     file = "./MC/results_i1_ML.RData")








# check goodness of fit
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    nu <- as.numeric(setup[3])
    d <- as.numeric(setup[2])
    
    X.est <- X.res[[i]]
    C.est <- C.res[[i]]
    X.tr  <- X.true[[i]]
    C.tr  <- C.true[[i]]
    
    # R^2_x 
    Rsq_x <- sapply(1:R, function(x) summary(lm(X.tr[, x] ~ X.est[, x]))$r.squared)
    Rsq_c <- sapply(1:R, function(x) summary(lm(C.tr[, x] ~ C.est[, x]))$r.squared)
    #Rsq_eta <- sapply(1:R, function(x) summary(lm(frac_diff(X.tr[, x], d) ~ frac_diff(X.est[, x], ))$r.squared)
}









