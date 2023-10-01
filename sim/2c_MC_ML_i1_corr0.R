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

# Wd, etc
setwd("/Users/tobias/Dokumente/Projekte/filtering unknown persistence/R/code")
source("./help functions/KF.R")

# Settings
n <- c(100, 200, 300)
d <- c(0.75, 1, 1.25, 1.75)
r <- c(1, 5, 10)
R <- 1000

# Strong persistence
ar <- c(-1.6, 0.8)
# Weaker persistemce
#ar <- c(-0.9, 0.4)
sig1.grid <- c(-0.5, -0.78, -0.86)
sig2.grid <- c(-0.5, 1.37, 2.15)
sig12.grid <- c(-1.10, -0.96, -0.85)


ar.1.grid <- -c(-1.8, -1.7, -1.6, -1.5, -1.4)
ar2.grid  <- -c(0.6, 0.7, 0.8, 0.9, 1)
gr.start <- expand.grid(sig1.grid, sig12.grid, sig2.grid, ar.1.grid, ar2.grid) %>% as.matrix()
gr.start <- gr.start[sapply(1:nrow(gr.start), function(i) toComp(c(gr.start[i, 4:5]))$stable) ,]

setups <- expand.grid(n, d, r)
colnames(setups) <- c("n", "d", "r")
setups <- setups[order(setups[, "n"]), ]

for ( i in 1:NROW(setups)){
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
    
    optfn <- function(n, y, x){
        # Estimate
        tryCatch({
            grid.st <- cbind(gr.start, sapply(1:nrow(gr.start), function(i) UC_opt_KF_i1_ML(gr.start[i, ], corr=TRUE,
                                                                                            y=y, START = 1, ll=TRUE,
                                                                                            nulim = c(1/100, 25),
                                                                                            deterministics = FALSE)))
            
            par0 <- grid.st[which.min(grid.st[,6]),-6]
            
            est <- optim(par=par0, START=1,
                         fn = UC_opt_KF_i1_ML, ll = TRUE, corr = TRUE,
                         y=y, method = "BFGS", nulim = c(1/100, 25))
            
            
            par <- est$par
            nu  <- mlogvech2mat(est$par[1:3])
            ar  <- par[-(1:3)]
            KS <- fUC_smooth(y, 1, nu, ar = -ar, corr=TRUE)
            KF <- fUC_comp(y, 1, nu, ar=-ar, corr=TRUE)
            
            # Calculate Rsq, etc
            SSR <- mean((x - KS$x)^2)
            SST <- mean((x - mean(x))^2)
            Rsq <- summary(lm(x ~ KS$x))$r.squared
            
            # Return results
            results <- c(nu[1,1], nu[2,2], nu[1,2], ar, SSR, Rsq)
            names(results) <- c("sigma_eta", "sigma_eps", "corr", "ar_1", "ar_2","SSR", "Rsq")
            
            
            return(results)
        }, error = function(e) return(rep(NA, 7))
        )
        
    }
    ### fUC part
    cl <- makeCluster(6)
    clusterExport(cl, c("optfn", "UC_opt_KF_i1_ML", "y", "fUC_comp", "ma_inf",
                        "embed0", "fUC_smooth", "x", "frac_diff", "gr.start"))
    clusterEvalQ(cl, library(fUCpack))
    clusterEvalQ(cl, library(KFAS))
    
    RESULTS <- parSapply(cl, 1:R, function(j) optfn(n, y[, j], x[, j]))
    stopCluster(cl)
    
    
    
    rownames(RESULTS) <- c("sigma_eta", "sigma_eps", "corr", "ar_1", "ar_2", "SSR", "Rsq")
    save(RESULTS, file = paste("./MC/Integer/ML/Sim_R", R, "_n", n, "_d", d, "_nu", nu, "_ar_high_corr.RData", sep=""))
}
