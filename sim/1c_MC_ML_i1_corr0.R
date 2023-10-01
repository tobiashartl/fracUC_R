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

# Wd, etc
setwd("/Users/tobias/Dokumente/Projekte/Filtering unknown persistence/R/code")
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
# define grid
r.grid <- log(c(1, 5, 10))
ar.1.grid <- -c(-1.8, -1.7, -1.6, -1.5, -1.4)
ar2.grid  <- -c(0.6, 0.7, 0.8, 0.9, 1)

gr.start <- expand.grid(r.grid, ar.1.grid, ar2.grid) %>% as.matrix()
gr.start <- gr.start[sapply(1:nrow(gr.start), function(i) toComp(c(gr.start[i, 2:3]))$stable) ,]


setups <- expand.grid(n, d, r)
colnames(setups) <- c("n", "d", "r")
setups <- setups[order(setups[, "n"]), ]

for ( i in 1:NROW(setups)){
    setup <- setups[ i, ]
    n <- as.numeric(setup[1])
    nu <- as.numeric(setup[3])
    d <- as.numeric(setup[2])
    
    #if(file.exists(file = paste("./MC/MC_2/Integer/Sim_R", R, "_n", n, "_d", d, "_nu", nu, "_ar_high_ML.RData", sep=""))) next
    set.seed(42)
    
    # Generate the data
    x <- frac_diff_multi(matrix(rnorm(n*R, mean=0, sd=1), n, R), -d)
    u <- matrix(rnorm(n*R, mean=0, sd=sqrt(nu)), n, R)
    c <- apply(u, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
    y <- x + c
    
    optfn <- function(n, y, x){
        # Estimate
        tryCatch({
            grid.st <- cbind(gr.start, sapply(1:nrow(gr.start), function(i) UC_opt_KF_i1(gr.start[i, ], 
                                                                                         y=y, START = 1, ll=TRUE,
                                                                                         nulim = c(1/100, 25),
                                                                                         deterministics = FALSE)))
            
            par0 <- grid.st[which.min(grid.st[,4]),-4]
            est <- optim(par=par0,
                         fn = UC_opt_KF_i1, ll = TRUE,
                         y=y, method = "BFGS", nulim = c(1/100, 100))
            
            
            par <- est$par
            nu  <- exp(par[1])
            ar  <- par[-1]
            KS <- fUC_smooth(y, 1, nu, ar = ar, corr=FALSE)
            KF <- fUC_comp(y, 1, nu, ar=ar)
            
            # Calculate Rsq, etc
            SSR <- mean((x - KS$x)^2)
            SST <- mean((x - mean(x))^2)
            Rsq <- summary(lm(x ~ KS$x))$r.squared
            
            # Return results
            results <- c(exp(par[1]), ar, SSR, Rsq)
            names(results) <- c("nu", "ar_1", "ar_2","SSR", "Rsq")
            
            
            return(results)
        }, error = function(e) return(rep(NA, 5))
        )
        
    }
    #optfn(n, y[,1], x[,1])
    ### fUC part
    cl <- makeCluster(6)
    clusterExport(cl, c("optfn", "UC_opt_KF_i1", "y", "fUC_comp", "ma_inf",
                        "embed0", "fUC_smooth", "x", "frac_diff", "gr.start"))
    clusterEvalQ(cl, library(fUCpack))
    clusterEvalQ(cl, library(KFAS))
    
    RESULTS <- parSapply(cl, 1:R, function(j) optfn(n, y[, j], x[, j]))
    stopCluster(cl)
    
    
    
    rownames(RESULTS) <- c("nu", "ar_1", "ar_2", "SSR", "Rsq")
    save(RESULTS, file = paste("./MC/Integer/ML/Sim_R", R, "_n", n, "_d", d, "_nu", nu, "_ar_high_ML.RData", sep=""))
}
