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
setwd("~/Dokumente/Projekte/filtering unknown persistence/R/code")


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
d.grid <- seq(0.75, 1.75, 0.25)
r.grid <- (c(1, 5, 10))
ar.1.grid <- c(-1.8, -1.7, -1.6, -1.5, -1.4)
ar2.grid  <- c(0.6, 0.7, 0.8, 0.9, 1)
corr_grid <- c(tcrossprod(-c(0.7, 0.8, 0.9), sqrt(r.grid)))

gr.start <- expand.grid(d.grid, r.grid, corr_grid, ar.1.grid, ar2.grid) %>% as.matrix()
gr.start <- gr.start[sapply(1:nrow(gr.start), function(i) toComp(-c(gr.start[i, 4:5]))$stable) ,]

setups <- expand.grid(n, d, r)
colnames(setups) <- c("n", "d", "r")
setups <- setups[order(setups[, "n"]), ]

for ( i in 1:NROW(setups)){
  setup <- setups[ i, ]
  n <- as.numeric(setup[1])
  nu <- as.numeric(setup[3])
  d <- as.numeric(setup[2])
  cat("Iteration ", i, "\n")
  set.seed(42)
  
  # Generate the data
  Q   <- matrix(c(1, -0.8*sqrt(nu), -0.8*sqrt(nu), nu), 2, 2)
  ETA <- EPS <- matrix(NA, nrow = n, ncol=R)
  
  for(j in 1:R){
    ERR <-   mvtnorm::rmvnorm(n, mean=c(0, 0), sigma = Q)
    ETA[, j] <- ERR[,1]
    EPS[, j] <- ERR[,2]
  }
  
  x <- frac_diff_multi(ETA, -d)
  c <- apply(EPS, 2, function(x) stats::filter(c(rep(0, length(ar)), x), filter = c(-ar), method = "recursive", sides = 1)[-(1:length(ar))])
  y <- yy<-x + c
  
  optfn <- function(n, y, x){
    tryCatch({
      # Estimate
      # check grid
      grid.st <- cbind(gr.start, sapply(1:nrow(gr.start), function(i) fUC_opt(gr.start[i, ], corr=TRUE,
                                                                              y=y, START = 1,
                                                                              nulim = c(1/100, 25),
                                                                              deterministics = FALSE)))
      
      par0 <- grid.st[which.min(grid.st[,6]),-6]
      
      (est <- optim(par=par0, control = list(), corr=TRUE,
                    fn = fUC_opt, quiet = TRUE, START = 1,
                    y=y, method = "BFGS", nulim = c(1/100, 25), d.int = c(0.5, 2.5),
                    deterministics = FALSE))
      
      # check for corner solutions: d->0, d->2, nu->1/100, nu->100
      j=1
      while(est$par[1] < .05 | est$par[1] > 1.95 | exp(est$par[2]) < .05 | exp(est$par[2]) > 24.95){
        if(j > 10) break
        j=j+1
        par0 <- grid.st[order(grid.st[,6]),-6][j , ]
        (est <- optim(par=par0, corr=TRUE,
                      fn = fUC_opt, quiet = TRUE, d.int = c(0.5, 2.5),
                      y=y, method = "BFGS", nulim = c(1/100, 25),
                      deterministics = FALSE, START = 1))
        
      }
      
      
      
      
      par <- est$par
      nu  <- (par[2])
      corr <- par[3]
      ar <- est$par[-(1:3)]
      #KF <- fUC_comp(y, par[1], nu, ar)
      #KF2 <- fUC_comp((1:length(y))^par[1], par[1], nu, ar)
      #mu  <- summary(lm(KF$v ~ -1 + KF2$v))$coefficients[1,1]
      Q <- matrix(c(1, corr, corr, nu), 2, 2)
      
      KS <- fUC_smooth(y, par[1], Q, ar = ar, corr = TRUE)
      Rsq <- summary(lm(x ~ KS$x))$r.squared
      
      
      # Return results
      results <- c(par[1], (par[2]), par[3], ar, Rsq)
      names(results) <- c("d", "nu", "corr", "ar_1", "ar_2", "Rsq")
      return(results)
    }, error = function(e) return(rep(NA, 6))
    )
    
  }
  
  #optfn(n, y[,3], x[,3])
  
  ### fUC part
  cl <- makeCluster(7)
  clusterExport(cl, c("optfn", "fUC_opt", "y", "fUC_comp", "ma_inf",
                      "embed0", "fUC_smooth", "x", "frac_diff", "lm", "n", "gr.start"))
  clusterEvalQ(cl, library(fUCpack))
  RESULTS <- parSapply(cl, 1:R, function(j){ 
    optfn(n, y[, j], x[, j])})
  stopCluster(cl)
  
  dEW_45 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2), alpha = .45))
  dEW_50 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2), alpha = .50))
  dEW_55 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2), alpha = .55))
  dEW_60 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2), alpha = .60))
  dEW_65 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2), alpha = .65))
  dEW_70 <- c(apply(y, 2, EW, type = 1, interval = c(0, 2), alpha = .70))
  
  RESULTS <- rbind(RESULTS, c(dEW_45), c(dEW_50), c(dEW_55), c(dEW_60), c(dEW_65), c(dEW_70))
  rownames(RESULTS) <- c("d", "nu", "corr", "ar_1", "ar_2", "Rsq", 
                         "d_45", "d_50", "d_55", "d_60", "d_65", "d_70")
  save(RESULTS, file = paste("./MC/Sim_R", R, "_n", n, "_d", d, "_nu", nu, "_ar_high_corr.RData", sep=""))
}

