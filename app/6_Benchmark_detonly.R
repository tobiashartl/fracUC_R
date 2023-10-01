# new benchmark
gc()
rm(list = ls())
library(fUCpack)
library(dplyr)
library(parallel)
library(CFFpack)


# Wd, etc
setwd("/Users/tobias/Dokumente/Projekte/filtering unknown persistence/R/")
source("./help functions/fUC_arma_approx.R")
source("./help functions/UC_dettrend_only.R")

# load data
data<- read.csv(file = "./Applications/temptest/noaa_SST.csv", skip = 4)[, c(1, 2)] %>%
    as.data.frame() 


y <- data[,2]
R <- 100
corr <- TRUE
eta  <- 0.025
plot(y)
START = 2
diffuse = T
seas = F
det = "trend"


for(p in c(4)){
    
    ### get starting values
    set.seed(42)
    
    # Q
   
    Sigma <- cbind(runif(R, -10, 0))
    Sigvec <- t((Sigma))
    
    
    # apply(Sigvec, 2, function(x) mlogvech2mat(x))
    
    if(p > 0){
        # Starting values for Phi
        as <- matrix(NA,R,p)
        as[,1] <- runif(R,0.5, 1)
        if (p>1)
        {
            as[,2:p] <- runif(R*(p-1),-0.7, 0.7)
            as <- t(apply(as,1,arfima::PacfToAR))
        }
        A <- -as
    }else{
        A <- NULL
    }
    if(p>0){
        if(!is.null(eta)){
            ev <- apply(-A, 1, function(x) max(abs(toComp(x)$eigv)))
            while(any(ev > 1-2*eta)){
                
                A[ev > 1-2*eta, 1] <-  runif(sum(ev > 1-2*eta),0.7, 1)
                if(p>1){
                    A[ev > 1-2*eta, 2:p] <- runif(sum(ev > 1-2*eta)*(p-1),-0.7, 0.7)
                    A[ev > 1-2*eta, ] <- -t(apply(A[ev > 1-2*eta, , drop=F],1,arfima::PacfToAR))
                }
                ev <- apply(-A, 1, function(x) max(abs(toComp(x)$eigv)))
            }
            
        } 
    }
    
    
    #Seav <- runif(R, -15, 0)
    #Irrv <- runif(R, 0, 7)
    d <- runif(R, 1/2, 2)
    
    START.val <- cbind(d, t(Sigvec), A)#, Irrv)
    RESULTS <- matrix(NA, nrow = R, ncol = ncol(START.val)+1)
    optfn <- function(par0, y){ tryCatch({
        est <- optim(par=par0, method = "BFGS", 
                     fn = UC_dettrend_opt_ML, 
                     y=y,  START = START,
                     nulim = c(0, Inf), quiet=T,
                     pq=c(p, 0),
                     diffuse = diffuse, return.det=F) # custom settings for ML

        
        cat("ll = ", est$value, ", par = ", est$par, "\n")  
        return(c(est$value, est$par))
    }, 
    error = function(e) return(rep(NA, 1 + ncol(START.val))))
    }
    
    #optfn(START.val[1,], y)
    
    system.time({
        cl <- makeCluster(6)
        clusterExport(cl, ls())
        clusterEvalQ(cl, library(ucminf))
        clusterEvalQ(cl, library(fUCpack))
        clusterEvalQ(cl, library(CFFpack))
        
        RESULTS <- parSapply(cl, 1:R, function(j) optfn(START.val[j,], y))
        stopCluster(cl)
    }
    )
    
    
    RESULTS <- t(RESULTS)
    save(RESULTS, file = paste("./Applications/temptest/Benchmark_ML_detonly_p", p, ".RData", sep=""))
    RESULTS.win <- RESULTS[which.min(RESULTS[,1]), ]
    
    
    cat("Results: det = ", det, ", p = ", p, "\n", "LL = ", RESULTS.win[1], ", theta = ", 
        RESULTS.win[-1], "\n")
}



for(p in c(2:12)){
    cat("p = ", p, "\n")
    load(file = paste("./Applications/temptest/Benchmark_ML_i1_cor_p", p, "_trend.RData", sep=""))
    RESULTS <- RESULTS[order(RESULTS[,1]),] %>%
        as.data.frame %>%
        subset(!is.na(V1)) 
    RESULTS$correlation <- apply(RESULTS[,2:4], 1, function(x) cov2cor(mlogvech2mat(x))[2,1])
    
    RESULTS.ar <- RESULTS[, 5:(4+p), drop=F]
    stable <- apply(RESULTS.ar, 1, function(x) all(abs(toComp(-x)$eigv) < 0.99))
    RESULTS <- RESULTS[stable, ]
    
    
    
    ndet = 1
    
    k <- ncol(RESULTS) - 2 + ndet
    
    RESULTS$BIC         <- k * log(length(y)) + 2*RESULTS[,1]
    RESULTS$AIC         <- 2 * k              + 2*RESULTS[,1]
    
    RESULTS[1:1,] %>% print
    
}

p <- 4
load(file = paste("./Applications/temptest/Benchmark_ML_i1_cor_p", p, "_trend.RData", sep=""))
theta <- RESULTS[which.min(RESULTS[,1]),-1]


# confindence bands
hessian <- optimHess(theta, 
                     fn = UC_opt_ML, y=y,
                     nulim = c(0, Inf), 
                     corr=TRUE, deterministics = T, START=2, diffuse = T, 
                     return.det = F, nu.opt=F, 
                     eta=NULL, Q.trans = "mlv", flip=F, neg=F, pq=c(p, 0), penalty.corr=F, 
                     seas=F, period=12, irregular=F)

J <- numDeriv::jacobian(mlogvech2mat, theta[1:3])
cov0 <- solve(hessian)
Cov <- bdiag(J[-2,], diag(p))%*%cov0%*%t(bdiag(J[-2, ], diag(p)))
theta.trans <- c(mlogvech2mat(theta[1:3])[c(TRUE, TRUE, FALSE, TRUE)], theta[-(1:3)])
se.trans <- sqrt(abs(diag(Cov)))
EST <- cbind(c(NA, theta.trans), c(NA, se.trans)) 
colnames(EST) <- c("par", "se") 
rownames(EST) <- c("d", "Q11", "Q21", "Q22", paste("ar_", 1:p, sep=""))
EST
EST[,1] / EST[,2]


Q <- mlogvech2mat(theta[1:3])


det <- UC_opt_ML(theta, y, nulim = c(0, Inf), 
                 corr=TRUE, deterministics = T, START=2, diffuse = T, 
                 return.det = TRUE, nu.opt=F, 
                 eta=NULL, Q.trans = "mlv", flip=F, neg=F, pq=c(p, 0), penalty.corr=F, 
                 seas=F, period=12, irregular=F)
n <- length(y)
Z <- matrix(frac_diff(rep(1, n), -1), n)
#ZZ <- embed0(rep(c(1, rep(0, 11)), length.out = length(y)), 12)
dettrend <- cbind(Z) %*% det
plot(y, type="l")
lines(dettrend + y[1], col="2")


TC <- fUC_KS_ARMA_approx(c(1, theta), y - Z %*% det, nulim = c(0, Inf), 
                         corr=TRUE, deterministics = "none", START=1, diffuse = T, 
                         return.det = F, nu.opt=F, d.int = c(0, 2.5), 
                         eta=NULL, Q.trans = "mlv", flip=F, neg=F, pq=c(p, 0), penalty.corr=F, 
                         seas=F, period=12, irregular=F)
plot(TC$x, type="l")
plot(TC$x + dettrend)
lines(y, col="2")
plot(ts(TC$c, start = c(1850, 1), frequency = 12))
abline(h=0)

ll <- fUC_opt_ML(theta = c(1, theta), 
                 y=y, nulim =c(0, Inf), quiet = TRUE, corr=TRUE, deterministics = "frac", 
                 START=2, diffuse =TRUE, return.det=F, nu.opt = F, d.int = c(0, 2.5), 
                 eta = NULL, Q.trans = "mlv", flip =F, neg=F)



k <- length(theta)+1
(BIC          <- k * log(length(y)) + 2*ll)
(AIC          <- 2 * k              + 2*ll)

SSR <- (fUC_comp((y - dettrend), 1, Q, theta[-(1:3)], corr=TRUE)$v[-1])^2 %>% sum(.)
add <- c(Q[2,2]/Q[1,1], Q[2,1]/Q[1,1], cov2cor(Q)[2,1], ll, SSR, AIC, BIC)

EST <- rbind(EST, cbind(add, NA))
rownames(EST) <- c("d", "Q11", "Q21", "Q22", paste("ar_", 1:p, sep=""),
                   "nu_1", "nu_2", "rho", "ll", "CSS", "AIC", "BIC")
colnames(EST)<- c("par", "se")
saveRDS(EST, file = "./Applications/temptest/Benchmark_i1_ML_Results.RDS")
data.plot <- data.frame(
    y=y,
    trend = TC$x + dettrend, 
    x = TC$x, 
    det = dettrend, 
    cycle = TC$c,
    time = seq(from = as.Date("1850-01-01"), to = as.Date("2023-07-01"),  by = "month"))
saveRDS(data.plot, file = "./Applications/temptest/Benchmark_i1_ML_TC.RDS")

