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

Qgen <- function(par){
    Q <- matrix(NA, 2, 2)
    Q[1,1] <- par[1]
    Q[2,2] <- par[2]
    Q[1,2] <- Q[2,1] <- par[3]
    return(cov2cor(Q)[2,1])
}




# CSS
results.list <- list()
statistics <- statisticsb <- matrix(NA, nrow(setups), 12)
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    try({
        load(paste("./MC/Sim_R1000_n", true.par[1], "_d", true.par[2], "_nu", true.par[3],
                   "_ar_high_corr.RData", sep=""))
        results.list[[i]] <- RESULTS[, colSums(is.na(RESULTS))==0 & RESULTS[1, ] < 1e+06]
        rho_hat <- apply(rbind(1, results.list[[i]][c("nu", "corr"), ]), 2, Qgen)
        
        
        MSE_d <- rowMeans((results.list[[i]][c("d", "d_45", "d_50", "d_55", "d_60", "d_65", "d_70"), ] - true.par$d)^2) %>%
            sqrt(.)
        #MSE_eta <- mean((results.list[[i]]["sigma_eta", ] - 1)^2)%>%
        #    sqrt(.)
        MSE_eps <- mean((results.list[[i]]["nu", ] - true.par$r)^2)%>%
            sqrt(.)
        MSE_rho <- mean((rho_hat + 0.8)^2)%>%
            sqrt(.)
        MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8))^2)%>%
            sqrt(.)
        
        
        Rsq   <- mean(results.list[[i]][c("Rsq"), ])
        
        statistics[i, ] <- c(MSE_d, MSE_eps, MSE_rho, MSE_ar, Rsq)
        
        
        MSE_d <- rowMeans((results.list[[i]][c("d", "d_45", "d_50", "d_55", "d_60", "d_65", "d_70"), ] - true.par$d)) 
        MSE_eps <- mean((results.list[[i]]["nu", ] - true.par$r))
        MSE_rho <- mean((rho_hat + 0.8))
        MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8)))
        
        statisticsb[i, ] <- c(MSE_d, MSE_eps, MSE_rho,  MSE_ar, Rsq)
    })
    
}


rownames(statistics) <- rownames(statisticsb) <- apply(setups[1:nrow(setups),], 1, function(x) paste(x, collapse="_"))
colnames(statistics) <- colnames(statisticsb) <- c("d", "d_45", "d_50", "d_55", "d_60", "d_65", "d_70", "nu", "rho", "ar1", "ar2", "Rsq")




# QML
results.list <- list()
statistics2 <- statistics2b <- matrix(NA, nrow(setups), 8)
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    try({
        load(paste("./MC/ML/Sim_R1000_n", true.par[1], "_d", true.par[2], "_nu", true.par[3],
                   "_ar_high_corr.RData", sep=""))
        results.list[[i]] <- RESULTS[, colSums(is.na(RESULTS))==0 & RESULTS[1, ] < 1e+06]
        rho_hat <- apply(results.list[[i]][c("Q11", "Q21", "Q22"), ], 2, function(x) cov2cor(mlogvech2mat(x))[2,1])
        sigma_eta <- apply(results.list[[i]][c("Q11", "Q21", "Q22"), ], 2, function(x) (mlogvech2mat(x))[1,1])
        sigma_eps <- apply(results.list[[i]][c("Q11", "Q21", "Q22"), ], 2, function(x) (mlogvech2mat(x))[2,2])
        
        
        MSE_d <- mean((results.list[[i]][c("d"), ] - true.par$d)^2) %>%
            sqrt(.)
        MSE_eta <- mean((sigma_eta - 1)^2)%>%
            sqrt(.)
        MSE_eps <- mean((sigma_eps - true.par$r)^2)%>%
            sqrt(.)
        MSE_nu <- mean((sigma_eps/sigma_eta - true.par$r)^2)%>%
            sqrt(.)
        MSE_rho <- mean((rho_hat + 0.8)^2)%>%
            sqrt(.)
        MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8))^2)%>%
            sqrt(.)
        
        
        Rsq   <- mean(results.list[[i]][c("Rsq"), ])
        
        statistics2[i, ] <- c(MSE_d, MSE_eta, MSE_eps, MSE_nu, MSE_rho, MSE_ar, Rsq)
        
        
        MSE_d <- mean((results.list[[i]][c("d"), ] - true.par$d)) 
        MSE_eta <- mean((sigma_eta - 1))
        MSE_eps <- mean((sigma_eps - true.par$r))
        MSE_nu <- mean((sigma_eps/sigma_eta - true.par$r))
        
        MSE_rho <- mean((rho_hat + 0.8))
        MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8)))
        
        statistics2b[i, ] <- c(MSE_d, MSE_eta, MSE_eps, MSE_nu, MSE_rho,  MSE_ar, Rsq)
    })
    
}

rownames(statistics2) <- rownames(statistics2b) <- apply(setups[1:nrow(setups),], 1, function(x) paste(x, collapse="_"))
colnames(statistics2) <- colnames(statistics2b) <- c("d_ML", "eta_ML", "eps_ML", "nu_ML", "rho_ML", "ar1_ML", "ar2_ML", "Rsq_ML")




# Load I(1) CSS estimates
results.list <- list()
statistics3 <- statistics3b <- matrix(NA, nrow(setups), 5)
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    try({
        load(paste("./MC/Integer/Sim_R1000_n", true.par[1], "_d", true.par[2], "_nu", true.par[3],
                   "_ar_high_corr.RData", sep=""))
        
        
        results.list[[i]] <- RESULTS[, RESULTS[1, ] < 1e+06]
        
        SSR   <- mean(results.list[[i]][c("SSR"), ])
        Rsq   <- mean(results.list[[i]][c("Rsq"), ])
        rho_hat <- apply(rbind(1, results.list[[i]][c("sigma_eps", "corr"), ]), 2, Qgen)
        
        MSE_eps <- mean((results.list[[i]]["sigma_eps", ] - true.par$r)^2)%>%
            sqrt(.)
        MSE_rho <- mean((rho_hat + 0.8)^2)%>%
            sqrt(.)
        MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(1.6, -0.8))^2)%>%
            sqrt(.)
        
        
        statistics3[i, ] <- c(MSE_eps, MSE_rho, MSE_ar, Rsq)
        
        MSE_eps <- mean((results.list[[i]]["sigma_eps", ] - true.par$r))
        MSE_rho <- mean((rho_hat + 0.8))
        MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(1.6, -0.8)))
        
        statistics3b[i, ] <- c(MSE_eps, MSE_rho, MSE_ar, Rsq)
        
    })
    
}
rownames(statistics3) <- rownames(statistics3b) <- apply(setups[1:nrow(setups),], 1, function(x) paste(x, collapse="_"))
colnames(statistics3) <- colnames(statistics3b) <- c("nu_i1", "rho_i1", "ar1_i1", "ar2_i1", "Rsq_i1")


#Load I(1) ML estimates
results.list <- list()
statistics4 <- statistics4b <- matrix(NA, nrow(setups), 7)
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    try({
        load(paste("./MC/Integer/ML/Sim_R1000_n", true.par[1], "_d", true.par[2], "_nu", true.par[3],
                   "_ar_high_corr.RData", sep=""))
        results.list[[i]] <- RESULTS[, colSums(is.na(RESULTS))==0 & RESULTS[1, ] < 1e+06]
        #MSE_d <- mean((results.list[[i]][c("d"), ] - true.par$d)^2)
        SSR   <- mean(results.list[[i]][c("SSR"), ])
        Rsq   <- mean(results.list[[i]][c("Rsq"), ])
        rho_hat <- apply(results.list[[i]][c("sigma_eta", "sigma_eps", "corr"), ], 2, Qgen)
        
        MSE_eta <- mean((results.list[[i]]["sigma_eta", ] - 1)^2)%>%
            sqrt(.)
        MSE_nu <- mean((results.list[[i]]["sigma_eps", ]/results.list[[i]]["sigma_eta", ] - true.par$r)^2)%>%
            sqrt(.)
        MSE_eps <- mean((results.list[[i]]["sigma_eps", ] - true.par$r)^2) %>% sqrt(.)
        
        MSE_rho <- mean((rho_hat + 0.8)^2)%>%
            sqrt(.)
        MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(1.6, -0.8))^2)%>%
            sqrt(.)
        
        
        
        statistics4[i, ] <- c(MSE_eta, MSE_eps, MSE_nu, MSE_rho, MSE_ar, Rsq)
        
        MSE_eta <- mean((results.list[[i]]["sigma_eta", ] - 1))
        MSE_eps <- mean((results.list[[i]]["sigma_eps", ]/results.list[[i]]["sigma_eta", ] - true.par$r))
        MSE_rho <- mean((rho_hat + 0.8))
        MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(1.6, -0.8)))
        MSE_eps <- mean((results.list[[i]]["sigma_eps", ] - true.par$r)) 
        MSE_nu <- mean((results.list[[i]]["sigma_eps", ]/results.list[[i]]["sigma_eta", ] - true.par$r))
        
        statistics4b[i, ] <- c(MSE_eta, MSE_eps, MSE_nu, MSE_rho, MSE_ar, Rsq)
        
    })
    
}
rownames(statistics4) <- rownames(statistics4b) <- apply(setups[1:nrow(setups),], 1, function(x) paste(x, collapse="_"))
colnames(statistics4) <- colnames(statistics4b) <- c("eta_i1_ML", "eps_i1_ML", "nu_i1_ML","rho_i1_ML", "ar1_i1_ML", "ar2_i1_ML", "Rsq_i1_ML")




statistics <- statistics %>% as.data.frame%>%tibble::rownames_to_column()
statistics2 <- statistics2 %>% as.data.frame %>% tibble::rownames_to_column()
statistics3 <- statistics3 %>% as.data.frame %>% tibble::rownames_to_column()
statistics4 <- statistics4 %>% as.data.frame %>% tibble::rownames_to_column()


statistics$n <- setups$n
statistics$d0 <- setups$d
statistics$r   <- setups$r
tab <- left_join(statistics, statistics2) %>%
    left_join(statistics3)%>%
    left_join(statistics4)%>%
    dplyr::select("n",  "r", "d0",
                  "d", "d_ML", "d_50", "d_55", "d_60", "d_65","d_70") %>%
    filter(n %in% c(100, 200, 300), # in the end: 100, 200, 300
           r %in% c(1, 5, 10))

tab1 <- tab[with(tab, order(n, r)),]
tab1$n[c(1:36)%in%c(2:12, 14:24, 26:36)] <- NA
tab1$r[!c(1:36)%in%seq(1, 36, by = 4)] <- NA




### Table 2: Bias d
statisticsb <- statisticsb %>% as.data.frame%>%tibble::rownames_to_column()
statistics2b <- statistics2b %>% as.data.frame%>%tibble::rownames_to_column()
statistics3b <- statistics3b %>% as.data.frame %>% tibble::rownames_to_column()
statistics4b <- statistics4b %>% as.data.frame %>% tibble::rownames_to_column()


statisticsb$n <- setups$n
statisticsb$d0 <- setups$d
statisticsb$r   <- setups$r

tab2 <- left_join(statisticsb, statistics3b) %>%
    left_join(statistics2b)%>%
    left_join(statistics4b)%>%
    select("n",  "r", "d0",
           "d", "d_ML", "d_50", "d_55","d_60", "d_65","d_70") %>%
    filter(n %in% c(100, 200, 300), # in the end: 100, 200, 300
           r %in% c(1, 5, 10), 
           d0 %in% c(0.75, 1, 1.25, 1.75))

tab2 <- tab2[with(tab2, order(n, r)),]
tab2$n[c(1:36)%in%c(2:12, 14:24, 26:36)] <- NA
tab2$r[!c(1:36)%in%seq(1, 36, by = 4)] <- NA

tabjoint <- cbind(tab1, tab2[, -(1:3)])
table <- tabjoint%>% 
    xtable::xtable(digits = c(0, 0, 0, 2, rep(3, 14)), include.rownames=FALSE)
print(table, include.rownames = FALSE)





tab3 <- left_join(statistics, statistics3) %>%
    left_join(statistics2)%>%
    left_join(statistics4)%>%
    select("n",  "r", "d0",
           "nu", "nu_i1",
           "rho", "rho_i1", 
           "ar1", "ar1_i1",
           "ar2", "ar2_i1") %>%
    filter(n %in% c(100, 200, 300), # in the end: 100, 200, 300
           r %in% c(1, 5, 10), 
           d0 %in% c(0.75, 1, 1.25, 1.75))


tab3 <- tab3[with(tab3, order(n, r)),]
tab3$n[c(1:36)%in%c(2:12, 14:24, 26:36)] <- NA
tab3$r[!c(1:36)%in%seq(1, 36, by = 4)] <- NA


tab3b <- left_join(statisticsb, statistics3b) %>%
    left_join(statistics2b)%>%
    left_join(statistics4b)%>%
    select("n",  "r", "d0",
           "nu", "nu_i1",
           "rho", "rho_i1", 
           "ar1", "ar1_i1",
           "ar2", "ar2_i1") %>%
    filter(n %in% c(100, 200, 300), # in the end: 100, 200, 300
           r %in% c(1, 5, 10), 
           d0 %in% c(0.75, 1, 1.25, 1.75))


tab3b <- tab3b[with(tab3b, order(n, r)),]
tab3b$n[c(1:36)%in%c(2:12, 14:24, 26:36)] <- NA
tab3b$r[!c(1:36)%in%seq(1, 36, by = 4)] <- NA

tabjoint <- cbind(tab3, tab3b[, -(1:3)])
table <- tabjoint%>% 
    xtable::xtable(digits = c(0, 0, 0, 2, rep(3, 16)), include.rownames=FALSE)
print(table, include.rownames = FALSE)



tab3c <- left_join(statistics, statistics3) %>%
    left_join(statistics2)%>%
    left_join(statistics4)%>%
    select("n",  "r", "d0",
           "eta_ML", "eta_i1_ML",
           "eps_ML", "eps_i1_ML",
           "rho_ML", "rho_i1_ML", 
           "ar1_ML", "ar1_i1_ML",
           "ar2_ML", "ar2_i1_ML") %>%
    filter(n %in% c(100, 200, 300), # in the end: 100, 200, 300
           r %in% c(1, 5, 10), 
           d0 %in% c(0.75, 1, 1.25, 1.75))


tab3c <- tab3c[with(tab3c, order(n, r)),]
tab3c$n[c(1:36)%in%c(2:12, 14:24, 26:36)] <- NA
tab3c$r[!c(1:36)%in%seq(1, 36, by = 4)] <- NA


tab3d <- left_join(statisticsb, statistics3b) %>%
    left_join(statistics2b)%>%
    left_join(statistics4b)%>%
    select("n",  "r", "d0",
           "eta_ML", "eta_i1_ML",
           "eps_ML", "eps_i1_ML",
           "rho_ML", "rho_i1_ML", 
           "ar1_ML", "ar1_i1_ML",
           "ar2_ML", "ar2_i1_ML") %>%
    filter(n %in% c(100, 200, 300), # in the end: 100, 200, 300
           r %in% c(1, 5, 10), 
           d0 %in% c(0.75, 1, 1.25, 1.75))


tab3d <- tab3d[with(tab3d, order(n, r)),]
tab3d$n[c(1:36)%in%c(2:12, 14:24, 26:36)] <- NA
tab3d$r[!c(1:36)%in%seq(1, 36, by = 4)] <- NA

tabjoint <- cbind(tab3c, tab3d[, -(1:3)])
table <- tabjoint%>% 
    xtable::xtable(digits = c(0, 0, 0, 2, rep(3, 20)), include.rownames=FALSE)
print(table, include.rownames = FALSE)




### Check goodness of fit for the different models
# CSS
load(file = "./MC/results_cor.RData")
X.hat.CSS <- X.res
C.hat.CSS <- C.res
# ML
load(file = "./MC/results_ML_corr.RData")
X.hat.ML <- X.res
C.hat.ML <- C.res
# CSS i1
load(file = "./MC/results_i1_corr.RData")
X.hat.i1 <- X.res
C.hat.i1 <- C.res
# CSS i1 ML
load(file = "./MC/results_i1_ML_corr.RData")
X.hat.i1.ML <- X.res
C.hat.i1.ML <- C.res
ar <- c(-1.6, 0.8)

X.true <- C.true <- list()
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
    
    
    X.true[[i]] <- x
    C.true[[i]] <- c
    cat("Iteration ", i, "\n")
}


calcRSQ <- function(x.true, x.est){
    Rsq <- mean(sapply(1:ncol(x.true), function(j) tryCatch(summary(lm(x.true[, j] ~ x.est[,j]))$r.squared,
                                                            error = function(e) return(NA))), na.rm = TRUE)
}

Rsq.CSS <- sapply(1:nrow(setups), function(j) calcRSQ(X.true[[j]], X.hat.CSS[[j]]))
Rsq.ML <- sapply(1:nrow(setups), function(j) calcRSQ(X.true[[j]], X.hat.ML[[j]]))
Rsq.i1 <- sapply(1:nrow(setups), function(j) calcRSQ(X.true[[j]], X.hat.i1[[j]]))
Rsq.i1.ML <- sapply(1:nrow(setups), function(j) calcRSQ(X.true[[j]], X.hat.i1.ML[[j]]))
Rsq.CSS.C <- sapply(1:nrow(setups), function(j) calcRSQ(C.true[[j]], C.hat.CSS[[j]]))
Rsq.ML.C <- sapply(1:nrow(setups), function(j) calcRSQ(C.true[[j]], C.hat.ML[[j]]))
Rsq.i1.C <- sapply(1:nrow(setups), function(j) calcRSQ(C.true[[j]], C.hat.i1[[j]]))
Rsq.i1.ML.C <- sapply(1:nrow(setups), function(j) calcRSQ(C.true[[j]], C.hat.i1.ML[[j]]))
tab4 <- cbind(setups[, c("n", "r", "d")], Rsq.CSS, Rsq.ML, Rsq.i1, Rsq.i1.ML, 
              Rsq.CSS.C, Rsq.ML.C, Rsq.i1.C, Rsq.i1.ML.C)

tab4 <- tab4[with(tab4, order(n, r)),]
tab4$n[c(1:36)%in%c(2:12, 14:24, 26:36)] <- NA
tab4$r[!c(1:36)%in%seq(1, 36, by = 4)] <- NA


table <- tab4%>% 
    xtable::xtable(digits = c(0, 0, 0, 2, rep(3, 8)), include.rownames=FALSE)
print(table, include.rownames = FALSE)