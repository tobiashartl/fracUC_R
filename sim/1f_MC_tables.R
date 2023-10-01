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


# Load CSS esimates
results.list <- list()
statistics <- statisticsb <- matrix(NA, nrow(setups), 11)
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    try({
        load(paste("./MC/Sim_R1000_n", true.par[1], "_d", true.par[2], "_nu", true.par[3],
                   "_ar_high_corr0.RData", sep=""))
        results.list[[i]] <- RESULTS
        MSE_d <- rowMeans((results.list[[i]][c("d", "d_45", "d_50", "d_55", "d_60", "d_65", "d_70"), ] - true.par$d)^2) %>%
            sqrt(.)
        MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8))^2)%>%
            sqrt(.)
        #SSR   <- mean(results.list[[i]][c("SSR"), ])
        Rsq   <- mean(results.list[[i]][c("Rsq"), ])
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - true.par$r)^2) %>%
            sqrt(.)
        statistics[i, ] <- c(MSE_d, MSE_nu, Rsq, MSE_ar)
        
        MSE_d <- apply((results.list[[i]][c("d", "d_45", "d_50", "d_55", "d_60", "d_65", "d_70"), ] - true.par$d), 1, mean) 
        MSE_ar <- -apply((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8)), 1, mean)
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - true.par$r)) 
        statisticsb[i, ] <- c(MSE_d, MSE_nu, Rsq, MSE_ar)
        
    })
    
}

rownames(statistics) <- rownames(statisticsb) <- apply(setups[,], 1, function(x) paste(x, collapse="_"))
colnames(statistics) <- colnames(statisticsb) <- c("d",
                                                   "d_45", "d_50", "d_55", "d_60", "d_65", "d_70",
                                                   "nu", "Rsq", "ar1", "ar2" )


# Load ML esimates
results.list <- list()
statistics2 <- statistics2b <- matrix(NA, nrow(setups), 5)
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    try({
        load(paste("./MC/ML/Sim_R1000_n", true.par[1], "_d", true.par[2], "_nu", true.par[3],
                   "_ar_high_corr0.RData", sep=""))
        results.list[[i]] <- RESULTS
        MSE_d <- mean((results.list[[i]][c("d"), ] - true.par$d)^2) %>%
            sqrt(.)
        MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8))^2)%>%
            sqrt(.)
        #SSR   <- mean(results.list[[i]][c("SSR"), ])
        Rsq   <- mean(results.list[[i]][c("Rsq"), ])
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - true.par$r)^2) %>%
            sqrt(.)
        statistics2[i, ] <- c(MSE_d, MSE_nu, Rsq, MSE_ar)
        
        
        MSE_d <- mean((results.list[[i]][c("d"), ] - true.par$d)) 
        MSE_ar <- -apply((results.list[[i]][c("ar_1", "ar_2"), ] - c(-1.6, 0.8)), 1, mean)
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - true.par$r)) 
        statistics2b[i, ] <- c(MSE_d, MSE_nu, Rsq, MSE_ar)
    })
    
}

rownames(statistics2) <- rownames(statistics2b) <- apply(setups[,], 1, function(x) paste(x, collapse="_"))
colnames(statistics2) <- colnames(statistics2b) <- c("d_ML",
                                                     "nu_ML", "Rsq_ML", "ar1_ML", "ar2_ML" )


# Load I(1) CSS estimates
results.list <- list()
statistics3 <- statistics3b <- matrix(NA, nrow(setups), 5)
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    try({
        load(paste("./MC/INTEGER/Sim_R1000_n", true.par[1], "_d", true.par[2], "_nu", true.par[3],
                   "_ar_high.RData", sep=""))
        results.list[[i]] <- RESULTS
        SSR   <- mean(results.list[[i]][c("SSR"), ])
        Rsq   <- mean(results.list[[i]][c("Rsq"), ])
        MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(1.6, -0.8))^2)%>%
            sqrt(.)
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - true.par$r)^2) %>% sqrt(.)
        statistics3[i, ] <- c(MSE_nu, SSR, Rsq, MSE_ar)
        
        MSE_ar <- apply((results.list[[i]][c("ar_1", "ar_2"), ] - c(1.6, -0.8)), 1, mean)
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - true.par$r)) 
        statistics3b[i, ] <- c(MSE_nu, SSR, Rsq, MSE_ar)
        
    })
    
}
rownames(statistics3) <- rownames(statistics3b) <- apply(setups[1:nrow(setups),], 1, function(x) paste(x, collapse="_"))
colnames(statistics3) <- colnames(statistics3b) <- c("nu_i1", "SSR_i1", "Rsq_i1", "ar1_i1", "ar2_i1")


# Load I(1) ML estimates
results.list <- list()
statistics4 <- statistics4b <- matrix(NA, nrow(setups), 5)
for(i in 1:nrow(setups)){
    true.par          <- setups[i, ]
    try({
        load(paste("./MC/INTEGER/ML/Sim_R1000_n", true.par[1], "_d", true.par[2], "_nu", true.par[3],
                   "_ar_high_ML.RData", sep=""))
        results.list[[i]] <- RESULTS
        SSR   <- mean(results.list[[i]][c("SSR"), ])
        Rsq   <- mean(results.list[[i]][c("Rsq"), ])
        MSE_ar <- rowMeans((results.list[[i]][c("ar_1", "ar_2"), ] - c(1.6, -0.8))^2)%>%
            sqrt(.)
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - true.par$r)^2) %>% sqrt(.)
        statistics4[i, ] <- c(MSE_nu, SSR, Rsq, MSE_ar)
        
        MSE_ar <- apply((results.list[[i]][c("ar_1", "ar_2"), ] - c(1.6, -0.8)), 1, mean)
        MSE_nu <- mean((results.list[[i]][c("nu"), ] - true.par$r)) 
        statistics4b[i, ] <- c(MSE_nu, SSR, Rsq, MSE_ar)
        
    })
    
}


rownames(statistics4) <- rownames(statistics4b) <- apply(setups[1:nrow(setups),], 1, function(x) paste(x, collapse="_"))
colnames(statistics4) <- colnames(statistics4b) <- c("nu_i1_ML", "SSR_i1_ML", "Rsq_i1_ML", "ar1_i1_ML", "ar2_i1_ML")




statistics <- statistics %>% as.data.frame%>%tibble::rownames_to_column()
statistics2 <- statistics2 %>% as.data.frame%>%tibble::rownames_to_column()
statistics3 <- statistics3 %>% as.data.frame %>% tibble::rownames_to_column()
statistics4 <- statistics4 %>% as.data.frame %>% tibble::rownames_to_column()


statistics$n <- setups$n
statistics$d0 <- setups$d
statistics$r   <- setups$r

### Table 1: RMSE d
tab1 <- left_join(statistics, statistics3) %>%
    left_join(statistics2)%>%
    left_join(statistics4)%>%
    dplyr::select("n",  "r", "d0",
                  "d", "d_ML", "d_50", "d_55","d_60", "d_65","d_70") %>%
    filter(n %in% c(100, 200, 300), # in the end: 100, 200, 300
           r %in% c(1, 5, 10), 
           d0 %in% c(0.75, 1, 1.25, 1.75))

tab1 <- tab1[with(tab1, order(n, r)),]
tab1$n[c(1:36)%in%c(2:12, 14:24, 26:36)] <- NA
tab1$r[!c(1:36)%in%seq(1, 36, by = 4)] <- NA
table <- tab1%>% 
    xtable::xtable(digits = c(0, 0, 0, 2, rep(3, 7)), include.rownames=FALSE)
print(table, include.rownames = FALSE)



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
    dplyr::select("n",  "r", "d0",
                  "d", "d_ML", "d_50", "d_55","d_60", "d_65","d_70") %>%
    filter(n %in% c(100, 200, 300), # in the end: 100, 200, 300
           r %in% c(1, 5, 10), 
           d0 %in% c(0.75, 1, 1.25, 1.75))

tab2 <- tab2[with(tab2, order(n, r)),]
tab2$n[c(1:36)%in%c(2:12, 14:24, 26:36)] <- NA
tab2$r[!c(1:36)%in%seq(1, 27, by = 4)] <- NA

tabjoint <- cbind(tab1, tab2[, -(1:3)])
table <- tabjoint%>% 
    xtable::xtable(digits = c(0, 0, 0, 2, rep(3, 14)), include.rownames=FALSE)
print(table, include.rownames = FALSE)



tab3 <- left_join(statistics, statistics3) %>%
    left_join(statistics2)%>%
    left_join(statistics4)%>%
    dplyr::select("n",  "r", "d0",
                  "nu", "nu_ML", "nu_i1","nu_i1_ML",
                  "ar1", "ar1_ML", "ar1_i1", "ar1_i1_ML",
                  "ar2", "ar2_ML", "ar2_i1", "ar2_i1_ML") %>%
    filter(n %in% c(100, 200, 300), # in the end: 100, 200, 300
           r %in% c(1, 5, 10), 
           d0 %in% c(0.75, 1, 1.25, 1.75))


tab3 <- tab3[with(tab3, order(n, r)),]
tab3$n[c(1:36)%in%c(2:12, 14:24, 26:36)] <- NA
tab3$r[!c(1:36)%in%seq(1, 36, by = 4)] <- NA


table <- tab3%>% 
    xtable::xtable(digits = c(0, 0, 0, 2, c(3, 3, 3, 3), rep(3, 8)), include.rownames=FALSE)
print(table, include.rownames = FALSE)


tab3 <- left_join(statisticsb, statistics3b) %>%
    left_join(statistics2b)%>%
    left_join(statistics4b)%>%
    dplyr::select("n",  "r", "d0",
                  "nu", "nu_ML", "nu_i1","nu_i1_ML",
                  "ar1", "ar1_ML", "ar1_i1", "ar1_i1_ML",
                  "ar2", "ar2_ML", "ar2_i1", "ar2_i1_ML") %>%
    filter(n %in% c(100, 200, 300), # in the end: 100, 200, 300
           r %in% c(1, 5, 10), 
           d0 %in% c(0.75, 1, 1.25, 1.75))


tab3 <- tab3[with(tab3, order(n, r)),]
tab3$n[c(1:36)%in%c(2:12, 14:24, 26:36)] <- NA
tab3$r[!c(1:36)%in%seq(1, 36, by = 4)] <- NA


table <- tab3%>% 
    xtable::xtable(digits = c(0, 0, 0, 2, c(3, 3, 3, 3), rep(3, 8)), include.rownames=FALSE)
print(table, include.rownames = FALSE)



### Check goodness of fit for the different models
# CSS
load(file = "./MC/results.RData")
X.hat.CSS <- X.res
C.hat.CSS <- C.res
# ML
load(file = "./MC/results_ML.RData")
X.hat.ML <- X.res
C.hat.ML <- C.res
# CSS i1
load(file = "./MC/results_i1.RData")
X.hat.i1 <- X.res
C.hat.i1 <- C.res
# CSS i1 ML
load(file = "./MC/results_i1_ML.RData")
X.hat.i1.ML <- X.res
C.hat.i1.ML <- C.res





calcRSQ <- function(x.true, x.est){
    Rsq <- mean(sapply(1:ncol(x.true), function(j) summary(lm(x.true[, j] ~ x.est[,j]))$r.squared))
}

Rsq.CSS <- sapply(1:nrow(setups), function(j) calcRSQ(X.true[[j]], X.hat.CSS[[j]]))
Rsq.ML <- sapply(1:nrow(setups), function(j) calcRSQ(X.true[[j]], X.hat.ML[[j]]))
Rsq.i1 <- sapply(1:nrow(setups), function(j) calcRSQ(X.true[[j]], X.hat.i1[[j]]))
Rsq.i1.ML <- sapply(1:nrow(setups), function(j) calcRSQ(X.true[[j]], X.hat.i1.ML[[j]]))
Rsq.CSS.C <- sapply(1:nrow(setups), function(j) calcRSQ(C.true[[j]], C.hat.CSS[[j]]))
Rsq.ML.C <- sapply(1:nrow(setups), function(j) calcRSQ(C.true[[j]], C.hat.ML[[j]]))
Rsq.i1.C <- sapply(1:nrow(setups), function(j) calcRSQ(C.true[[j]], C.hat.i1[[j]]))
Rsq.i1.ML.C <- sapply(1:nrow(setups), function(j) calcRSQ(C.true[[j]], C.hat.i1.ML[[j]]))
tab4 <- cbind(setups[, c("n", "r", "d")], Rsq.CSS, Rsq.ML, Rsq.i1, Rsq.i1.ML, Rsq.CSS.C, Rsq.ML.C, Rsq.i1.C, Rsq.i1.ML.C)

tab4 <- tab4[with(tab4, order(n, r)),]
tab4$n[c(1:36)%in%c(2:12, 14:24, 26:36)] <- NA
tab4$r[!c(1:36)%in%seq(1, 36, by = 4)] <- NA


table <- tab4%>% 
    xtable::xtable(digits = c(0, 0, 0, 2, rep(3, 8)), include.rownames=FALSE)
print(table, include.rownames = FALSE)

