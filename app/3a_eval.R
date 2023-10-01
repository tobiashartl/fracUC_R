#Packages
gc()
rm(list = ls())
library(fUCpack)
library(dplyr)
library(parallel)
library(CFFpack)


# Wd, etc
setwd("/Users/tobias/Dokumente/Projekte/filtering unknown persistence/R/code")
source("./help functions/fUC_arma_approx.R")

# load data
data<- read.csv(file = "./app/noaa_SST.csv", skip = 4)[, c(1, 2)] %>%
    as.data.frame() 


y <- data[,2]
corr <- TRUE
START = 2
diffuse = T
det <- "frac"
p   <- 4
est.final <- readRDS(file = paste("./app/NOAA_SST_est_", det, "_p", p, "_hess.RDS", sep=""))


# =============================================================================
# Parameter analysis
# =============================================================================
theta <- est.final$par
(d <- theta[1])
(Q <- mlogvech2mat(theta[2:4])) %>% cov2cor
(ar   <- est.final$par[-(1:4)]) 
toComp(-ar)


H <- est.final$hessian
J <- numDeriv::jacobian(mlogvech2mat, est.final$par[2:4])
cov0 <- solve(est.final$hessian)
Cov <- bdiag(1, J[-2,], diag(p))%*%cov0%*%t(bdiag(1, J[-2, ], diag(p)))
theta.trans <- c(est.final$par[1], mlogvech2mat(est.final$par[2:4])[c(TRUE, TRUE, FALSE, TRUE)], est.final$par[-(1:4)])
se.trans <- sqrt(abs(diag(Cov)))
EST <- cbind(theta.trans, se.trans) 
colnames(EST) <- c("par", "se") 
rownames(EST) <- c("d", "Q11", "Q21", "Q22", paste("ar_", 1:p, sep=""))
EST
EST[,1] / EST[,2]
EST[1,1] + c(-qt(.975, 2083)*EST[1,2], qt(.975, 2083)*EST[1,2])
(EST[4,1]/EST[2,1])


ar_coef <- -EST[5:(4+p),1]
plot(ARMAtoMA(ar_coef, lag.max=100))
toComp(ar_coef)
# add further relevant terms
det      <- fUC_opt_ML(theta, y, nulim = c(0, Inf), 
                       corr=TRUE, deterministics = "frac", START = 1, 
                       diffuse = TRUE, return.det = TRUE, nu.opt = FALSE)
det.trend <- frac_diff(rep(1, length(y)), -theta[1])*det

k <- length(est.final$par)+1
(BIC          <- k * log(length(y)) + 2*est.final$value)
(AIC          <- 2 * k              + 2*est.final$value)

SSR <- (fUC_comp((y - det.trend), d, Q, ar, corr=TRUE)$v[-1])^2 %>% sum(.)
add <- c(Q[2,2]/Q[1,1], Q[2,1]/Q[1,1], cov2cor(Q)[2,1], est.final$value, SSR, AIC, BIC)
EST <- rbind(EST, cbind(add, NA))
rownames(EST) <- c("d", "Q11", "Q21", "Q22", paste("ar_", 1:p, sep=""),
                   "nu_1", "nu_2", "rho", "ll", "CSS", "AIC", "BIC")
colnames(EST)<- c("par", "se")
saveRDS(EST, file = "./app/fUC_ML_Results.RDS")



# =============================================================================
# Trend-cycle decomposition
# =============================================================================


plot(y, type="l")
lines(det.trend, col="2")


TC      <- fUC_smooth(y-det.trend, d, Q, ar, corr=TRUE)
plot(y, type="l")
lines(det.trend + TC$x, col="2")


# forecasting exercise: 
t.end <- 2083 + 80*12
x.pred <- c(TC$x, rep(NA, 50*12))
coefs  <- -frac_diff(c(1, rep(0, t.end)), EST[1,1])[-1]
for(tt in 2084:t.end){
    x.pred[tt] <- x.pred[(tt-1):1] %*% coefs[1:(tt-1)]
}

det.trend.pred <- frac_diff(rep(1, t.end), -EST[1,1])*det
plot(det.trend.pred)
y_longrun <- det.trend.pred + x.pred

plot(ts(y_longrun, start = c(1850, 1), frequency = 12))
lines(ts(det.trend.pred, start = c(1850, 1), frequency=12), col=3)
lines(ts(y, start = c(1850, 1), freq = 12), col="2")

data.plot <- data.frame(
    y=y,
    trend = TC$x + det.trend, 
    x = TC$x, 
    det = det.trend, 
    cycle = TC$c,
    time = seq(from = as.Date("1850-01-01"), to = as.Date("2023-07-01"),  by = "month"))

saveRDS(data.plot, file = "./app/fUC_ML_TC.RDS")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbbPalette <- wesanderson::wes_palette("Darjeeling1", n=5)

# nino periods: 1950 onwards: https://origin.cpc.ncep.noaa.gov/products/analysis_monitoring/ensostuff/ONI_v5.php
# before: 

library(ggplot2)
gg1 <- ggplot(data.plot, aes(x = time, y = y)) +
    # la nina
    geom_rect(xmin=as.Date("1871-06-01"), xmax=as.Date("1871-07-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1871-12-01"), xmax=as.Date("1873-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1874-01-01"), xmax=as.Date("1876-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1878-10-01"), xmax=as.Date("1879-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1879-05-01"), xmax=as.Date("1879-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1879-10-01"), xmax=as.Date("1880-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1881-08-01"), xmax=as.Date("1881-09-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1881-12-01"), xmax=as.Date("1882-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1882-03-01"), xmax=as.Date("1882-07-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1883-01-01"), xmax=as.Date("1883-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1886-05-01"), xmax=as.Date("1887-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1889-07-01"), xmax=as.Date("1891-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1892-03-01"), xmax=as.Date("1895-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1898-02-01"), xmax=as.Date("1898-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1898-08-01"), xmax=as.Date("1898-10-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1899-01-01"), xmax=as.Date("1899-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1901-11-01"), xmax=as.Date("1901-11-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1903-08-01"), xmax=as.Date("1903-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1903-11-01"), xmax=as.Date("1904-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1906-09-01"), xmax=as.Date("1907-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1908-03-01"), xmax=as.Date("1911-07-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1912-07-01"), xmax=as.Date("1912-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1913-04-01"), xmax=as.Date("1913-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1915-11-01"), xmax=as.Date("1915-12-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1916-02-01"), xmax=as.Date("1918-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1921-03-01"), xmax=as.Date("1921-10-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1922-02-01"), xmax=as.Date("1922-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1922-12-01"), xmax=as.Date("1923-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1924-05-01"), xmax=as.Date("1925-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1933-06-01"), xmax=as.Date("1934-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1938-04-01"), xmax=as.Date("1939-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1942-07-01"), xmax=as.Date("1943-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1945-08-01"), xmax=as.Date("1945-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1947-09-01"), xmax=as.Date("1947-10-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    # la nina new
    geom_rect(xmin=as.Date("1949-08-01"), xmax=as.Date("1950-07-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1950-11-01"), xmax=as.Date("1951-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1954-05-01"), xmax=as.Date("1954-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1954-07-01"), xmax=as.Date("1956-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1964-04-01"), xmax=as.Date("1965-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1968-01-01"), xmax=as.Date("1968-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1970-07-01"), xmax=as.Date("1972-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1973-05-01"), xmax=as.Date("1974-07-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1974-10-01"), xmax=as.Date("1976-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1981-02-01"), xmax=as.Date("1982-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1983-10-01"), xmax=as.Date("1984-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1984-05-01"), xmax=as.Date("1984-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1984-10-01"), xmax=as.Date("1985-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1988-05-01"), xmax=as.Date("1989-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1995-08-01"), xmax=as.Date("1996-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1997-01-01"), xmax=as.Date("1997-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1998-07-01"), xmax=as.Date("2001-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2005-11-01"), xmax=as.Date("2006-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2007-07-01"), xmax=as.Date("2008-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2008-11-01"), xmax=as.Date("2009-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2010-06-01"), xmax=as.Date("2011-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2011-08-01"), xmax=as.Date("2012-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2016-08-01"), xmax=as.Date("2016-12-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2017-10-01"), xmax=as.Date("2018-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2020-08-01"), xmax=as.Date("2021-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2021-09-01"), xmax=as.Date("2023-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    
    # el nino hist
    geom_rect(xmin=as.Date("1877-06-01"), xmax=as.Date("1878-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1884-05-01"), xmax=as.Date("1884-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1884-09-01"), xmax=as.Date("1884-09-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1885-01-01"), xmax=as.Date("1885-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1885-05-01"), xmax=as.Date("1885-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1885-09-01"), xmax=as.Date("1885-12-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1888-03-01"), xmax=as.Date("1889-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1896-06-01"), xmax=as.Date("1897-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1899-08-01"), xmax=as.Date("1900-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1901-01-01"), xmax=as.Date("1901-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1902-05-01"), xmax=as.Date("1903-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1904-08-01"), xmax=as.Date("1906-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1911-10-01"), xmax=as.Date("1912-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1913-09-01"), xmax=as.Date("1913-09-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1913-11-01"), xmax=as.Date("1914-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1914-07-01"), xmax=as.Date("1915-07-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1918-07-01"), xmax=as.Date("1919-10-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1920-02-01"), xmax=as.Date("1920-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1923-08-01"), xmax=as.Date("1924-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1925-08-01"), xmax=as.Date("1926-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1929-07-01"), xmax=as.Date("1929-10-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1930-02-01"), xmax=as.Date("1931-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1932-04-01"), xmax=as.Date("1932-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1939-08-01"), xmax=as.Date("1939-10-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1940-01-01"), xmax=as.Date("1942-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1948-03-01"), xmax=as.Date("1948-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    # el nino
    geom_rect(xmin=as.Date("1951-06-01"), xmax=as.Date("1952-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1953-02-01"), xmax=as.Date("1954-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1957-04-01"), xmax=as.Date("1958-07-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1958-11-01"), xmax=as.Date("1959-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1963-06-01"), xmax=as.Date("1964-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1965-06-01"), xmax=as.Date("1966-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1968-07-01"), xmax=as.Date("1968-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1968-10-01"), xmax=as.Date("1969-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1969-08-01"), xmax=as.Date("1970-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1972-05-01"), xmax=as.Date("1973-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1976-09-01"), xmax=as.Date("1977-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1977-09-01"), xmax=as.Date("1978-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1979-11-01"), xmax=as.Date("1980-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1982-05-01"), xmax=as.Date("1983-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1986-09-01"), xmax=as.Date("1988-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1991-06-01"), xmax=as.Date("1992-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1993-04-01"), xmax=as.Date("1993-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1994-09-01"), xmax=as.Date("1995-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1997-05-01"), xmax=as.Date("1998-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2002-06-01"), xmax=as.Date("2003-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2004-08-01"), xmax=as.Date("2005-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2006-09-01"), xmax=as.Date("2007-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2009-08-01"), xmax=as.Date("2010-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2014-11-01"), xmax=as.Date("2015-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2015-03-01"), xmax=as.Date("2016-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2018-10-01"), xmax=as.Date("2019-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2019-11-01"), xmax=as.Date("2019-12-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2023-06-01"), xmax=as.Date("2023-07-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    labs(x = "time", y="")+
    geom_line(color="black") + 
    ggtitle(bquote("Trend temperature anomalies")) + 
    geom_line(mapping = aes(x = time, y = trend), color = cbbPalette[1], linetype = "twodash") +
    #    geom_line(mapping = aes(x = time, y = tau), color = cbbPalette[1], linetype = "twodash") +
    theme_classic()
#gg1


gg2 <- ggplot(data.plot, aes(x = time, y = cycle)) +
    # la nina
    geom_rect(xmin=as.Date("1871-06-01"), xmax=as.Date("1871-07-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1871-12-01"), xmax=as.Date("1873-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1874-01-01"), xmax=as.Date("1876-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1878-10-01"), xmax=as.Date("1879-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1879-05-01"), xmax=as.Date("1879-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1879-10-01"), xmax=as.Date("1880-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1881-08-01"), xmax=as.Date("1881-09-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1881-12-01"), xmax=as.Date("1882-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1882-03-01"), xmax=as.Date("1882-07-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1883-01-01"), xmax=as.Date("1883-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1886-05-01"), xmax=as.Date("1887-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1889-07-01"), xmax=as.Date("1891-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1892-03-01"), xmax=as.Date("1895-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1898-02-01"), xmax=as.Date("1898-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1898-08-01"), xmax=as.Date("1898-10-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1899-01-01"), xmax=as.Date("1899-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1901-11-01"), xmax=as.Date("1901-11-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1903-08-01"), xmax=as.Date("1903-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1903-11-01"), xmax=as.Date("1904-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1906-09-01"), xmax=as.Date("1907-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1908-03-01"), xmax=as.Date("1911-07-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1912-07-01"), xmax=as.Date("1912-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1913-04-01"), xmax=as.Date("1913-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1915-11-01"), xmax=as.Date("1915-12-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1916-02-01"), xmax=as.Date("1918-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1921-03-01"), xmax=as.Date("1921-10-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1922-02-01"), xmax=as.Date("1922-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1922-12-01"), xmax=as.Date("1923-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1924-05-01"), xmax=as.Date("1925-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1933-06-01"), xmax=as.Date("1934-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1938-04-01"), xmax=as.Date("1939-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1942-07-01"), xmax=as.Date("1943-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1945-08-01"), xmax=as.Date("1945-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1947-09-01"), xmax=as.Date("1947-10-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    # la nina new
    geom_rect(xmin=as.Date("1949-08-01"), xmax=as.Date("1950-07-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1950-11-01"), xmax=as.Date("1951-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1954-05-01"), xmax=as.Date("1954-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1954-07-01"), xmax=as.Date("1956-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1964-04-01"), xmax=as.Date("1965-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1968-01-01"), xmax=as.Date("1968-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1970-07-01"), xmax=as.Date("1972-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1973-05-01"), xmax=as.Date("1974-07-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1974-10-01"), xmax=as.Date("1976-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1981-02-01"), xmax=as.Date("1982-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1983-10-01"), xmax=as.Date("1984-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1984-05-01"), xmax=as.Date("1984-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1984-10-01"), xmax=as.Date("1985-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1988-05-01"), xmax=as.Date("1989-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1995-08-01"), xmax=as.Date("1996-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1997-01-01"), xmax=as.Date("1997-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("1998-07-01"), xmax=as.Date("2001-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2005-11-01"), xmax=as.Date("2006-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2007-07-01"), xmax=as.Date("2008-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2008-11-01"), xmax=as.Date("2009-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2010-06-01"), xmax=as.Date("2011-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2011-08-01"), xmax=as.Date("2012-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2016-08-01"), xmax=as.Date("2016-12-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2017-10-01"), xmax=as.Date("2018-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2020-08-01"), xmax=as.Date("2021-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    geom_rect(xmin=as.Date("2021-09-01"), xmax=as.Date("2023-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#d1e5f0", alpha=0.15, col="#d1e5f0") +
    
    # el nino hist
    geom_rect(xmin=as.Date("1877-06-01"), xmax=as.Date("1878-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1884-05-01"), xmax=as.Date("1884-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1884-09-01"), xmax=as.Date("1884-09-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1885-01-01"), xmax=as.Date("1885-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1885-05-01"), xmax=as.Date("1885-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1885-09-01"), xmax=as.Date("1885-12-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1888-03-01"), xmax=as.Date("1889-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1896-06-01"), xmax=as.Date("1897-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1899-08-01"), xmax=as.Date("1900-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1901-01-01"), xmax=as.Date("1901-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1902-05-01"), xmax=as.Date("1903-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1904-08-01"), xmax=as.Date("1906-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1911-10-01"), xmax=as.Date("1912-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1913-09-01"), xmax=as.Date("1913-09-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1913-11-01"), xmax=as.Date("1914-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1914-07-01"), xmax=as.Date("1915-07-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1918-07-01"), xmax=as.Date("1919-10-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1920-02-01"), xmax=as.Date("1920-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1923-08-01"), xmax=as.Date("1924-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1925-08-01"), xmax=as.Date("1926-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1929-07-01"), xmax=as.Date("1929-10-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1930-02-01"), xmax=as.Date("1931-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1932-04-01"), xmax=as.Date("1932-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1939-08-01"), xmax=as.Date("1939-10-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1940-01-01"), xmax=as.Date("1942-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1948-03-01"), xmax=as.Date("1948-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    # el nino
    geom_rect(xmin=as.Date("1951-06-01"), xmax=as.Date("1952-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1953-02-01"), xmax=as.Date("1954-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1957-04-01"), xmax=as.Date("1958-07-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1958-11-01"), xmax=as.Date("1959-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1963-06-01"), xmax=as.Date("1964-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1965-06-01"), xmax=as.Date("1966-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1968-07-01"), xmax=as.Date("1968-08-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1968-10-01"), xmax=as.Date("1969-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1969-08-01"), xmax=as.Date("1970-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1972-05-01"), xmax=as.Date("1973-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1976-09-01"), xmax=as.Date("1977-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1977-09-01"), xmax=as.Date("1978-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1979-11-01"), xmax=as.Date("1980-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1982-05-01"), xmax=as.Date("1983-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1986-09-01"), xmax=as.Date("1988-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1991-06-01"), xmax=as.Date("1992-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1993-04-01"), xmax=as.Date("1993-06-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1994-09-01"), xmax=as.Date("1995-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("1997-05-01"), xmax=as.Date("1998-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2002-06-01"), xmax=as.Date("2003-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2004-08-01"), xmax=as.Date("2005-02-28"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2006-09-01"), xmax=as.Date("2007-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2009-08-01"), xmax=as.Date("2010-03-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2014-11-01"), xmax=as.Date("2015-01-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2015-03-01"), xmax=as.Date("2016-04-30"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2018-10-01"), xmax=as.Date("2019-05-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2019-11-01"), xmax=as.Date("2019-12-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    geom_rect(xmin=as.Date("2023-06-01"), xmax=as.Date("2023-07-31"), 
              ymin=-Inf, ymax=Inf, fill="#fddbc7", alpha=0.15, col="#fddbc7") +
    labs(x = "time", y="")+
    geom_line(color="black") + 
    ggtitle(bquote("Cyclical temperature anomalies")) + 
    #    geom_line(mapping = aes(x = time, y = tau), color = cbbPalette[1], linetype = "twodash") +
    theme_classic()
#gg2



