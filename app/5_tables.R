# final tables
#Packages
gc()
rm(list = ls())
library(fUCpack)
library(dplyr)
library(CFFpack)
library(xtable)


# Wd, etc
setwd("/Users/tobias/Dokumente/Projekte/filtering unknown persistence/R/code")
source("./help functions/fUC_arma_approx.R")

# load data
data<- read.csv(file = "./app/noaa_SST.csv", skip = 4)[, c(1, 2)] %>%
    as.data.frame() 


y <- data[,2]
p   <- 4
#est.final <- readRDS(file = paste("./Applications/temptest/NOAA_SST_est_", "frac", "_p", p, "_hess.RDS", sep=""))



# lag selection:
#           frac        I(1)        I(2)
# AIC       4           2           2
# BIC       4           4           2




# load estimation results
Est.id <- readRDS(file = "./app/fUC_ML_Results.RDS")
Est.i1 <- readRDS(file = "./app/Benchmark_i1_ML_Results.RDS")
Est.i2 <- readRDS(file = "./app/Benchmark_i2_ML_Results.RDS")

Est.id["d", 1] + c(qnorm(0.975)*Est.id["d", 2], -qnorm(0.975)*Est.id["d", 2])
Est.id["d", 1] + c(qnorm(0.995)*Est.id["d", 2], -qnorm(0.995)*Est.id["d", 2])
tab.all <- data.frame(est.id = Est.id[,1],
                      se.id = Est.id[,2],
                      est.i1 = Est.i1[,1],
                      se.i1 = Est.i1[,2],
                      est.i2 = Est.i2[,1],
                      se.i2 = Est.i2[,2]
)

tab.all["ll", ] <- - tab.all["ll", ]

colnames(tab.all) = c("Estimate", "Std. Error", "Estimate", "Std. Error", "Estimate", "Std. Error")
mdat <- matrix(c(4, rep(4,3),rep(3,4), rep(3, (3)), 5, 5, 6, 6),
               nrow = 15, ncol=7, byrow=F)
tab <- xtable::xtable(tab.all, display = c("s",rep("G", ncol(tab.all))), digits=mdat)

print.xtable(tab, type = "latex")






# figures
TC.id <- readRDS(file = "./app/fUC_ML_TC.RDS")
TC.i1 <- readRDS(file = "./app/Benchmark_i1_ML_TC.RDS")
TC.i2 <- readRDS(file = "./app/Benchmark_i2_ML_TC.RDS")

hp <- mFilter::hpfilter(y, freq = 14400, type = "lambda")

# Trend
data.plot <- data.frame(
    y = TC.id$y,
    trend.id  = TC.id$trend,
    cycle.id  = TC.id$cycle,
    trend.i1  = TC.i1$trend,
    cycle.i1  = TC.i1$cycle,
    trend.i2  = TC.i2$trend,
    cycle.i2  = TC.i2$cycle,
    trend.hp  = hp$trend,
    cycle.hp  = y - hp$trend, 
    time      = TC.id$time
)






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
    geom_line(mapping = aes(x = time, y = trend.id), color = cbbPalette[1], linetype = "twodash") +
    geom_line(mapping = aes(x = time, y = trend.i1), color = cbbPalette[2], linetype = "longdash") +
    geom_line(mapping = aes(x = time, y = trend.i2), color = cbbPalette[3], linetype = "F1") +
    geom_line(mapping = aes(x = time, y = trend.hp), color = "#af8dc3", linetype = "dashed") +
    
    #    geom_line(mapping = aes(x = time, y = tau), color = cbbPalette[1], linetype = "twodash") +
    theme_classic()
#gg1


gg2 <- ggplot(data.plot, aes(x = time, y = cycle.i1 - cycle.id)) +
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
    labs(x = "(b)", y="")+
    geom_line(color='black', linetype = "solid") + 
    ylim(-.25, .2) + 
    #ggtitle(bquote("Cyclical temperature anomalies")) + 
    #geom_line(mapping = aes(x = time, y = cycle.i1), color = cbbPalette[2], linetype = "longdash") +
    #geom_line(mapping = aes(x = time, y = cycle.i2), color = cbbPalette[3], linetype = "F1") +
    #geom_line(mapping = aes(x = time, y = cycle.hp), color = "#af8dc3", linetype = "dashed") +
    
    #    geom_line(mapping = aes(x = time, y = tau), color = cbbPalette[1], linetype = "twodash") +
    theme_classic()
#gg1
gg3 <- ggplot(data.plot, aes(x = time, y = cycle.i2 - cycle.id)) +
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
    labs(x = "(c)", y="")+
    geom_line(color='black', linetype = "solid") + 
    ylim(-.25, .2) + 
    #    geom_line(mapping = aes(x = time, y = tau), color = cbbPalette[1], linetype = "twodash") +
    theme_classic()
#gg1

gg4 <- ggplot(data.plot, aes(x = time, y = cycle.hp - cycle.id)) +
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
    labs(x = "(d)", y="")+
    geom_line(color='black', linetype = "solid") + 
    ylim(-.25, .2) + 
    # ggtitle(bquote("Cyclical temperature anomalies")) + 
    # geom_line(mapping = aes(x = time, y = cycle.i1), color = cbbPalette[2], linetype = "longdash") +
    # geom_line(mapping = aes(x = time, y = cycle.i2), color = cbbPalette[3], linetype = "F1") +
    # geom_line(mapping = aes(x = time, y = cycle.hp), color = "#af8dc3", linetype = "dashed") +
    # 
    #    geom_line(mapping = aes(x = time, y = tau), color = cbbPalette[1], linetype = "twodash") +
    theme_classic()
#gg1


gg5 <- ggplot(data.plot, aes(x = time, y = cycle.id)) +
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
    labs(x = "(a)", y="")+
    geom_line(color='black', linetype = "solid") + 
    ggtitle(bquote("Cyclical temperature anomalies")) + 
    # geom_line(mapping = aes(x = time, y = cycle.i1), color = cbbPalette[2], linetype = "longdash") +
    # geom_line(mapping = aes(x = time, y = cycle.i2), color = cbbPalette[3], linetype = "F1") +
    # geom_line(mapping = aes(x = time, y = cycle.hp), color = "#af8dc3", linetype = "dashed") +
    # 
    #    geom_line(mapping = aes(x = time, y = tau), color = cbbPalette[1], linetype = "twodash") +
    theme_classic()
#gg1



KF_id <- fUC_comp(y - TC.id$det, Est.id[1,1], matrix(c(Est.id[2:3, 1], Est.id[3:4, 1]), 2, 2),
                  Est.id[5:8,1], corr=TRUE)

# checks
acf_v <- acf(KF_id$v, type = "correlation", plot = T, demean=F, lag.max=48)
data.plotting <- data.frame(autocorrelation = acf_v$acf[-1],
                            lag = 1:48)

ic_alpha= function(alpha, acf_res){
    return(qnorm((1 + (1 - alpha))/2)/sqrt(acf_res$n.used))
}
hline1 <- ic_alpha(0.05, acf_v)
hline2 <- ic_alpha(0.01, acf_v)

plot_v <- ggplot(data.plotting, aes(x = lag, y = autocorrelation)) +
    labs(x = "Lag", y="ACF")+
    geom_hline(aes(yintercept=0)) + 
    geom_segment(mapping = aes(xend = lag, yend=0))+
    ggtitle(bquote("I(d): Autocorrelation in "~v[t])) + 
    geom_hline(aes(yintercept = hline1), linetype = 2, color = cbbPalette[1])+
    geom_hline(aes(yintercept = -hline1), linetype = 2, color = cbbPalette[1])+
    geom_hline(aes(yintercept = hline2), linetype = 2, color = cbbPalette[2])+
    geom_hline(aes(yintercept = -hline2), linetype = 2, color = cbbPalette[2])+
    ylim(c(-0.06, 0.2))+
    theme_classic()





KF_i1 <- fUC_comp(y - TC.i1$det, 1, matrix(c(Est.i1[2:3, 1], Est.i1[3:4, 1]), 2, 2),
                  Est.i1[5:8,1], corr=TRUE)
acf_vi1 <- acf(KF_i1$v, type = "correlation", plot = T, demean=F, lag.max=48)
data.plotting.i1 <- data.frame(autocorrelation = acf_vi1$acf[-1],
                            lag = 1:48)

hline1 <- ic_alpha(0.05, acf_vi1)
hline2 <- ic_alpha(0.01, acf_vi1)

plot_v.i1 <- ggplot(data.plotting.i1, aes(x = lag, y = autocorrelation)) +
    labs(x = "Lag", y="")+
    geom_hline(aes(yintercept=0)) + 
    geom_segment(mapping = aes(xend = lag, yend=0))+
    ggtitle(bquote("I(1): Autocorrelation in "~v[t])) + 
    geom_hline(aes(yintercept = hline1), linetype = 2, color = cbbPalette[1])+
    geom_hline(aes(yintercept = -hline1), linetype = 2, color = cbbPalette[1])+
    geom_hline(aes(yintercept = hline2), linetype = 2, color = cbbPalette[2])+
    geom_hline(aes(yintercept = -hline2), linetype = 2, color = cbbPalette[2])+
    ylim(c(-0.06, 0.2))+
    theme_classic()

KF_i2 <- fUC_comp(y - TC.i2$det, 2, matrix(c(Est.i2[2:3, 1], Est.i2[3:4, 1]), 2, 2),
                  Est.i2[5:8,1], corr=TRUE)
acf_vi2 <- acf(KF_i2$v, type = "correlation", plot = T, demean=F, lag.max=48)
data.plotting.i2 <- data.frame(autocorrelation = acf_vi2$acf[-1],
                               lag = 1:48)

hline1 <- ic_alpha(0.05, acf_vi2)
hline2 <- ic_alpha(0.01, acf_vi2)

plot_v.i2 <- ggplot(data.plotting.i2, aes(x = lag, y = autocorrelation)) +
    labs(x = "Lag", y="")+
    geom_hline(aes(yintercept=0)) + 
    geom_segment(mapping = aes(xend = lag, yend=0))+
    ggtitle(bquote("I(2): Autocorrelation in "~v[t])) + 
    geom_hline(aes(yintercept = hline1), linetype = 2, color = cbbPalette[1])+
    geom_hline(aes(yintercept = -hline1), linetype = 2, color = cbbPalette[1])+
    geom_hline(aes(yintercept = hline2), linetype = 2, color = cbbPalette[2])+
    geom_hline(aes(yintercept = -hline2), linetype = 2, color = cbbPalette[2])+
    ylim(c(-0.06, 0.2))+
    theme_classic()
ggcbind <- ggarrange(plot_v, plot_v.i1, plot_v.i2, ncol = 3)


