
# This script processes study-level data scraped from the meta-meta-analyses.s

rm(list = ls())

########################### PRELIMINARIES ########################### 

library(dplyr)
library(readxl)
library(robumeta)
library(Hmisc)
library(data.table)
library(testthat)
library(here)
library(MetaUtility)
library(msm)

raw.data.dir = here("Empirical estimates of confounding bias/Raw data")
prepped.data.dir = here("Empirical estimates of confounding bias/Prepped data")
code.dir = here()


########################### BUN 2020 ########################### 

# can't scrape because they don't have study-level data
#  no luck emailing them

########################### GOLDER 2011 ########################### 

# data scraped from page 6 forest plot
# ratios are OR_RCT / OR_obs 
setwd(raw.data.dir)
setwd("Golder data")

dg = read_xlsx("Golder forest plot scraped.xlsx")
expect_equal( nrow(dg), 58 )


dg2 = parse_CI_string(string = dg$string, sep = ",")

# invert so the ratios are OR_obs / OR_RCT
dg2$yi = 1/dg2$yi
dg2$hi = 1/dg2$hi

dg3 = scrape_meta(type = "RR",
            est = dg2$yi,
            hi = dg2$hi)

dg3$study = "Golder 2011"


########################### SHIKATA 2006 ########################### 

setwd(raw.data.dir)
setwd("Shikata data")

# scraped from Table 2 (pg 674)
ds = read_xlsx("Shikata table scraped.xlsx")
expect_equal( nrow(ds), 52 )

# this one reported the ORs for RCTs and observationals separately, so we will 
#  calculate the ratios ourselves
# get variances of log-ratios from RCTs and observationals
#  so we can use the delta method for the ratio
RCT.var = scrape_meta( type = "RR",
                      est = ds$est.RCT,
                      hi = parse_CI_string(string = ds$ci.RCT, sep = "-")$hi )$vyi


obs.var = scrape_meta( type = "RR",
                       est = ds$est.obs,
                       hi = parse_CI_string(string = ds$ci.obs, sep = "-")$hi )$vyi



ds2 = data.frame( est.obs = ds$est.obs,
                  est.RCT = ds$est.RCT,
                  obs.var = obs.var,
                  RCT.var = RCT.var,
                  yi = log( ds$est.obs / ds$est.RCT ) )


# get variance of log-ratio (obs vs. RCT) using delta method
ds2 = ds2 %>% rowwise() %>%
  mutate( vyi = deltamethod( ~ log(x1) - log(x2),
                             mean = c( est.obs, est.RCT ), cov = diag( c(obs.var, RCT.var) ) )^2 )


# sanity check: check first entry, using numbers straight from table
expect_equal( ds2$yi[1],
              log(0.59 / 0.65) )


v1 = scrape_meta( est = 0.59, hi = .85 )$vyi
v2 = scrape_meta( est = 0.65, hi = 1.4 )$vyi

expect_equal( ds2$vyi[1],
              deltamethod( ~ log(x1) - log(x2),
                           mean = c( 0.59, 0.65 ), cov = diag( c(v1, v2) ) )^2 )



ds2$study = "Shikata 2006"

########################### IOANNIDIS 2001 ########################### 

# could digitize the scatterplot, but doesn't have variances


########################### FRANKLIN 2020 (RCT-DUPLICATE) ########################### 

# data provided by authors
setwd(raw.data.dir)
setwd("Franklin data")

dfr = read_xlsx("Franklin table.xlsx")
expect_equal( nrow(dfr), 10 )


dfr2 = parse_CI_string(string = df$string, sep = ",")

est.RCT = parse_CI_string(string = dfr$`RCT result`, sep = ",")$yi
est.obs = parse_CI_string(string = dfr$`RWE results`, sep = ",")$yi

# get variances of log-ratios from RCTs and observationals
#  so we can use the delta method for the ratio
RCT.var = scrape_meta( type = "RR",
                       est = est.RCT,
                       hi = parse_CI_string(string = dfr$`RCT result`, sep = ",")$hi )$vyi


obs.var = scrape_meta( type = "RR",
                       est = est.obs,
                       hi = parse_CI_string(string = dfr$`RWE results`, sep = ",")$hi )$vyi



dfr2 = data.frame( est.obs = est.obs,
                  est.RCT = est.RCT,
                  obs.var = obs.var,
                  RCT.var = RCT.var,
                  yi = log( est.obs / est.RCT ) )


# get variance of log-ratio using delta method
dfr2 = dfr2 %>% rowwise() %>%
  mutate( vyi = deltamethod( ~ log(x1) - log(x2),
                             mean = c( est.obs, est.RCT ), cov = diag( c(obs.var, RCT.var) ) )^2 )


# sanity check: check first entry, using numbers straight from table
expect_equal( dfr2$yi[1],
              log(0.82 / 0.87) )


v1 = scrape_meta( est = 0.82, hi = .87 )$vyi
v2 = scrape_meta( est = 0.87, hi = .97 )$vyi

expect_equal( dfr2$vyi[1],
              deltamethod( ~ log(x1) - log(x2),
                           mean = c( .82, .87 ), cov = diag( c(v1, v2) ) )^2 )



dfr2$study = "Franklin 2020"




########################### MERGE AND SAVE RESULTS ########################### 

d = rbind( dg3 %>% select( c("study", "yi", "vyi") ),
           ds2 %>% select( c("study", "yi", "vyi") ),
           dfr2 %>% select( c("study", "yi", "vyi") ) )

setwd(prepped.data.dir)
write.csv(d, "prepped_merged_data.csv")



