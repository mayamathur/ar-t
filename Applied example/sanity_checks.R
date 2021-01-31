
rm( list = ls() )

####################### SET UP ####################### 

code.dir = "~/Dropbox/Personal computer/Independent studies/2020/Metasens tutorial/Linked to OSF (Metasens tutorial)/Code (git)"
private.data.dir = "~/Dropbox/Personal computer/Independent studies/2020/Metasens tutorial/Data (private)"
results.dir = "~/Dropbox/Personal computer/Independent studies/2020/Metasens tutorial/Linked to OSF (Metasens tutorial)/Results from R"
overleaf.dir = "~/Dropbox/Apps/Overleaf/Metasens tutorial"

setwd(code.dir)
source("helper.R")

setwd(private.data.dir)
dfs = list( read.csv("flegal_prepped.csv"),
            read.csv("gbc_prepped.csv"),
            read.csv("kodama_prepped.csv") )

library(MetaUtility)
library(ggplot2)
library(dplyr)
library(metafor)
library(testthat)
library(boot)

# use local versions of EValue code
detach(package:EValue, unload = TRUE)
setwd("~/Dropbox/Personal computer/Independent studies/EValue package/evalue_package_git/EValue/R")
source("EValue.R")


####################### SPOT-CHECK KODAMA #######################

##### Set Up #####
dat = dfs[[3]]
meta.name = "Kodama"
ql = list( log(1), log(1.5), log(1.1) )  # different from others because otherwise Phat is too close to 1
tail = c( "below", "above", "above" )
take.exp = TRUE
r = 0.10
varB.ratio = c(0, 0.8)
Bmax = 3
breaks.x1 = seq(1, 3, .5)
est.type = "RR"

setwd(code.dir)
source("general_analysis.R")

##### Naive Meta-Analysis #####
( meta = rma.uni( yi = dat$yi, 
                  vi = dat$vi, 
                  method = "REML", 
                  knha = TRUE ) )

mu = meta$b
t2 = meta$tau2
mu.lo = meta$ci.lb
mu.hi = meta$ci.ub
mu.se = meta$se
mu.pval = meta$pval

##### Choose q Such That Phat in [0.15, 0.85] for both homo and hetero bias #####
.q = log(1.1)
( Phat = confounded_meta( q = .q,
                        muB = log(1.5),
                        sigB = sqrt(t2) * .8,
                        #sigB = 0,
                        yr = mu,
                        vyr = mu.se^2,
                        t2 = t2,
                        vt2 = meta$se.tau2^2,
                        tail = "above" ) )

##### That and Ghat: Parametric #####
r = .15
( cm = confounded_meta( q = .q,
                 r = r,
                 yr = mu,
                 vyr = mu.se^2,
                 t2 = t2,
                 vt2 = meta$se.tau2^2,
                 tail = "above" ) )
# calibrated version is pretty similar :)


##### E-Value #####
# E-value: 2.51 to shift to 1.1
( evals = evalues.RR( est = exp(mu),
                    lo = exp(mu.lo),
                    hi = exp(mu.hi),
                    true = exp(.q) ) )


# expect this to be close to Ghat(q = mu, r = 0.50) (i.e., shift the mean to q)
confounded_meta( q = log(.q),
                 r = 0.50,
                 muB = log(mu),
                 sigB = 0,
                 yr = mu,
                 vyr = mu.se^2,
                 t2 = t2,
                 vt2 = meta$se.tau2^2,
                 tail = "above" )

# I think the huge discrepancy reflects the impact of allowing for heterogeneous bias
#  in this one: compare the two plots, for example

##### 
TG = That_causal( .q = .q,  # on log-RR scale (transformed if needed)  #~~TEMP ONLY
                  .r = r,
                  .B.vec = seq( 1, 5, .001 ),  # RR scale
                  .calib = dat$calib, # on log-RR scale (transformed if needed)
                  .tail = "above",
                  .causative = (mu>0),
                  
                  .give.CI = TRUE,
                  .R = 500,
                  .dat = dat,
                  .calib.name = "calib" )
# **Ghat: 3.07










evalues.RR()