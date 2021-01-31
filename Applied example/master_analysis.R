
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



# # choose bias SD by looking at a particular cross-section of the plot (muB)
# B_pct( muB = log(2), sigB = 0.01 )
  
digits = 2
varB.ratio = c(0, 0.8)
boot.reps = 2000  # ~~ increase later!

# start results df from scratch
if( exists("resE") ) rm("resE")


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

####################### NORMALITY ASSUMPTION ####################### 

# investigate normality for each meta-analysis
# all seem nicely normal :)
for ( dat in dfs ){
  
  dat$calib = calib_ests( yi = dat$yi,
                          sei = sqrt(dat$vi) )
  
  print( shapiro.test(dat$calib) )
  
  qqnorm(dat$calib)
  qqline(dat$calib, col = "red")
}

# shared across analyses
r = 0.15
take.exp = TRUE
varB.ratio = c(0, 0.8)
q = log(1.1)
boot.reps = 2000

####################### SIMPLE ANALYSES - KODAMA ####################### 

rm( list = ls()[ !ls() %in% c("dfs",
                              "r",
                              "take.exp", 
                              "varB.ratio",
                              "q",
                              "boot.reps") ])

# use local versions of EValue code
detach(package:EValue, unload = TRUE)
setwd("~/Dropbox/Personal computer/Independent studies/EValue package/evalue_package_git/EValue/R")
source("EValue.R")

##### Kodama #####
dat = dfs[[3]]
meta.name = "Kodama"
tail = "above"
muB.phat = log(1.5)


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
sigB = varB.ratio * sqrt(t2)

# to report
exp(mu)
exp(mu.lo)
exp(mu.hi)
sqrt(t2)

##### Naive Phat #####
( Phat.para = MetaUtility::prop_stronger( q = q, 
                                          M = mu,
                                          t2 = t2,
                                          se.M = mu.se,
                                          se.t2 = meta$se.tau2,
                                          tail = tail, 
                                          estimate.method = "parametric",
                                          ci.method = "parametric",
                                          dat = dat, 
                                          R = boot.reps, 
                                          yi.name = "yi",
                                          vi.name = "vi" ) )

( Phat.calib = MetaUtility::prop_stronger( q = q, 
                                         tail = tail, 
                                         estimate.method = "calibrated",
                                         ci.method = "calibrated",
                                         dat = dat, 
                                         R = boot.reps, 
                                         yi.name = "yi",
                                         vi.name = "vi" ) )



##### Choose muhatB.phat and q Such That Phat in [0.15, 0.85] for both homo and hetero bias #####
( PhatB = confounded_meta( q = q,
                          muB = muB.phat,
                          sigB = sigB[2],
                          #sigB = 0,
                          yr = mu,
                          vyr = mu.se^2,
                          t2 = t2,
                          vt2 = meta$se.tau2^2,
                          tail = tail ) )

# here, hetero bias is conservative; compare to homo case:
( PhatB = confounded_meta( q = q,
                           muB = muB.phat,
                           sigB = sigB[1],
                           #sigB = 0,
                           yr = mu,
                           vyr = mu.se^2,
                           t2 = t2,
                           vt2 = meta$se.tau2^2,
                           tail = tail ) )

##### That and Ghat: Parametric #####
( cm = confounded_meta( q = q,
                        r = r,
                        yr = mu,
                        vyr = mu.se^2,
                        t2 = t2,
                        vt2 = meta$se.tau2^2,
                        tail = tail ) )

##### E-Value #####
# E-value: 2.84 to shift to null
( evals = evalues.RR( est = exp(mu),
                      lo = exp(mu.lo),
                      hi = exp(mu.hi),
                      true = 1 ) )

# # sanity check: compare Ghat with r = 0.5 (3.27) to non-null E-value (2.51)
# expq = 1.1
# confounded_meta( q = log(expq),
#                  muB = muB.phat,
#                  sigB = 0,
#                  r = .5,
#                  yr = mu,
#                  vyr = mu.se^2,
#                  t2 = t2,
#                  vt2 = meta$se.tau2^2,
#                  tail = "above" )
# 
# evalues.RR( est = exp(mu),
#             lo = exp(mu.lo),
#             hi = exp(mu.hi),
#             true = expq ) 


####################### SIMPLE ANALYSES - FLEGAL ####################### 

rm( list = ls()[ !ls() %in% c("dfs",
                              "r",
                              "take.exp", 
                              "varB.ratio",
                              "q",
                              "boot.reps") ])

# use local versions of EValue code
detach(package:EValue, unload = TRUE)
setwd("~/Dropbox/Personal computer/Independent studies/EValue package/evalue_package_git/EValue/R")
source("EValue.R")


dat = dfs[[1]]
meta.name = "Flegal"
q = log(.9)
tail = c("below")
Bmax = 1.3
muB.phat = log(1.03)

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

sigB = varB.ratio * sqrt(t2)

# to report
exp(mu)
exp(mu.lo)
exp(mu.hi)
sqrt(t2)

##### Naive Phat #####
( Phat.para = MetaUtility::prop_stronger( q = q, 
                                          M = mu,
                                          t2 = t2,
                                          se.M = mu.se,
                                          se.t2 = meta$se.tau2,
                                          tail = tail, 
                                          estimate.method = "parametric",
                                          ci.method = "parametric",
                                          dat = dat, 
                                          R = boot.reps, 
                                          yi.name = "yi",
                                          vi.name = "vi" ) )

( Phat.calib = MetaUtility::prop_stronger( q = q, 
                                           tail = tail, 
                                           estimate.method = "calibrated",
                                           ci.method = "calibrated",
                                           dat = dat, 
                                           R = boot.reps, 
                                           yi.name = "yi",
                                           vi.name = "vi" ) )

# and effects in the other direction (above 1.1)
# 9% (95% CI: [4%, 14%]) 
( MetaUtility::prop_stronger( q = log(1.1), 
                              M = mu,
                              t2 = t2,
                              se.M = mu.se,
                              se.t2 = meta$se.tau2,
                              tail = "above", 
                              estimate.method = "calibrated",
                              ci.method = "calibrated",
                              dat = dat, 
                              R = boot.reps, 
                              yi.name = "yi",
                              vi.name = "vi" ) )


##### Choose q Such That Phat in [0.15, 0.85] for both homo and hetero bias #####
( PhatB = confounded_meta( q = q,
                           muB = muB.phat,
                           sigB = sigB[2],
                           #sigB = 0,
                           yr = mu,
                           vyr = mu.se^2,
                           t2 = t2,
                           vt2 = meta$se.tau2^2,
                           tail = tail ) )

# here, homo bias is conservative; compare to homo case:
( PhatB = confounded_meta( q = q,
                           muB = muB.phat,
                           sigB = sigB[1],
                           #sigB = 0,
                           yr = mu,
                           vyr = mu.se^2,
                           t2 = t2,
                           vt2 = meta$se.tau2^2,
                           tail = tail ) )

##### That and Ghat: Parametric #####
( cm = confounded_meta( q = q,
                        r = r,
                        yr = mu,
                        vyr = mu.se^2,
                        t2 = t2,
                        vt2 = meta$se.tau2^2,
                        tail = tail ) )
# calibrated version is pretty similar :)

##### E-Value #####
# E-value: 1.35 to shift to null
( evals = evalues.RR( est = exp(mu),
                      lo = exp(mu.lo),
                      hi = exp(mu.hi),
                      true = 1 ) )


####################### SIMPLE ANALYSES - GBC ####################### 


rm( list = ls()[ !ls() %in% c("dfs",
                              "r",
                              "take.exp", 
                              "varB.ratio",
                              "q",
                              "boot.reps") ])

# use local versions of EValue code
detach(package:EValue, unload = TRUE)
setwd("~/Dropbox/Personal computer/Independent studies/EValue package/evalue_package_git/EValue/R")
source("EValue.R")

dat = dfs[[2]]
meta.name = "GBC"
tail = c( "above" )
muB.phat = log(1.03)
q = log(1.1)



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

sigB = varB.ratio * sqrt(t2)

# to report
exp(mu)
exp(mu.lo)
exp(mu.hi)
sqrt(t2)

##### Naive Phat #####
( Phat.para = MetaUtility::prop_stronger( q = q, 
                                          M = mu,
                                          t2 = t2,
                                          se.M = mu.se,
                                          se.t2 = meta$se.tau2,
                                          tail = tail, 
                                          estimate.method = "parametric",
                                          ci.method = "parametric",
                                          dat = dat, 
                                          R = boot.reps, 
                                          yi.name = "yi",
                                          vi.name = "vi" ) )

( Phat.calib = MetaUtility::prop_stronger( q = q, 
                                           tail = tail, 
                                           estimate.method = "calibrated",
                                           ci.method = "calibrated",
                                           dat = dat, 
                                           R = boot.reps, 
                                           yi.name = "yi",
                                           vi.name = "vi" ) )

# proportion in the other direction: 2% [0%, 8%]
MetaUtility::prop_stronger( q = log(.9), 
                            M = mu,
                            t2 = t2,
                            se.M = mu.se,
                            se.t2 = meta$se.tau2,
                            tail = "below", 
                            estimate.method = "calibrated",
                            ci.method = "calibrated",
                            dat = dat, 
                            R = boot.reps, 
                            yi.name = "yi",
                            vi.name = "vi" )



##### Choose q Such That Phat in [0.15, 0.85] for both homo and hetero bias #####
( PhatB = confounded_meta( q = q,
                           muB = muB.phat,
                           sigB = sigB[2],
                           #sigB = 0,
                           yr = mu,
                           vyr = mu.se^2,
                           t2 = t2,
                           vt2 = meta$se.tau2^2,
                           tail = tail ) )

# here, homo bias is conservative; compare to homo case:
( PhatB = confounded_meta( q = q,
                           muB = muB.phat,
                           sigB = sigB[1],
                           #sigB = 0,
                           yr = mu,
                           vyr = mu.se^2,
                           t2 = t2,
                           vt2 = meta$se.tau2^2,
                           tail = tail ) )

##### That and Ghat: Parametric #####
( cm = confounded_meta( q = .q,
                        r = r,
                        yr = mu,
                        vyr = mu.se^2,
                        t2 = t2,
                        vt2 = meta$se.tau2^2,
                        tail = tail ) )

##### E-Value #####
# E-value: 1.35 to shift to null
( evals = evalues.RR( est = exp(mu),
                      lo = exp(mu.lo),
                      hi = exp(mu.hi),
                      true = 1 ) )




####################### COMPLICATED ANALYSES ####################### 


# maybe use same Bmax throughout so that the plots are more clear?

##### Flegal #####
dat = dfs[[1]]
meta.name = "Flegal"
ql = list( log(1), log(.9) )
tail = c( "above", "below" )
take.exp = TRUE
r = 0.10
varB.ratio = c(0, 0.8)
Bmax = 1.3
breaks.x1 = seq(1, Bmax, .1)
est.type = "HR"

setwd(code.dir)
source("general_analysis.R")

##### GBC #####

dat = dfs[[2]]
meta.name = "GBC"
ql = list( log(1), log(1.1) )
tail = c( "below", "above" )
take.exp = TRUE
r = 0.10
varB.ratio = c(0, 0.8)
Bmax = 1.5
breaks.x1 = seq(1, Bmax, .1)
est.type = "HR"

setwd(code.dir)
source("general_analysis.R")


##### Kodama #####
dat = dfs[[3]]
meta.name = "Kodama"
ql = list( log(1),
           #log(1.5),
           log(1.1) )  # different from others because otherwise Phat is too close to 1
tail = c( "below", "above", "above" )
take.exp = TRUE
r = 0.10
varB.ratio = c(0, 0.8)
Bmax = 3
breaks.x1 = seq(1, 3, .5)
est.type = "RR"

setwd(code.dir)
source("general_analysis.R")



# # with a different choice of Bmax and 
# dat = dfs[[3]]
# meta.name = "Kodama"
# ql = list( log(1.1) )  # different from others because otherwise Phat is too close to 1
# tail = c( "above" )
# take.exp = TRUE
# r = 0.10
# varB.ratio = c(0, 0.8)
# Bmax = 3
# breaks.x1 = seq(1, 3, .5)
# est.type = "RR"
# 
# setwd(code.dir)
# source("general_analysis.R")
# 
# # choose a good threshold for Kodama
# Phat_causal( .q = log(1.1),
#           .B = 2,
#           .calib = calib_ests( dfs[[3]]$yi, sqrt(dfs[[3]]$vi) ), # assumed on log scale
#           .tail = "above",
#           .causative = TRUE,
#           
#           .give.CI = FALSE,
#           .R = 500,
#           .dat = dfs[[3]] )


##### Look At Aggregated Results #####
resE

setwd(results.dir)
write.csv(resE, "all_results.csv")
  