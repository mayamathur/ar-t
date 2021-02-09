
# Goal: Try to break the conservative interpreation of NPPhat
#  i.e., assume homogeneous Bhomo and calculate Phat
#  then treat that as a lower *bound* on the real Phat
#  if in fact the bias were heterogeneous but <= Bhomo

# aiming to get a situation in which the calibrated estimates
#  are not sufficiently shrunk because of the heterogeneous bias
#  st some calibrated estimates remain > q

library(MetaUtility)
library(testthat)
library(metafor)


############################## SIMULATE DATA ############################## 

k = 500

# true causal mean is less than q with no heterogeneity
# so true Phat = 0
muT = 0.8
q = 1.2
VT = 0.3
# study-level population causal effects
thT = rnorm( mean = muT, sd = sqrt(VT), n = k )
summary(thT)

# assumed homogeneous bias 
Bhomo = 0.4

# true bias in each study (very heterogeneous)
muB.star = 0.2
Bi.star = rnorm( mean = muB.star, sd = 5, n = k ) 
# truncate it
Bi.star = pmin( Bhomo, Bi.star )
# optionally, make it always > 0 (biased upward)
Bi.star = pmax( 0, Bi.star)
summary(Bi.star)

# draw studies
se = runif( n = k, min = 0.1, max = 0.5 )
#se = rep(0.0001,k)  # n -> infinity for all studiess
eps = rnorm( n = k, mean = 0, sd = se )
# population confounded effects with hetero bias
thC = thT + Bi.star

# study-level confounded estimates
yiC = thC + eps


# TARGET: population P>q for causal effects
( truth = mean( thT > q ) )

# sample target for Phat if we had bias distribution correct: population P>q for confounded effects
mean( thC > q )

# look at distributions
hist(thT)
hist(thC)
hist(yiC)


############################## TRUE PHATS ############################## 

# sample estimate with WRONG (homo) bias distribution
# using naive t2 as in NPPhat, which is too large
# see the naive t2
meta1 = rma.uni( yi = yiC, sei = se, method = "DL" )
meta1$b; mean(thC)  # should be close
sqrt(meta1$tau2); sd(Bi.star)  # should be pretty big because of bias


##### Hetero Calibration ######
# benchmark: sample estimate if we had the CORRECT bias distribution
# get the corrected t2 FIRST
yiT = yiC - Bi.star
meta2 = rma.uni( yi = yiT, sei = se, method = "DL" )
meta2$tau2  # should be small

if (meta2$tau2 == 0){
  ens2 = rep(meta2$b, k)
  } else {
    # bias-correct the meta-analysis itself, then calibrate
    # note use of yiT for the residual
    ens2 = c(meta2$b) + sqrt( c(meta2$tau2) / ( c(meta2$tau2) + se^2 ) ) * ( yiT - c(meta2$b) )
  } 


# c.f. another (equivalent?) way that's more useful for comparing to homo calibration:
# use confounded point estimate but correct heterogeneity when calibrating
# then bias-correct each calibrated estimate at the end
# this is 
meta1$b - mean(Bi.star); meta2$b  # should be similar

# **HETERO CALIBRATED ESTIMATES
ensHet = ( c(meta1$b) - mean(Bi.star) + sqrt( c(meta2$tau2) / ( c(meta2$tau2) + se^2 ) ) * ( yiC - Bi.star - c(meta1$b) + mean(Bi.star) ) ) 

hist(ensHet)
hist(ens2)
expect_equal( ens2, ensHet, tol = 0.001 )  # yes, seems equivalent
# *critical: we need to subtract Bi.star WITHIN the residuals, not at the end
#  otherwise the residuals are overdispersed

( PhatHetero = mean( ensHet > q ) ); truth


##### Homo Calibration #####

#ens4 = calib_ests( yi = yiC, sei = se ) - Bhomo

# **HOMO CALIBRATED ESTIMATES
ensHomo = ( c(meta1$b) - Bhomo + sqrt( c(meta1$tau2) / ( c(meta1$tau2) + se^2 ) ) * ( yiC - Bhomo - c(meta1$b) + Bhomo ) ) 

hist(ensHomo)

( PhatHomo = mean(ensHomo > q) )


##### 2021-2-8 Calibration #2 #####

# using mean bias,  but wrong residual
ensHomo2 = ( c(meta1$b) - mean(Bi.star) + sqrt( c(meta2$tau2) / ( c(meta2$tau2) + se^2 ) ) * ( yiC - c(meta1$b) ) ) 

( PhatHomo2 = mean(ensHomo2 > q) )
# **yes, unfortunately this is too large because the residuals are too variable


##### 2021-2-8 Calibration #2 #####

# try to fix the variance of residuals
# changed denom of calibration term to the CONFOUNDED heterogeneity to cancel out when taking
#  var of this whole thing
ensHomo3 = ( c(meta1$b) - mean(Bi.star) + sqrt( c(meta2$tau2) / ( c(meta1$tau2) + se^2 ) ) * ( yiC - c(meta1$b) ) ) 

# should be close: 
sqrt( c(meta2$tau2) / ( c(meta1$tau2) + se^2 ) ) * ( yiC - c(meta1$b) ); VT

( PhatHomo3 = mean(ensHomo3 > q) )


##### Toward Theory #####

# look at the yi that are candidates to be above vs. below q with different calibration #####

mean( yiC < c(meta1$b) - mean(Bi.star) )  # these will increase when calibrating but remain <q regardless of how we calibrate (because the *corrected* mean estimate is < q)
mean( yiC > c(meta1$b) - mean(Bi.star) & yiC < q )  # this will decrease when calibrating but again remain <q regardless of how we calibrate
mean( yiC > q )  # these are the important ones that could become <q depending on how we calibrate



# of these candidates, how often are they >q under homo calibration but <q under hetero calibration?
ind = which( yiC > q )
mean( yiT[ind] > q )
mean( ensHomo[ind] > q )
mean( ensHet[ind] > q )

# are they EVER larger with homo calibration?
mean( ensHomo[ind] > ensHet[ind] )
# **yes, sometimes


# check intermediate theory
c1 = sqrt( c(meta1$tau2) / ( c(meta1$tau2) + se^2 ) )
c2 = sqrt( c(meta2$tau2) / ( c(meta2$tau2) + se^2 ) )

expect_equal( ensHomo - ensHet,
              ( muB.star - Bhomo ) + (c1 - c2)*yiC - (c1 - c2)*c(meta1$b) + c2*Bi.star - c2*c(muB.star),
              tol = 0.001 )

summary( (c2 - c1) / (c2-1) )


