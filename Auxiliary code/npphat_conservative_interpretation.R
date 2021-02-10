
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


# generate from a uniform mixture with endpoints [-b, -a] and [a, b]
#  but shifted so that the grand mean is mu
runif2 = function(n,
                  mu,
                  V) {
  # calculate lower limit for positive distribution, a
  # arbitrary, but can't have too large or else it;s impossible to find
  #  a valid b
  a = sqrt(V)/2  
  
  # calculate upper endpoint for positive distribution, b
  b = abs( 0.5 * ( sqrt(3) * sqrt( 4*V - a^2 ) - a ) )
  
  # prior to mean shift
  components = sample(1:2,
                      prob=c(0.5, 0.5),
                      size=n,
                      replace=TRUE)
  
  mins = c( -b, a )
  maxes = c( -a, b )
  
  samples = runif(n=n,
                  min=mins[components],
                  max=maxes[components])
  
  # mean-shift them
  samples = samples + mu
  
  return( list(x = samples,
               a = a, 
               b = b) )
}


############################## SIMULATE DATA ############################## 

k = 500

# need muT < q to be in the correct case
muT = 1
q = 1.2
VT = 0.25
# study-level population causal effects
thT = rnorm( mean = muT, sd = sqrt(VT), n = k )
summary(thT)

# true bias in each study (very heterogeneous)
muB.star = 0.2

# assumed homogeneous bias 
Bhomo = muB.star


##### Draw Bias #####
# # truncated bias
# Bi.star = rnorm( mean = muB.star, sd = 5, n = k ) 
# # truncate it
# Bi.star = pmin( Bhomo, Bi.star )
# # optionally, make it always > 0 (biased upward)
# Bi.star = pmax( 0, Bi.star)
# summary(Bi.star)

# # non-truncated exponential bias
# Bi.star = rexp( n = k, rate = 1/muB.star )
# # hist(Bi.star)
# summary(Bi.star)

# **SAVE: non-truncated bimodal bias
# this is a great example because ensHomo is really biased, but ensHomo3 is not 
# bias heterogeneity is in terms of true heterogeneity 
Bi.star = runif2( n = k, mu = muB.star, V = VT * 1.2)$x
hist(Bi.star)

# bias correlated with true effects and exponential
# shifted to have desired muB
beta = 1.3
x = rexp( n = k, rate = 1/(thT*beta) )
Bi.star = x - mean(x) + muB.star
hist(Bi.star)
summary(Bi.star)
mean(Bi.star); muB.star


##### Draw Studies #####
se = runif( n = k, min = 0.1, max = .1 )
#se = rep(0.0001,k)  # n -> infinity for all studiess
eps = rnorm( n = k, mean = 0, sd = se )
# population confounded effects with hetero bias
thC = thT + Bi.star
hist(thC)

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
# **yes, unfortunately this is too large because the residuals are too variable



##### ***2021-2-8 Calibration #2 #####

# **this one works, but obviously the numerator is unknown, so in practice we have to bound it
# try to fix the variance of residuals
# changed denom of calibration term to the CONFOUNDED heterogeneity to cancel out when taking
#  var of this whole thing
ensHomo3 = ( c(meta1$b) - mean(Bi.star) + sqrt( c(meta2$tau2) / ( c(meta1$tau2) + se^2 ) ) * ( yiC - c(meta1$b) ) ) 

# should be close: 
sqrt( c(meta2$tau2) / ( c(meta1$tau2) + se^2 ) ) * ( yiC - c(meta1$b) ); VT

( PhatHomo3 = mean(ensHomo3 > q) )


# this really seems to work pretty well


##### ***2021-2-9 Calibration #4 #####
# now use our upper bound on VB in the numerator
# suppose we actually guess the upper bound correctly

# sanity check
expect_equal( var(thT),
              var(thC) - var(Bi.star) - 2 * cor(thT, Bi.star) * sd(thT) * sd(Bi.star) )

# upper bound setting correlation to 1
UB = var(Bi.star)
UB.rho = 0.5
var(thT); var(thC) - UB - 2 *UB.rho* sd(thT) * sqrt(UB)

# conservative heterogeneity estimate
t2t.cons = min(0, var(thC) - UB - 2 *UB.rho* sd(thT) * sqrt(UB) )

ensHomo4 = ( c(meta1$b) - mean(Bi.star) + sqrt( c( ) / ( c(meta1$tau2) + se^2 ) ) * ( yiC - c(meta1$b) ) ) 

# should be close: 
sqrt( c(meta2$tau2) / ( c(meta1$tau2) + se^2 ) ) * ( yiC - c(meta1$b) ); VT

( PhatHomo3 = mean(ensHomo3 > q) )


# this really seems to work pretty well





