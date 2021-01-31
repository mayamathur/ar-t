

####################### ** CV FITNESS & MORTALITY ####################### 

# using the ones with adjusted RRs; would be interesting to compare to the 
# crude RRs! 
rr = c( 1.47,
       1.53,
       1.25,
       1.22,
       1.92,
       1.25,
       2.23,
       1.38,
       1.42,
       2.58,
       1.83,
       2.38,
       2.09,
       2.17,
       1.89,
       2.96 )
hi = c( 1.55,
        1.85,
        1.47,
        2.57,
        2.59,
        1.8,
        3.3, 
        2.83,
        1.87,
        4.09,
        2.73,
        4.14,
        2.96,
        4.30,
        2.86,
        4.46 )

# already have RRs
library(MetaUtility)
d = MetaUtility::scrape_meta( type = "RR",
                              est = rr, 
                              hi = hi,
                              sqrt = FALSE )

library(metafor)
( meta = rma.uni( yi = d$yi,
                  vi = d$vyi,
                  method = "REML",
                  knha = TRUE ) )

exp(meta$b)
# similar to reported (1.70), but not exact


q = log(1.1)

# for Shapiro test (seems okay :) )
MetaUtility::prop_stronger(q = q,
                           M = meta$b,
                           t2 = meta$tau2,
                           se.t2 = meta$se.tau2,
                           tail = "above",
                           dat = d, 
                           yi.name = "yi",
                           vi.name = "vyi" )


sens_plot( type = "line",
           q = q, 
           yr = meta$b,
           vyr = meta$se^2,
           t2 = meta$tau2,
           vt2 = meta$se.tau2^2,
           Bmax = log(2.5) )

# y_R^t = meta$b - log(1.5)
# here q < y_R^t
# y_R^c > 0
# so Phat under homogeneity is a *lower* bound and it will only look more robust
#  with heterogeneity
confounded_meta( 
  q = q, 
  r = 0.1,
  muB = log(1.5),
  sigB = 0,
  yr = meta$b,
  vyr = meta$se^2,
  t2 = meta$tau2,
  vt2 = meta$se.tau2^2,
  tail = "above" )
# this bias factor of 1.5 corresponds to RR strengths of 2.4
# Ghat = 3.46


# E-value for point estimate: 3.06 and 2.27
evalues.RR(est = meta$b, hi = meta$ci.ub)



####################### SMOKING & DIABETES ####################### 

rr = c(1.50,
       1.2,
       1.88,
       2.38,
       .82,
       1.42,
       1.47,
       1.62,
       2.74,
       1.63,
       1.13,
       1.74,
       1.3,
       2.47,
       1.24,
       1.35,
       1.06,
       3.74,
       1.5,
       1.46,
       1.31,
       1.94,
       2.15,
       1.62,
       1.65)
hi = c(2.1,
       1.8,
       3.02,
       7.40,
       1.3,
       1.83,
       1.92,
       2.59,
       7.13,
       1.77,
       1.19,
       2.43,
       1.47,
       9.3,
       1.66,
       1.54,
       1.3,
       12.91,
       2.1,
       1.66,
       1.62,
       3.25,
       3.11,
       2.5,
       2.13)

# already have RRs
library(MetaUtility)
d = MetaUtility::scrape_meta( type = "RR",
                              est = rr, 
                              hi = hi,
                              sqrt = FALSE )

library(metafor)
( meta = rma.uni( yi = d$yi,
                  vi = d$vyi,
                  method = "REML",
                  knha = TRUE ) )

exp(meta$b)
# similar to reported, but not exact


q = log(1.1)

# for Shapiro test (seems okay :) )
MetaUtility::prop_stronger(q = q,
                           M = meta$b,
                           t2 = meta$tau2,
                           se.t2 = meta$se.tau2,
                           tail = "above",
                           dat = d, 
                           yi.name = "yi",
                           vi.name = "vyi" )


sens_plot( type = "line",
           q = q, 
           yr = meta$b,
           vyr = meta$se^2,
           t2 = meta$tau2,
           vt2 = meta$se.tau2^2,
           Bmax = log(2) )

# y_R^t = meta$b - log(1.25)
# here q < y_R^t
# y_R^c > 0
# so Phat under homogeneity is a *lower* bound and it will only look more robust
#  with heterogeneity
confounded_meta( 
  q = q, 
  r = 0.1,
  muB = log(1.25),
  sigB = 0,
  yr = meta$b,
  vyr = meta$se^2,
  t2 = meta$tau2,
  vt2 = meta$se.tau2^2,
  tail = "above" )
# this bias factor of 1.25 corresponds to RR strengths of 1.81
# Ghat = 2.44


# E-value for point estimate
evalues.RR(est = exp(meta$b), hi = exp(meta$ci.ub))
# 2.19

# I think this example could work, although it's not *super* robust




####################### SEX ABUSE & PELVIC PAIN ####################### 
# unclear whether should be counted as a "rare" outcome; e.g.:
# https://www.ncbi.nlm.nih.gov/pubmed/24658485
or = c(3.33,
       8,
       7.33,
       1.25,
       2.91,
       2.06,
       1.68,
       5.38,
       1.52,
       3.22)
hi = c(9.04,
       52.99,
       38.88,
       2.90,
       8.37,
       4.17,
       6.65,
       12.76,
       4.10,
       7.36)

# will conservatively treat as rare

library(MetaUtility)
d = MetaUtility::scrape_meta( type = "RR",
             est = or, 
             hi = hi,
             sqrt = FALSE )

library(metafor)
( meta = rma.uni( yi = d$yi,
                vi = d$vyi,
                method = "REML",
                knha = TRUE ) )

exp(meta$b)
# similar to reported, but not exact


q = log(1.1)

# for Shapiro test (seems okay :) )
MetaUtility::prop_stronger(q = q,
                           M = meta$b,
                           t2 = meta$tau2,
                           se.t2 = meta$se.tau2,
                           tail = "above",
                           dat = d, 
                           yi.name = "yi",
                           vi.name = "vyi" )


sens_plot( type = "line",
          q = q, 
          yr = meta$b,
          vyr = meta$se^2,
          t2 = meta$tau2,
          vt2 = meta$se.tau2^2)

# y_R^t = meta$b - log(2)
# here q < y_R^t
# y_R^c > 0
# so Phat under homogeneity is a *lower* bound and it will only look more robust
#  with heterogeneity
confounded_meta( 
           q = q, 
           muB = log(2),
           sigB = 0,
           yr = meta$b,
           vyr = meta$se^2,
           t2 = meta$tau2,
           vt2 = meta$se.tau2^2,
           tail = "above" )
# with bias factor of 2 in each study, 77% above RR = 1.1 but CI is not informative at all 

# CONCLUSION: TOO FEW STUDIES FOR INFORMATIVE CI HERE



