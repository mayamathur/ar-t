

rm( list = ls() )

####################### SET UP ####################### 

library(MetaUtility)
library(ggplot2)
library(dplyr)
library(metafor)
library(testthat)
library(boot)
#library(EValue)
library(here)
library(data.table)

#@use dev versions of code that have heterogeneous Tmin, Gmin
detach("package:EValue", unload = TRUE)
setwd("~/Dropbox/Personal computer/Independent studies/R packages/EValue package (git)/evalue_package/EValue/R")
source("meta-analysis.R")

code.dir = here("Applied example/Code")
data.dir = here("Applied example/data")
results.dir = here("Applied example/Results from R")
overleaf.dir = "~/Dropbox/Apps/Overleaf/AR-T (Ann Rev tutorial on bias in meta-analyses)/R_objects"

setwd(code.dir)
source("helper.R")

setwd(data.dir)
d = fread("kodama_prepped.csv")

# rounding digits
digits = 2

# no sci notation
options( scipen = 99 )

# should we redo the slow, bootstrapped plot?
redo.bootstrapped.plot = FALSE

####################### NAIVE META-ANALYSIS ####################### 

# got through each meta-analysis and fit naive analysis
# also check normality, relevant only for parametric sensitivity analyses that use
#   heterogeneous bias

##### Naive meta-analysis #####
( meta = rma.uni( yi = d$yi, 
                  vi = d$vi, 
                  slab = d$study,  # for use in forest plot later
                  measure = "RR",  # ditto
                  method = "REML", 
                  knha = TRUE ) )

mu = meta$b
t2 = meta$tau2
mu.lo = meta$ci.lb
mu.hi = meta$ci.ub
mu.se = meta$se
mu.pval = meta$pval

print(meta)

# to report
cat("\n***Naive estimate (RR or HR) and CI:", paste( round( exp(mu), digits ),
                                                     " [",
                                                     round( exp(mu.lo), digits ),
                                                     ", ",
                                                     round( exp(mu.hi), digits ),
                                                     "]",
                                                     sep = "") )  

update_result_csv( name = "Naive meta est",
                   value = round( exp(mu), digits ) )

update_result_csv( name = "Naive meta tau",
                   value = round( sqrt(t2), digits ) )


update_result_csv( name = "Naive meta lo",
                   value = round( exp(mu.lo), digits ) )


update_result_csv( name = "Naive meta hi",
                   value = round( exp(mu.hi), digits ) )


update_result_csv( name = "Naive meta pval",
                   value = format.pval( mu.pval, eps=0.001 ) )


##### Normality #####
d$calib = calib_ests( yi = d$yi,
                      sei = sqrt(d$vi) )

cat("\n\nShapiro test: ")
print( shapiro.test(d$calib) )

qqnorm(d$calib)
qqline(d$calib, col = "red")


####################### SENSITIVITY ANALYSES ####################### 

###### E-value for point estimate ######
evals = evalue( est = RR( exp(meta$b) ),
                lo = RR( exp(meta$ci.lb) ),
                hi = RR( exp(meta$ci.ub) ) )

# save results
update_result_csv( name = "Evalue est",
                   value = round( evals[ "E-values", "point" ], digits ) )

update_result_csv( name = "Evalue CI",
                   value = round( evals[ "E-values", "lower" ], digits ) )


###### Naive Phat ######
# uncorrected estimate of effects with RR > 1.1
# Tmin and Gmin here also give amount of bias and of confounding required
#   to reduce to less than 15% the percentage of meaningfully large effects
cm = confounded_meta( method = "calibrated",
                      q = log(1.1),
                      r = 0.15,
                      tail = "above",
                      muB = 0,
                      sigB = 0,
                      d = d,
                      yi.name = "yi",
                      vi.name = "vi" )

# save results
update_result_csv( name = "Naive perc gt 1.1",
                   value = round( 100*cm$Est[cm$Value == "Prop"], 0 ) )

# no CI because at the ceiling


###### Homogeneous Bias ######
cmHomo = confounded_meta( method = "calibrated",
                      q = log(1.1),
                      r = 0.15,
                      tail = "above",
                      
                      dat = d,
                      yi.name = "yi",
                      vi.name = "vi",
                      R = 1000 )

# save for plot
ThatHomo = cmHomo$Est[cmHomo$Value == "Tmin"]
GhatHomo = cmHomo$Est[cmHomo$Value == "Gmin"]

# save results
update_result_csv( name = "Ghat homogeneous",
                   value = round( GhatHomo, 2 ) )

update_result_csv( name = "Ghat homogeneous lo",
                   value = round( cmHomo$CI.lo[cmHomo$Value == "Gmin"], 2 ) )

update_result_csv( name = "Ghat homogeneous hi",
                   value = round( cmHomo$CI.hi[cmHomo$Value == "Gmin"], 2 ) )


###### Heterogeneous Bias ######

cmHetero = confounded_meta( method = "parametric",
                      q = log(1.1),
                      r = 0.15,
                      tail = "above",
                      muB = log(1.5),
                      sigB = sqrt( 0.8 * meta$tau2 ),
                      yr = meta$b,
                      vyr = meta$se^2,
                      t2 = meta$tau2,
                      vt2 = meta$se.tau2^2 )


# save for plot
#@c.f. homogeneous bias with parametric estimation: Tmin = 1.93, Gmin = 3.27
ThatHetero = cmHetero$Est[cmHetero$Value == "Tmin"]
GhatHetero = cmHetero$Est[cmHetero$Value == "Gmin"]

# save results
update_result_csv( name = "Ghat hetero",
                   value = round( GhatHetero, 2 ) )

update_result_csv( name = "Ghat hetero lo",
                   value = round( cmHetero$CI.lo[cmHetero$Value == "Gmin"], 2 ) )

update_result_csv( name = "Ghat hetero hi",
                   value = round( cmHetero$CI.hi[cmHetero$Value == "Gmin"], 2 ) )


####################### PLOTS #######################

##### Forest Plot #####

for ( dir in c(results.dir, overleaf.dir) ) {
  setwd(dir)
  png( filename = "kodama_forest.png" )
  forest(meta,
         showweights = FALSE,
         header = TRUE,
         xlab = "Estimate (RR)",
         mlab = "Pooled estimate",
         transf = exp,
         steps = 9,
         order = "obs",
         col = "orange",
         refline = 1 )
  dev.off()
}

##### Homogeneous Bias Line Plot #####

if ( redo.bootstrapped.plot == TRUE ) {
  p = sens_plot( method = "calibrated",
                 type = "line",
                 
                 q = log(1.1),
                 tail = "above",
                 Bmax = log(2.5),
                 breaks.x1 = seq( 1, 2.5, .25 ),
                 
                 dat = d,
                 yi.name = "yi",
                 vi.name = "vi",
                 R = 1000 )
  
  # customize
  p = p + geom_hline( yintercept = 0.15, lty = 2, color = "red" ) +
    geom_vline( xintercept = ThatHomo, lty = 2, color = "red" )
  
}

my_ggsave( name = "kodama_plot_homo.pdf",
           width = 6,
           height = 5 )


##### Heterogeneous Bias Line Plot #####

sens_plot( method = "parametric",
           type = "line",
           q = log(1.1),
           tail = "above",
           Bmax = log(2.5),
           breaks.x1 = seq( 1, 2.5, .25 ),
           
           muB = log(1.5),
           sigB = sqrt( 0.8 * meta$tau2 ),
           yr = meta$b,
           vyr = meta$se^2,
           t2 = meta$tau2,
           vt2 = meta$se.tau2^2 )

p2 = last_plot()

# customize
p2 = p2 + geom_hline( yintercept = 0.15, lty = 2, color = "red" ) +
  geom_vline( xintercept = ThatHetero, lty = 2, color = "red" )

my_ggsave( name = "kodama_plot_hetero.pdf",
           width = 6,
           height = 5 )

















