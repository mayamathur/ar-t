

if(exists("sigB")) rm(sigB)

##### Regular Meta-Analysis #####
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

# add calibrated estimates to dataset
dat$calib = calib_ests( yi = dat$yi,
                        sei = sqrt(dat$vi) )

##### Check for Bad Arguments Regarding SigB #####

# should not provide both of these

if ( exists("sigB") ) {
  if ( any( !is.na(sigB) ) & any( !is.na(varB.ratio) ) ) stop("Only provide one of sigB and varB.ratio, not both")
  
  if ( any( !is.na(sigB) ) & any(sigB > t2) ) stop( paste( "Some sigB > naive t2 of",
                                                                            round(t2,2),
                                                                            sep = " ") )
}



# calculate absolute sigB if given ratio 
if ( any( !is.na(varB.ratio) ) ) sigB = sqrt(varB.ratio * t2)

print( paste( meta.name, ": using sigB = ", sigB ) )

# browser()
# # sanity check:
# confounded_meta( q = log(1.5),
#                  r = 0.10,
#                  muB = log(1.2),
#                  sigB = 0,
#                  yr = mu,
#                  vyr = mu.se^2,
#                  t2 = t2,
#                  vt2 = meta$se.tau2^2,
#                  tail = "above" )


##### Forest Plot #####

# forest(meta,
#        slab = dat$study,
#        xlab = "Point estimate")

##### Phat Naive - Both Parametric and Calibrated #####
# Phat from calibrated estimates
Phat.l = lapply( ql,
                 FUN = function(.q) {
                   
                   # get the appropriate entry of tail argument
                   ind = which( unlist(ql) == .q )
                   
                   Phat.calib = MetaUtility::prop_stronger( q = .q, 
                                                            tail = tail[ind], 
                                                            estimate.method = "calibrated",
                                                            ci.method = "calibrated",
                                                            dat = dat, 
                                                            R = boot.reps, 
                                                            yi.name = "yi",
                                                            vi.name = "vi" )
                   
                   Phat.para = MetaUtility::prop_stronger( q = .q, 
                                                           M = mu,
                                                           t2 = t2,
                                                           se.M = mu.se,
                                                           se.t2 = meta$se.tau2,
                                                           tail = tail[ind], 
                                                           estimate.method = "parametric",
                                                           ci.method = "parametric",
                                                           dat = dat, 
                                                           R = boot.reps, 
                                                           yi.name = "yi",
                                                           vi.name = "vi" )
                   
                   return( data.frame( Phat.calib = Phat.calib$est,
                                       Phat.calib.lo = Phat.calib$lo,
                                       Phat.calib.hi = Phat.calib$hi,
                                       Phat.para = Phat.para$est,
                                       Phat.para.lo = Phat.para$lo,
                                       Phat.para.hi = Phat.para$hi ) )
                 } )  # end lapply


Phat.df = do.call( rbind, 
                   Phat.l )

# note that this always reports only the calibrated estimates
Phat.df$string = paste( round( 100*Phat.df$Phat.calib,
                               digits = 0 ),
                        format_CI( 100*Phat.df$Phat.calib.lo,
                                   100*Phat.df$Phat.calib.hi,
                                   digits = 0 ),
                        sep = " " )


##### That and Ghat - Calibrated #####

TG.l = lapply( X = ql,
               FUN = function(.q){
                 
                 # index of the q that we're working on
                 ind = which( unlist(ql) == .q )
                 
                 # get both That and Ghat
                 TG = That_causal( .q = .q,  # on log-RR scale (transformed if needed)  #~~TEMP ONLY
                                   .r = r,
                                   .B.vec = seq( 1, Bmax, .001 ),  # RR scale
                                   .calib = dat$calib, # on log-RR scale (transformed if needed)
                                   .tail = tail[ind],
                                   .causative = (mu>0),
                                   
                                   .give.CI = TRUE,
                                   .R = boot.reps,
                                   .dat = dat,
                                   .calib.name = "calib" )
                 
                 return( data.frame( That = paste( round( TG$Est[ TG$Stat == "That" ], digits ),
                                                   format_CI(  TG$lo[ TG$Stat == "That" ],
                                                               TG$hi[ TG$Stat == "That" ],
                                                               digits) ),
                                     
                                     Ghat = paste( round( TG$Est[ TG$Stat == "Ghat" ], digits ),
                                                   format_CI(  TG$lo[ TG$Stat == "Ghat" ],
                                                               TG$hi[ TG$Stat == "Ghat" ],
                                                               digits) ) ) )
                 
               } )

TG.df = do.call( rbind,
                 TG.l )


##### Phat(B) Plots #####

# here we need to make a decision
# 1. always do calibrated plot with sigma^2B=0 so that we don't have to cut off plot
# 2. if effects are normal, also do parametric plot with sigma^2B>0, but cut off plot at 0.15 and 0.85

for (i in length(ql)) {
  
  # ASSUMES ONLY 2 VALUES OF SIGB AND THAT THE SECOND ONE IS HETEROGENEOUS
  # homogeneous bias vs. parametric Phat, showing full range of Phats
  sens_plot_calib( q = ql[[i]],
                   B.vec = seq(1, 1.3, .01),
                   tail = tail[[i]],
                   causative = (mu>0),
                   est.type = est.type,
                   boot.reps = boot.reps,
                   dat = dat )  # from analyze.R
  
  setwd(results.dir)
  ggsave( filename = paste( meta.name, "_", "sens_calib_expq", exp(ql[[i]]), "_homoB.png", sep="" ),
          width = 5,
          height = 4,
          units = "in")
  
  # heterogeneous bias vs. parametric Phat, but cutting off extreme Phats
  sens_plot_para( type = "line",
                  q = ql[[i]],
                  Bmin = log(1),
                  Bmax = log(Bmax),
                  sigB = sigB[2],
                  breaks.x1 = breaks.x1,
                  tail = tail[[i]],
                  
                  yr = meta$b,
                  vyr = meta$se^2,
                  t2 = meta$tau2,
                  vt2 = meta$se.tau2^2,
                  
                  est.type = est.type )
  
  setwd(results.dir)
  ggsave( filename = paste( meta.name, "_", "sens_para_expq", exp(ql[[i]]), "_heteroB.png", sep="" ),
          width = 5,
          height = 4,
          units = "in")
  
}


##### E-value for Point Estimate #####

Eval.l = lapply( X = ql,
                 FUN = function(.q){
                   evals = suppressMessages( EValue::evalues.RR( est = exp(mu),
                                                                 lo = exp(mu.lo),
                                                                 hi = exp(mu.hi),
                                                                 true = exp(.q) ) )
                   
                   return( data.frame( Eval.est = evals[2,1],
                                       # grab the CI E-value that isn't NA
                                       Eval.ci = evals[2,2:3][ !is.na(evals[2,2:3]) ] ) )
                   
                 } )

Eval.df = do.call( rbind,
                   Eval.l )


##### Put Results in Dataframe #####
if (take.exp == TRUE) {
  est = exp(mu)
  lo = exp(mu.lo)
  hi = exp(mu.hi)
} else {
  est = mu
  lo = mu.lo
  hi = mu.hi
}

est.string = paste( round( est, digits ),
                    format_CI( lo, 
                               hi,
                               digits),
                    sep = " " )

tau.string = round( sqrt(t2), digits)

new.row = data.frame( Meta = meta.name,
                      k = nrow(dat),
                      Est = est.string,
                      Pval = format_stat(mu.pval, cutoffs = c(.1, .0001) ),
                      Tau = tau.string
)

# tail is now just for string purposes
# tail = rep("above", length(unlist(ql)))
# tail[unlist(ql) < 0] = "below"
if (take.exp == TRUE) q.vec = exp(unlist(ql)) else q.vec = unlist(ql)
Phat.names = paste( "Percent ", tail, " ", q.vec, sep = "" )
# new.row[, Phat.names ] = NA

new.row[ , Phat.names ] = Phat.df$string

# That and Ghat
That.names = paste( "That ", tail, " ", q.vec, sep = "" )
Ghat.names = paste( "Ghat ", tail, " ", q.vec, sep = "" )

# bm
new.row[ , That.names ] = TG.df$That
new.row[ , Ghat.names ] = TG.df$Ghat

# E-value for point estimate and CI
Eval.est.names = paste( "Eval est ", q.vec, sep = "" )
Eval.CI.names = paste( "Eval CI ", q.vec, sep = "" )
new.row[ , Eval.est.names ] = round( Eval.df$Eval.est, digits )
new.row[ , Eval.CI.names ] = round( Eval.df$Eval.ci, digits )

# resE is, or will become, a global variable
if ( !exists("resE") ){
  resE <<- new.row
} else {
  library(plyr)
  resE <<- rbind.fill(resE, new.row)
  detach("package:plyr", unload=TRUE)
}