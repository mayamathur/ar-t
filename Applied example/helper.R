
# NOT IN USE ANYMORE
# works with the global variable resE
analyze_one_meta = function( dat,
                             
                             meta.name,
                             take.exp,
                             normal,  # are true effects normal?
                             
                             ql,  # log scale; given as list
                             tail,  # vector with same length as ql
                             

                             boot.reps = 2000,
                             
                             r, # for That and Ghat
                             
                             Bmax,
                             sigB = NA,  # absolute values for sigB
                             varB.ratio = NA,  # alternatively, sigB as a proportion of t2
                             breaks.x1,
                             est.type,  # for labelling y-axis
                             
                             digits = 2) {
  
  # TEST ONLY
  dat = dfs[[1]]
  meta.name = "Flegal"
  ql = list( log(1), log(.9) )
  tail = c( "above", "below" )
  take.exp = TRUE
  boot.reps = 500
  digits = 2
  #sigB = c(0, .5)
  r = 0.10
  varB.ratio = c(0, 0.8)
  Bmax = 1.3
  breaks.x1 = seq(1, Bmax, .1)
  est.type = "HR"
  
  
  
  
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
  if ( any( !is.na(sigB) ) & any( !is.na(varB.ratio) ) ) stop("Only provide one of sigB and varB.ratio, not both")
  
  if ( any( !is.na(sigB) ) & any(sigB > t2) ) stop( paste( "Some sigB > naive t2 of",
                                     round(t2,2),
                                     sep = " ") )
  
  # calculate absolute sigB if given ratio 
  if ( any( !is.na(varB.ratio) ) ) sigB = sqrt(varB.ratio * t2)
  
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
  
  # string is for just the parametric estimate
  # ~~~ give this an "if" so that it gives the parametric one if normal and not too extreme
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
} 




################################ MISCELLANEOUS ################################

# give 2.5th, 97.5th, 10th, and 90h percentiles of B on RR scale
# given its assumed parameters
B_pct = function( pcts = c( .025, .975, .1, .9 ),
                  muB,  # log scale
                  sigB ) {
  
  data.frame( pct = pcts,
              B = exp( qnorm( p = pcts,
                              mean = muB,
                              sd = sigB ) ) )
}




# for reproducible manuscript-writing
# adds a row to the file "stats_for_paper" with a new statistic or value for the manuscript
# optionally, "section" describes the section of code producing a given result
update_result_csv = function( name,
                              section = NA,
                              value = NA,
                              print = FALSE ) {
  setwd(results.dir)
  
  new.rows = data.frame( name,
                         value = as.character(value),
                         section = as.character(section) )
  
  # to avoid issues with variable types when overwriting
  new.rows$name = as.character(new.rows$name)
  new.rows$value = as.character(new.rows$value)
  new.rows$section = as.character(new.rows$section)
  
  
  if ( "stats_for_paper.csv" %in% list.files() ) {
    res = read.csv( "stats_for_paper.csv",
                    stringsAsFactors = FALSE,
                    colClasses = rep("character", 3 ) )
    
    # if this entry is already in the results file, overwrite the
    #  old one
    if ( all(name %in% res$name) ) res[ res$name %in% name, ] = new.rows
    else res = rbind(res, new.rows)
  }
  
  if ( !"stats_for_paper.csv" %in% list.files() ) {
    res = new.rows
  }
  
  write.csv( res, 
             "stats_for_paper.csv",
             row.names = FALSE,
             quote = FALSE )
  
  # also write to Overleaf
  setwd(overleaf.dir)
  write.csv( res, 
             "stats_for_paper.csv",
             row.names = FALSE,
             quote = FALSE )
  
  if ( print == TRUE ) {
    View(res)
  }
}

############################### EFFECT-SIZE CONVERSIONS ############################### 

smd_to_logRR = function(smd) {
  0.91 * smd
}


logHR_to_logRR = function(logRR){
  log( ( 1 - 0.5^sqrt( exp(logRR) ) ) / ( 1 - 0.5^sqrt( 1 / exp(logRR) ) ) )
}

############################### FNS FOR PHAT CAUSAL ############################### 

# define transformation in a way that is monotonic over the effective range of B (>1)
# to avoid ggplot errors
g = Vectorize( function(x) {
  # # this is needed only for That_causal in case we pass something that is NA
  # if ( is.na(x) | is.null(x) ) return(NA)
  if (x < 1) return( x / 1e10 )
  x + sqrt( x^2 - x )
} )

###### Phat after shifting by bias factor B and using calibrated estimates #####
# ** assumes the dataset has "yi" and "vi" coded as log-RRs

# .q = threshold on log-RR scale
# .B = bias factor on *RR* scale
# .calib = the calibrated estimates on log-RR scale
# .tail = "above" q or "below" q
# .causative = was the meta-analytic RR > 1 (causative) or < 1 (protective)?
# .R = number of bootstrap iterates to run
# .dat = the dataset for bootstrapping
Phat_causal = function( .q,
                        .B,
                        .calib, # assumed on log scale
                        .tail,
                        .causative,
                        
                        .give.CI = TRUE,
                        .R = 2000,
                        .dat = NA ) {
  
  
  # confounding-adjusted calibrated estimates
  # adjust upward or downward based on whether the naive estimate was
  #  causative or protective
  if ( .causative == TRUE ) {
    calib.t = .calib - log(.B)
  } else {
    calib.t = .calib + log(.B)
  }
  
  # confounding-adjusted Phat
  if ( .tail == "above" ) Phat.t = mean( calib.t > .q )
  if ( .tail == "below" ) Phat.t = mean( calib.t < .q )
  
  
  if ( .give.CI == FALSE ) {
    return(Phat.t)
    
  } else {
    boot.res = suppressWarnings( boot( data = .dat,
                                       parallel = "multicore",
                                       R = .R, 
                                       statistic = Phat_causal_bt,
                                       # below arguments are being passed to Phat_causal_bt
                                       .q = .q,
                                       .B = .B,
                                       .tail = .tail,
                                       .causative = .causative ) )
    
    bootCIs = boot.ci(boot.res,
                      type="bca",
                      conf = 0.95 )
    
    lo = bootCIs$bca[4]
    hi = bootCIs$bca[5]
    SE = sd(boot.res$t)
    
    return( data.frame( Est = Phat.t,
                        SE = SE,
                        lo = lo, 
                        hi = hi ) )
  }
}



###### Simplified version of above for boot to call #####

# assumes the dataset has variables called "yi" and "vi"
Phat_causal_bt = function( original,
                           indices,
                           .q,
                           .B,
                           .tail,
                           .causative ) {
  
  # bootstrapped dataset
  b = original[indices,]
  
  # get new calibrated estimates
  b$calib = calib_ests( yi = b$yi,
                        sei = sqrt(b$vi) )
  
  phatb = Phat_causal( .q = .q, 
                       .B = .B,
                       .calib = b$calib, 
                       .tail = .tail,
                       .causative,
                       .give.CI = FALSE)
  return(phatb)
}


############################### FNS FOR THAT AND GHAT ############################### 

##### That and Ghat from grid search of Phat values #####

# for each of a vector of bias factors, calculates Phat causal and then finds the one
#  that's closest to threshold proportion, .r
That_causal = function( .q,  # on log-RR scale (transformed if needed)
                        .r,
                        .B.vec,  # RR scale
                        .calib, # on log-RR scale (transformed if needed)
                        .tail,
                        .causative,
                        
                        .give.CI = TRUE,
                        .R = 2000,
                        .dat,
                        .calib.name ) {
  
  # list of bias factor values
  Bl = as.list(.B.vec)
  
  Phat.t.vec = lapply( Bl,
                       FUN = function(B) Phat_causal( .q = .q, 
                                                      .B = B,
                                                      .calib = .calib,
                                                      .tail = .tail,
                                                      .causative = .causative,
                                                      .give.CI = FALSE ) )
  
  res = data.frame( B = .B.vec,
                    Phat.t = unlist(Phat.t.vec) )
  
  # find the value of B that makes Phat causal as close as possible to chosen threshold, r
  That = res$B[ which.min( abs( res$Phat.t - .r ) ) ]
  # transform to G-scale
  Ghat = g(That)
  
  
  if ( .give.CI == FALSE ) {
    
    return( data.frame( That, Ghat ) )
    
  } else {
    boot.res = suppressWarnings( boot( data = .dat,
                                       parallel = "multicore",
                                       R = .R, 
                                       statistic = That_causal_bt,
                                       # below arguments are being passed to get_stat
                                       .calib.name = .calib.name,
                                       .q = .q,
                                       .r = .r,
                                       .B.vec = .B.vec,
                                       .tail = .tail,
                                       .causative = .causative ) )
    
    #browser()
    
    tryCatch({
      bootCIs = boot.ci(boot.res,
                        type="bca",
                        conf = 0.95 )
      
      lo = bootCIs$bca[4]
      hi = bootCIs$bca[5]
      SE = sd(boot.res$t)
      Note = NA
      
      # this can happen if all resampled stats are the same (but no error)
      if ( is.null(bootCIs) ){
        # FOR SOME BIZARRE REASON, THIS LOOP (EVEN WHEN ENTERED) ACTUALLY ASSIGNS THE FIRST 3 TO BE NULL
        #  AND THE LAST TO BE NA
        lo <<- NA
        hi <<- NA
        SE <<- NA
        Note <<- "BCa failed"
      }
      
    }, err = function(error){
      lo <<- NA
      hi <<- NA
      SE <<- NA
      Note <<- err$message
    })
    
    if ( is.null(lo) ) lo2 = c(NA, NA) else lo2 = c(lo, g(lo))
    if ( is.null(hi) ) hi2 = c(NA, NA) else hi2 = c(hi, g(hi))

    return( data.frame( Stat = c("That", "Ghat"),
                        Est = c(That, Ghat),
                        SE = c(SE, NA),  # ~~~ for MM's future reference: for latter, could replace with delta method if SE were actually wanted
                        lo = lo2, 
                        hi = hi2,
                        Meaning = c("Bias factor required", "Confounding strength required"),
                        Boot.note = Note ) )
  }
}


##### Simplified version of the above for boot to call #####
That_causal_bt = function( original,
                           indices, 
                           .calib.name,
                           .causative,
                           .q,
                           .r,
                           .B.vec,
                           .tail ) {
  b = original[indices,]
  
  Bl = as.list(.B.vec)
  
  # calculate Phat for a vector of B
  Phat.t.vec = unlist( lapply( Bl,
                               FUN = function(B) Phat_causal( .q = .q, 
                                                              .B = B,
                                                              .calib = b[[.calib.name]],
                                                              .tail = .tail,
                                                              .causative = .causative,
                                                              .give.CI = FALSE ) ) )
  
  
  That = .B.vec[ which.min( abs( Phat.t.vec - .r ) ) ]
  return(That)
}



############################### SENSITIVITY PLOT FNS ############################### 

###### Modified from Evalue package #####

# changed breaks
# also added est.type argument to relabel y-axis
sens_plot_para = function( type,
                           q,
                           muB=NA,
                           Bmin=log(1),
                           Bmax=log(5),
                           sigB=0,
                           yr,
                           vyr=NA,
                           t2,
                           vt2=NA,
                           tail,
                           breaks.x1=NA,
                           breaks.x2=NA,
                           est.type = "RR",
                           CI.level=0.95 ) {
  
  ##### Check for Bad Input ######
  if ( type=="dist" ) {
    
    if( is.na(muB) ) stop("For type='dist', must provide muB")
    
    if ( ( length(muB) > 1 ) | ( length(sigB) > 1 ) ) {
      stop( "For type='dist', muB and sigB must be length 1")
    }
  }
  
  if ( type=="line" ) {
    
    if ( is.na(vyr) | is.na(vt2) ) {
      message( "No confidence interval because vyr or vt2 is NULL")
    }
  }
  
  ##### Distribution Plot ######
  if ( type=="dist" ) {
    
    # simulate confounded distribution
    reps = 10000
    RR.c = exp( rnorm( n=reps, mean=yr, sd=sqrt(t2) ) )
    
    # simulate unconfounded distribution
    Mt = ifelse( yr > 0, yr - muB, yr + muB )
    RR.t = exp( rnorm( n=reps, mean=Mt, sd=sqrt(t2-sigB^2) ) )
    
    # get reasonable limits for X-axis
    x.min = min( quantile(RR.c, 0.01), quantile(RR.t, 0.01) )
    x.max = max( quantile(RR.c, 0.99), quantile(RR.t, 0.99) )
    
    temp = data.frame( group = rep( c( "Observed", "True" ), each = reps ),
                       val = c( RR.c, RR.t ) )
    
    colors=c("black", "orange")
    p = ggplot2::ggplot( data=temp, aes(x=temp$val, group=temp$group ) ) +
      geom_density( aes( fill=temp$group ), alpha=0.4 ) +
      theme_bw() + xlab("Study-specific relative risks") +
      ylab("") + guides(fill=guide_legend(title=" ")) +
      scale_fill_manual(values=colors) +
      geom_vline( xintercept = exp(q), lty=2, color="red" ) +
      scale_x_continuous( limits=c(x.min, x.max), breaks = seq( round(x.min), round(x.max), 0.5) ) +
      ggtitle("Observed and true relative risk distributions")
    
    graphics::plot(p)
  }
  
  ##### Line Plot ######
  if ( type=="line" ) {
    
    # get mean bias factor values for a bunch of different B's
    t = data.frame( B = seq(Bmin, Bmax, .01), phat = NA, lo = NA, hi = NA )
    t$eB = exp(t$B)
    
    for ( i in 1:dim(t)[1] ) {
      # r is irrelevant here
      cm = confounded_meta(q, r=0.10, muB=t$B[i], sigB,
                           yr, vyr, t2, vt2,
                           CI.level=CI.level,
                           tail = tail)
      t$phat[i] = cm$Est[ cm$Value=="Prop" ]
      t$lo[i] = cm$CI.lo[ cm$Value=="Prop" ]
      t$hi[i] = cm$CI.hi[ cm$Value=="Prop" ]
    }
    
    # cut off values for which we don't want to calculate a parametric CI
    t = t %>% filter( phat > 0.15 & phat < 0.85 )
    # also resent Bmin and Bmax to cut off plot
    Bmin = min(t$B)
    Bmax = max(t$B)
      
      
    # compute values of g for the dual X-axis
    if ( any( is.na(breaks.x1) ) ) breaks.x1 = seq( exp(Bmin), exp(Bmax), .5 )
    if ( any( is.na(breaks.x2) ) ) breaks.x2 = round( breaks.x1 + sqrt( breaks.x1^2 - breaks.x1 ), 2)
    
    # define transformation in a way that is monotonic over the effective range of B (>1)
    # to avoid ggplot errors
    g = Vectorize( function(x) {
      if (x < 1) return( x / 1e10 )
      x + sqrt( x^2 - x )
    } )
    

    p = ggplot2::ggplot( t, aes(x=t$eB, y=t$phat ) ) + theme_bw() +
      scale_y_continuous( limits=c(0,1),
                          breaks=seq(0, 1, .1) ) +
      
      # scale_x_continuous(  
      #                      
      #                      sec.axis = sec_axis( ~ g(.),  # confounding strength axis
      #                                           name = "Minimum strength of both confounding RRs" ) ) +
     
      scale_x_continuous(  limits = c(min(breaks.x1),
                                      max(breaks.x1)),
                           breaks = breaks.x1,
                           sec.axis = sec_axis( ~ g(.),  # confounding strength axis
                                                name = "Minimum strength of both confounding RRs",
                                                breaks=breaks.x2 ) ) +
      geom_line(lwd=1.2) +
      xlab("Bias factor (RR scale)") +
      ylab( paste( ifelse( tail == "above",
                           paste( "Estimated proportion of studies with true", est.type, " >", round( exp(q), 3 ) ),
                           paste( "Estimated proportion of studies with true", est.type, " <", round( exp(q), 3 ) ) ) ) )
    
    # can't compute a CI if the bounds aren't there
    no.CI = any( is.na(t$lo) ) | any( is.na(t$hi) )
    
    if ( no.CI ) graphics::plot(p)
    else p + ggplot2::geom_ribbon( aes(ymin=t$lo, ymax=t$hi), alpha=0.15 )
    
  }
}

# # test it
# # bm
# dat = dfs[[1]]
# meta = rma.uni( yi = dat$yi, 
#                     vi = dat$vi, 
#                     method = "REML", 
#                     knha = TRUE )
# 
# breaks.x1 = seq(1, 1.20, .05)
# breaks.x2 = g(breaks.x1)
# sens_plot_para( type = "line",
#                 q = log(1),
#                 Bmin = log(1),
#                 Bmax = log(2),
#                 sigB = .1,
#                 breaks.x1 = breaks.x1,
#                 
#                 yr = meta$b,
#                 vyr = meta$se^2,
#                 t2 = meta$tau2,
#                 vt2 = meta$se.tau2^2,
#                 
#                 est.type = "HR" )


sens_plot_calib = function(q,
                           B.vec,
                           tail,
                           causative,
                           est.type,
                           boot.reps = 2000,
                           dat,
                           yi.name = "yi",
                           vi.name = "vi",
                           lower.x.jump = .1,  # axis tick mark spacings
                           upper.x.jump = .2 ){
  
  dat$calib = calib_ests( yi = dat$yi,
                          sei = sqrt(dat$vi) )
  
  
  Bl = as.list(B.vec)
  Bmax = max(B.vec)
  
  # again pass threshold and calibrated estimates on log-RR scale
  Phat.t.vec = lapply( Bl,
                       FUN = function(.B) Phat_causal( .q = q, 
                                                       .B = .B,
                                                       .calib = dat$calib,
                                                       .tail = tail,  # always use this because we're only considering q=0
                                                       .give.CI = FALSE,
                                                       .causative = causative ) )
  
  
  res = data.frame( B = B.vec,
                    Est = unlist(Phat.t.vec) )
  
  ##### Selective Bootstrapping #####
  
  # look at just the values of B at which Phat jumps
  #  this will not exceed the number of point estimates in the meta-analysis
  res.short = res[ diff(res$Est) != 0, ]
  
  # bootstrap a CI for each entry in res.short
  temp = res.short %>% rowwise() %>%
    do( Phat_causal( .q = q, 
                     .B = .$B,
                     .calib = dat$calib,
                     .tail = tail,
                     .give.CI = TRUE,
                     .dat = dat,
                     .R = boot.reps,
                     .causative = causative ) )
  
  # merge this with the full-length res dataframe, merging by Phat itself
  res = merge( res, temp, by.x = "Est", by.y = "Est")
  
  #browser()
  
  ##### Make Plot #####
  ggplot( data = res,
          aes( x = B,
               y = Est ) ) +
    theme_bw() +
    
    # # proprtion "r" line
    # geom_hline( yintercept = r, 
    #             lty = 2,
    #             color = "red" ) +
    
    # # That line
    # geom_vline( xintercept = That$Est[ That$Stat == "That" ], 
    #             lty = 2,
    #             color = "black" ) +
    
    scale_y_continuous( limits=c(0,1), breaks=seq(0, 1, .1)) +
    scale_x_continuous(  breaks = seq(1, Bmax, lower.x.jump),
                         sec.axis = sec_axis( ~ g(.),  # confounding strength axis
                                              name = "Minimum strength of both confounding RRs",
                                              breaks = seq(1, g(Bmax), upper.x.jump )) ) +
    geom_line(lwd=1.2) +
    
    xlab("Hypothetical bias factor in all studies (RR scale)") +
    
    ylab( paste( ifelse( tail == "above",
                         paste( "Estimated proportion of studies with true", est.type, " >", round( exp(q), 3 ) ),
                         paste( "Estimated proportion of studies with true", est.type, " <", round( exp(q), 3 ) ) ) ) ) +
    
    geom_ribbon( aes(ymin=res$lo, ymax=res$hi), alpha=0.15, fill = "black" ) 
}


# # test it
# sens_plot_calib( q = log(.9),
#                  B.vec = seq(1, 1.3, .01),
#                  tail = "below",
#                  causative = FALSE,
#                  est.type = "HR",
#                  boot.reps = 200,
#                  dat = dfs[[1]] )  # from analyze.R






