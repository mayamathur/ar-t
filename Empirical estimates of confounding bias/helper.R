

quick_Phat = function(.q, .yi, .vyi) {

  .dat = data.frame( yi = .yi, vyi = .vyi )
  
  Phat = prop_stronger(q = .q,
                       yi.name = "yi",
                       vi.name = "vyi",
                       dat = .dat,
                       # choose tail depending on sign of q
                       tail = ifelse(.q < log(1), "below", "above") )
  
  return( data.frame( Phat = paste( round( 100 * Phat$est ),
                                    "% [",
                                    round( 100 * Phat$lo ),
                                    "%, ",
                                    round( 100 * Phat$hi ),
                                    "%]",
                                    sep = "" )
                      ) )
  
}


