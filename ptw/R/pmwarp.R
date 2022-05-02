pmwarp <- function (ref, samp, optim.crit, init.coef, try = FALSE,
                    mode = c("forward", "backward"), alg = c("ptw", "sptw"),
                    trwdth, trwdth.res, smooth.param, ndx=40, ...)
{
  mode <- match.arg(mode)
  alg <- match.arg(alg, c("ptw","sptw"))
  
  ## Multiply coefficients to prevent them from becoming too 
  ## small for numeric precision.
  if (alg == "ptw"){
    n <- length(init.coef)
    ncr <- ncol(ref)
    time <- (1:ncr)/ncr
    B <- matrix(time, nrow = ncr, ncol = n)
    B <- t(apply(B, 1, cumprod))/B
    a <- init.coef * ncr^(0:(n - 1))
  } else if (alg == "sptw"){
    if (is.matrix(samp)){
      m <- max(dim(samp), dim(ref))
    } else m <- max(length(samp), length(ref))
    time <- 1:m
    bdeg = 3 # degree of the B-splines
    B <- cbind(time, JOPS::bbase(time, 1, max(time), ndx, bdeg))
  
    n = ncol(B) # number of basis functions
    # If a not defined in function call, initialize warping coefficients of the
    # B-splines to 0 and the coefficient of the linear basis function to 1
    # if (missing(a))
    a <- c(1,rep(0,ndx+bdeg))
  }

  if (optim.crit == "RMS" & smooth.param > 0) {
    samp.sm <- t(apply(samp, 1, difsm, smooth.param))
    ref.sm <- t(apply(ref, 1, difsm, smooth.param))
  }
  
  if (!try) { # perform optimization
    switch(optim.crit,
           RMS = {
             if (smooth.param > 0) {
               Opt <- optim(a, RMS, NULL, ref.sm, samp.sm, B, mode = mode, ...)
             } else {
               Opt <- optim(a, RMS, NULL, ref, samp, B, mode = mode, ...)
             }},
           WCC = {
             wghts <- 1 - (0:trwdth)/trwdth
             ref.acors <- apply(ref, 1, wac, trwdth = trwdth, wghts = wghts)
             Opt <- optim(a, WCC, NULL, ref, samp, B,
                          trwdth = trwdth, wghts = wghts,
                          ref.acors = ref.acors, mode = mode, ...)
           })
    
    a <- c(Opt$par)
    v <- Opt$value

    ## if the optimization is done with a different smoothing or a
    ## different triangle, the final value for the optimization
    ## criterion is recalculated using the "original" data
    if ((optim.crit == "RMS" && smooth.param > 0) ||
        (optim.crit == "WCC" && trwdth != trwdth.res)) {
      v <- switch(optim.crit,
                  RMS = RMS(a, ref, samp, B, mode = mode),
                  WCC = WCC(a, ref, samp, B, trwdth.res, mode = mode))
    }
  }

  ## calculate, or possibly re-calculate, quality of current solution
  if (try) {
    if (optim.crit == "WCC") {
      v <- WCC(a, ref, samp, B, trwdth.res, mode = mode)
    } else {
      if (smooth.param > 0) {
        v <- RMS(a, ref.sm, samp.sm, B, mode = mode)
      } else {      
        v <- RMS(a, ref, samp, B, mode = mode)
      }
    }
  }

  ## back-transform coefficients
  w <- B %*% a
  a <- a/ncr^(0:(n-1))
  
  list(w = w, a = a, v = v)
} 
