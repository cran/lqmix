lqr.fit = function(y,x.fixed,namesFix,qtl,nObs,verbose){

  # initial settings
  # ****************
  namesFix = gsub(":", ".", namesFix)

  # number of parameters for the longitudinal process
  pf = ncol(x.fixed)

  # ---- parameter estimation ----
  # ***********************************
  mod0 = suppressWarnings(rq(y ~ .-1, data = data.frame(x.fixed)))

  # estimated fixed coefficients
  betaf = mod0$coefficients

  Xbeta = x.fixed %*% betaf
  linear.predictor = c(Xbeta)

  # estimated scale parameter
  scale = mean(rho(x = (y - linear.predictor), qtl = qtl))

  # compute densities
  # ******************
  resid = y - linear.predictor
  fit  = dal(y, linear.predictor, scale, qtl = qtl)
  lk = sum(log(fit))

  sigmaErr = sqrt(varAL(scale, qtl))

  if(verbose){
    cat("------------|-------------|\n")
    cat("  iteration |      lk     |\n")
    cat("------------|-------------|\n")
    cat(sprintf("%11g", c(0, lk)), "\n", sep = " | ")
    cat("------------|-------------|\n")
  }

  npar = pf+1

  aic = -2*lk + (2*npar)
  bic = -2*lk + npar *(log(nObs))

  # arrange output
  names(betaf) = namesFix

  res = list()
  res$betaf = betaf
  res$scale = scale
  res$sigma.e = sigmaErr
  res$lk = lk
  res$npar = npar
  res$AIC = aic
  res$BIC = bic
  res$qtl = qtl
  res$nobs = nObs

  res$pf = pf
  return(res)
}

