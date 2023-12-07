lqmixTC.fit = function(y, x.fixed, namesFix, x.random, namesRan, id, G, qtl, n, T, Ti, nObs,
                       order.time, ranInt, ranSlope, fixed, start,eps, maxit, parInit, verbose=TRUE,seed=NULL,parallel){


  # initial settings
  # ****************
  namesFix = gsub(":", ".", namesFix)
  namesRan = gsub(":", ".", namesRan)

  fixInt = ifelse(ranInt == 0, 1, 0)


  # number of parameters for the longitudinal process
  if(fixed | fixInt) pf = ncol(x.fixed) else pf = 0
  pr = ncol(x.random)


  # ---- parameter initialization ----
  # ***********************************

  if(fixed | fixInt){

    mod1 = lm(y ~ .-1, data = data.frame(x.fixed,x.random))
    betaf = mod1$coef[1:pf]
    names(betaf) = namesFix
    mod0 = lm(y ~ .-1, data = data.frame(x.random,x.fixed))

  }else{
    mod0 = lm(y ~ .-1, data = data.frame(x.random,x.fixed))
    betaf = NULL
    }

  scale = mean(rho(mod0$residuals, qtl = qtl))
  if(start == 0){

    # deterministic start for parameters related to latent variables
    sumMod0 = matrix(suppressWarnings(summary(mod0)$coef[1:pr,]),nrow = pr)
    tmp = matrix(cbind(sumMod0[,1] - 2 * sumMod0[,2], sumMod0[,1] + 2 * sumMod0[,2]), ncol = 2)
    betar = apply(tmp, 1, function(xx){seq(xx[1],xx[2], length = G)})
    colnames(betar) = namesRan

    # ----- mixture probs -----
    pg = rep(1/G, G)

  }else if(start == 1){
    if(!is.null(seed)) set.seed(seed)
    # random start for parameters related to latent variables
    betar = sapply(mod0$coefficients[1:pr], function(xx){sort(xx + rnorm(G,0,1))})
    colnames(betar) = namesRan

    pg = rep(1/G, G)
    num.pg = abs(pg + rnorm(G, 0, 0.4))
    pg = num.pg / sum(num.pg)

  } else{

    if(is.null(unlist(parInit))) stop("Initial parameter values must be given in input")
    # start from given values
    betaf = parInit$betaf
    betar = parInit$betarTC

    pg = parInit$pg
    scale = parInit$scale

  }


  # ---- build design matrices and vectors -- (from nObs to nObs*G) ----
  # *********************************************************************
  onesG = rep(1, G)

  # response variable
  yg = onesG %x% y

  # fixed intercept and covariates
  if(fixed | (ranInt == 0)){
    xf = onesG %x% x.fixed
    colnames(xf) = namesFix
  }else xf = NULL


  # compute densities
  # ******************
  if(fixed) Xbeta = matrix(x.fixed %*% betaf, nObs, G) else Xbeta = 0
  Zbg = x.random %*% t(betar)

  linear.predictor = c(Xbeta +  Zbg)

  resid = (yg - linear.predictor)
  resid = array(resid, c(nObs, G))
  minres = apply(resid^2, c(1), min)
  minres = array(minres %x% t(onesG), c(nObs,G))
  wgt = resid^2 == minres

  fitg  = matrix(dal(yg, linear.predictor, scale, qtl = qtl), nObs,G)
  fig = exp(rowsum(log(fitg), group=id))
  fig = pmax(fig, 10^-300)

  # compute likelihood
  out = lkComputeTC(pg=pg, fig=fig)
  lk = out$lk; li = out$li; lig=out$lig

  # start EM
  # *********

  if(verbose){
    cat("------------|-------------|-------------|-------------|\n")
    cat("  iteration |      G      |      lk     |   (lk-lko)  |\n")
    cat("------------|-------------|-------------|-------------|\n")
    cat(sprintf("%11g", c(0, G, lk, NA)), "\n", sep = " | ")
  }

  iter = 0; lk0 = lk
  while (((lk - lk0) > eps | iter == 0) & (iter <= maxit)){

    iter = iter +1; lk0 = lk

    # E step
    # *******
    post = postComputeTC(lig=lig, fig=fig, pg=pg, G=G, Ti=Ti, order.time=order.time)
    wig = post$wig
    Wig = post$Wig

    # M-step
    # ******
    pg = apply(wig, 2, sum)/n
    yg2 = c(yg - c(Zbg))
    if(fixed | fixInt){
      mod = suppressWarnings(rq(yg2 ~.-1, weights = c(Wig), data = data.frame(xf), tau = qtl))
      betaf = mod$coefficients
    }else betaf = NULL

    if(!is.null(betaf)) Xbeta = matrix(x.fixed %*% betaf, nObs, G) else Xbeta=0

    yg2 = matrix(yg - c(Xbeta), nObs, G)
    for (g in 1:G){
      mod = suppressWarnings(rq(c(yg2[,g]) ~ .-1, weights = Wig[,g], tau = qtl, data = data.frame(x.random)))
      betar[g,] = mod$coefficients
    }

    Zbg = x.random %*% t(betar)
    linear.predictor = c(Xbeta +  Zbg)
    scale = sum(c(Wig) * c(rho(x = yg - linear.predictor, qtl = qtl))) / sum(Wig)

    # compute densities
    # ******************
    fitg  = matrix(dal(yg, linear.predictor, scale, qtl = qtl), nObs,G)
    fig = exp(rowsum(log(fitg),group=id))
    fig = pmax(fig, 10^-300)

    # compute likelihood
    out = lkComputeTC(pg=pg, fig=fig)
    lk = out$lk; li = out$li; lig=out$lig

    if(verbose)
      if(iter/10 == floor(iter/10)) cat(sprintf("%11g", c(iter, G, lk, (lk -lk0))), "\n", sep = " | ")
  }

  if(verbose){
    if(iter/10 > floor(iter/10)) cat(sprintf("%11g", c(iter, G, lk, (lk -lk0))), "\n", sep = " | ")
    cat("------------|-------------|-------------|-------------|\n")
  }


  sigmaErr = sqrt(varAL(scale, qtl))
  npar = (pr)*G + pf + (G-1)+1

  aic = -2*lk + (2*npar)
  bic = -2*lk + npar *(log(n))

  # arrange output
  if(!is.null(betaf)) names(betaf) = namesFix
  colnames(betar) = namesRan

  rownames(betar) = paste("Comp", 1:nrow(betar), sep = "")
  names(pg) = paste("Comp", 1:nrow(betar), sep = "")

  res = list()
  res$betaf = betaf
  res$betarTC = betar
  res$pg = pg
  res$scale = scale
  res$sigma.e = sigmaErr
  res$lk = lk
  res$npar = npar
  res$AIC = aic
  res$BIC = bic
  res$qtl = qtl
  res$G = G
  res$nsbjs = n
  res$nobs = nObs

  return(res)
}

