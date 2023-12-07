lqmixTV.fit = function(y, x.fixed, namesFix, x.random, namesRan, sbj.obs, time.obs, observed, m, qtl, n, T, Ti, nObs,
                       order.time, ranInt, ranSlope, fixed, start,eps, maxit, parInit, verbose=TRUE,seed=NULL){


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
    mod0 = lm(y ~ .-1, data = data.frame(x.random))
    betaf = NULL
  }
  scale = mean(rho(mod0$residuals, qtl = qtl))

  if(start == 0){

    # deterministic start for parameters related to latent variables
    sumMod0 = matrix(suppressWarnings(summary(mod0)$coef[1:pr,]),nrow = pr)
    tmp = matrix(cbind(sumMod0[,1] - 2 * sumMod0[,2], sumMod0[,1] + 2 * sumMod0[,2]), ncol = 2)
    betar = apply(tmp, 1, function(xx){seq(xx[1],xx[2], length = m)})
    colnames(betar) = namesRan

    # -- initial probs --
    delta = rep(1/m, m)

    # -- transition probs --
    s1 = 9
    Gamma = matrix(1,m,m)
    diag(Gamma) = s1+1
    Gamma = Gamma/(m+s1)


  }else if(start == 1){
    if(!is.null(seed)) set.seed(seed)

    # random start for parameters related to latent variables
    betar = sapply(mod0$coefficients[1:pr], function(xx){sort(xx + rnorm(m,0,1))})
    colnames(betar) = namesRan
    delta = rep(1/m, m)
    num.delta = abs(delta + rnorm(m, 0, 0.4))
    delta = num.delta / sum(num.delta)

    s1 = abs(9*rnorm(1))
    Gamma = matrix(1,m,m)
    diag(Gamma) = s1+1
    Gamma = Gamma/(m+s1)

  }else{

    if(is.null(unlist(parInit))) stop("Initial parameter values must be given in input")

    # start from given values
    betaf = parInit$betaf
    betar = parInit$betarTV

    delta = parInit$delta
    Gamma = parInit$Gamma

    scale = parInit$scale
  }

  # ---- build design matrices and vectors -- (from nObs to nObs*m) ----
  # *********************************************************************
  onesm = rep(1, m)

  # response variable
  yh = onesm %x% y

  # fixed intercept and covariates
  if(fixed | (ranInt == 0)){
    xf = onesm %x% x.fixed
    colnames(xf) = namesFix
  }else xf = NULL


  # compute densities
  # ******************
  if(fixed) Xbeta = matrix(x.fixed %*% betaf, nObs, m) else Xbeta = 0
  Wbetah = x.random %*% t(betar)

  linear.predictor = c(Xbeta + Wbetah)

  resid = (yh - linear.predictor)
  resid = array(resid, c(nObs, m))
  minres = apply(resid^2, c(1), min)
  minres = array(minres %x% t(onesm), c(nObs,m))
  wgt = resid^2 == minres

  fith  = matrix(dal(yh, linear.predictor, scale, qtl = qtl), nObs,m)

  # compute log-likelihood
  out = lkComputeTV_intermittent(delta, Gamma, fith, m, sbj.obs, time.obs, n, T, observed)
  lk = out$lk; li = out$li; A = out$A


  # start EM
  # *********

  if(verbose){
    cat("------------|-------------|-------------|-------------|\n")
    cat("  iteration |      m      |      lk     |   (lk-lko)  |\n")
    cat("------------|-------------|-------------|-------------|\n")
    cat(sprintf("%11g", c(0, m, lk, NA)), "\n", sep = " | ")
  }

  iter = 0; lk0 = lk
  while (((lk - lk0) > eps | iter == 0) & (iter <= maxit)){

    iter = iter +1; lk0 = lk

    # E step
    # *******
    post = postComputeTV_intermittent(A, li, delta, Gamma, fith, m, sbj.obs, time.obs, observed)
    uSingle = post$uSingle
    uCouple = post$uCouple

    # M-step
    # ******
    delta = colMeans(uSingle[time.obs == 1, ])

    num = apply(uCouple, c(2,3), sum, na.rm = T)
    den = apply(uCouple, c(2), sum, na.rm = T) %o% rep(1, m)
    Gamma = num/den

    yh2 = c(yh - c(Wbetah))
    if(fixed | fixInt){
      mod = suppressWarnings(rq(yh2 ~.-1, weights = c(uSingle), data = data.frame(xf), tau = qtl))
      betaf = mod$coefficients
    }else betaf = NULL

    if(!is.null(betaf)) Xbeta = matrix(x.fixed %*% betaf, nObs, m) else Xbeta=0

    yh2 = matrix(yh - c(Xbeta), nObs, m)
    for(h in 1:m){
      mod = suppressWarnings(rq(c(yh2[,h]) ~ .-1, weights = uSingle[,h], tau = qtl, data = data.frame(x.random)))
      betar[h,] = mod$coefficients
    }
    Wbetah = x.random %*% t(betar)

    linear.predictor = c(Xbeta +  Wbetah)
    scale = sum(c(uSingle) * c(rho(x = yh - linear.predictor, qtl = qtl))) / sum(uSingle)

    # compute densities
    # ******************
    fith  = matrix(dal(yh, linear.predictor, scale, qtl = qtl), nObs,m)
    fih = exp(rowsum(log(fith),group=sbj.obs))

    # compute log-likelihood
    out = lkComputeTV_intermittent(delta, Gamma, fith, m, sbj.obs, time.obs, n, T, observed)
    lk = out$lk; li = out$li; A = out$A


    if(verbose)
      if(iter/10 == floor(iter/10)) cat(sprintf("%11g", c(iter, m, lk, (lk -lk0))), "\n", sep = " | ")
  }
  if(verbose){
    if(iter/10 > floor(iter/10)) cat(sprintf("%11g", c(iter, m, lk, (lk -lk0))), "\n", sep = " | ")
    cat("------------|-------------|-------------|-------------|\n")
  }

  sigmaErr = sqrt(varAL(scale, qtl))
  npar = (pr)*m + pf + (m-1)*(m+1) +1

  aic = -2*lk + (2*npar)
  bic = -2*lk + npar *(log(n))

  # arrange output
  if(fixed) names(betaf) = namesFix
  colnames(betar) = namesRan

  rownames(betar) = paste("St", 1:nrow(betar), sep="")
  names(delta) = paste("St", 1:nrow(betar), sep="")
  rownames(Gamma) = paste("fromSt", 1:nrow(betar), sep="")
  colnames(Gamma) = paste("toSt", 1:nrow(betar), sep="")

  res = list()
  res$betaf = betaf
  res$betarTV = betar
  res$delta = delta
  res$Gamma = Gamma
  res$scale = scale
  res$sigma.e = sigmaErr
  res$lk = lk
  res$npar = npar
  res$AIC = aic
  res$BIC = bic
  res$qtl = qtl
  res$m = m
  res$nsbjs = n
  res$nobs = nObs
  return(res)
}
