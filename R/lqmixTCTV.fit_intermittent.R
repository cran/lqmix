lqmixTCTV.fit = function(y, x.fixed, namesFix, x.randomTC,x.randomTV, namesRanTC, namesRanTV, sbj.obs, time.obs, observed, m,G, qtl, n, T, Ti, nObs,
                         order.time, ranInt, ranSlopeTC,ranSlopeTV, fixed, start,eps, maxit, parInit, verbose=TRUE,seed=NULL){


  # initial settings
  # ****************

  # initial settings
  # ****************
  namesFix = gsub(":", ".", namesFix)
  namesRanTC = gsub(":", ".", namesRanTC)
  namesRanTV = gsub(":", ".", namesRanTV)

  fixInt = ifelse(ranInt == 0, 1, 0)
  ranIntTC = ifelse(ranInt == 1, 1, 0)
  ranIntTV = ifelse(ranInt == 2, 1, 0)

  # number of parameters for the longitudinal process
  if(fixed | fixInt) pf = ncol(x.fixed) else pf = 0
  prTC = ncol(x.randomTC)
  prTV = ncol(x.randomTV)


  # ---- parameter initialization ----
  # ***********************************

  if(fixed | fixInt){
    mod1 = lm(y ~ .-1, data = data.frame(x.fixed,x.randomTC, x.randomTV))
    betaf = mod1$coef[1:pf]
    mod0TC = lm(y ~ .-1, data = data.frame(x.randomTC,x.fixed,x.randomTV))
    mod0TV = lm(y ~ .-1, data = data.frame(x.randomTV,x.fixed,x.randomTC))
  }else{
    mod0TC = lm(y ~ .-1, data = data.frame(x.randomTC,x.randomTV))
    mod0TV = lm(y ~ .-1, data = data.frame(x.randomTV,x.randomTC))
    betaf = NULL
  }
  scale = mean(rho(mod0TV$residuals, qtl = qtl))

  if(start == 0){

    # deterministic start for parameters related to latent variables
    sumMod0TC = matrix(suppressWarnings(summary(mod0TC)$coef[1:prTC,]),nrow = prTC)
    sumMod0TV = matrix(suppressWarnings(summary(mod0TV)$coef[1:prTV,]),nrow = prTV)
    # sumMod0TC = suppressWarnings(summary(mod0TC)$coef)
    # sumMod0TV = suppressWarnings(summary(mod0TV)$coef)

    tmp = matrix(cbind(sumMod0TC[,1] - 2 * sumMod0TC[,2], sumMod0TC[,1] + 2 * sumMod0TC[,2]), ncol=2)
    betarTC = matrix(apply(tmp, 1, function(xx){seq(xx[1],xx[2], length = G)}), nrow = G) # matrix of size (G*prTV)
    colnames(betarTC) = namesRanTC

    tmp = matrix(cbind(sumMod0TV[,1] - 2 * sumMod0TV[,2], sumMod0TV[,1] + 2 * sumMod0TV[,2]), ncol=2)
    betarTV = matrix(apply(tmp, 1, function(xx){seq(xx[1],xx[2], length = m)}), nrow = m) # matrix of size (m*prTV)
    colnames(betarTV) = namesRanTV

    # -- initial probs --
    delta = rep(1/m, m)

    # -- transition probs --
    s1 = 9
    Gamma = matrix(1,m,m)
    diag(Gamma) = s1+1
    Gamma = Gamma/(m+s1)

    # mixture probs ----
    pg = rep(1/G, G)

  }else if(start == 1){
    if(!is.null(seed)) set.seed(seed)

    betarTC = matrix(sapply(mod0TC$coef[1:prTC], function(xx){sort(xx + rnorm(G,0,1))}), nrow = G)
    colnames(betarTC) = namesRanTC

    betarTV = matrix(sapply(mod0TV$coef[1:prTV], function(xx){sort(xx + rnorm(m,0,1))}), nrow = m)
    colnames(betarTV) = namesRanTV

    delta = rep(1/m, m)
    num.delta = abs(delta + rnorm(m, 0, 0.4))
    delta = num.delta / sum(num.delta)

    s1 = abs(9*rnorm(1))
    Gamma = matrix(1,m,m)
    diag(Gamma) = s1+1
    Gamma = Gamma/(m+s1)

    pg = runif(G)

    pg = pg/sum(pg)

  }else{

    if(is.null(unlist(parInit))) stop("Initial parameter values must be given in input")
    # start from given values
    betaf = parInit$betaf
    betarTC = parInit$betarTC
    betarTV = parInit$betarTV

    delta = parInit$delta
    Gamma = parInit$Gamma

    pg = parInit$pg

    scale = parInit$scale
  }

  # ---- build design matrices and vectors -- (from nObs to nObs*m) ----
  # *********************************************************************
  onesm = rep(1, m)
  onesG = rep(1, G)
  onesmG = rep(1, m*G)

  # response variable
  yhg = array(y, c(nObs, m, G))

  # fixed intercept and covariates
  if(fixed | fixInt){
    xf = onesmG %x% x.fixed
    colnames(xf) = namesFix
  }else x.fixed = xf = NULL


#   # # random covariates
  xrTC = rep(1, m) %x% x.randomTC
  xrTV = rep(1, G) %x% x.randomTV


  # compute densities
  # ******************
  if(!is.null(betaf)) Xbeta = array(x.fixed %*% betaf, c(nObs, m, G))  else Xbeta = 0
  Zbg = aperm(array(x.randomTC %*% t(betarTC), c(nObs, G, m)), c(1,3,2))
  Wbetah = array(x.randomTV %*% t(betarTV), c(nObs, m, G))

  linear.predictor = c(Xbeta + Zbg + Wbetah)

  resid = (yhg - linear.predictor)
  resid = array(resid, c(nObs, m, G))
  minres = apply(resid^2, 1, min)
  minres = array(minres %x% t(onesm) %x% t(onesG), c(nObs, m, G))
  wgt = resid^2 == minres

  fithg  = array(dal(yhg, linear.predictor, scale, qtl = qtl), c(nObs,m, G))

  # compute log-likelihood
  out = lkComputeTCTV(delta, Gamma, pg, fithg, m, G, sbj.obs, time.obs, n, T, observed)

  lk = out$lk; li = out$li; A = out$A

  # start EM
  # *********

  if(verbose){
    cat("------------|-------------|-------------|-------------|-------------|\n")
    cat("  iteration |      m      |      G      |      lk     |   (lk-lko)  |\n")
    cat("------------|-------------|-------------|-------------|-------------|\n")
    cat(sprintf("%11g", c(0, m, G, lk, NA)), "\n", sep = " | ")
  }

  iter = 0; lk0 = lk
  while (((lk - lk0) > eps | iter == 0) & (iter <= maxit)){
    iter = iter +1; lk0 = lk

    # E step
    # *******
    post = postComputeTCTV(A, li, delta, Gamma, pg, fithg, m, G, sbj.obs, time.obs, n, T, observed)

    uSingle = post$uSingle
    uCouple = post$uCouple
    etag = post$etag
    wgt = post$wgt

    # M-step
    # ******

    delta = colMeans(uSingle[time.obs == 1, ])

    num = apply(uCouple, c(2,3), sum, na.rm = T)
    den = apply(uCouple, c(2), sum, na.rm = T) %o% rep(1, m)
    Gamma = num/den
    pg = colMeans(etag)


    ## estimate fixed parameters in the longitudinal data model
    yhg2 = c(yhg - c(Zbg + Wbetah))

    if(fixed | fixInt){
      mod = rq(yhg2 ~.-1, weights = c(wgt), data = data.frame(xf), tau = qtl)
      betaf = mod$coef
    }else betaf = NULL
    if(!is.null(betaf)) Xbeta = array(x.fixed %*% betaf, c(nObs, m, G)) else Xbeta = 0

    # estimate TC parameters in the longitudinal data model
    yhg2 = yhg - c(Xbeta + Wbetah)
    for (g in 1:G){
      mod = rq(c(yhg2[,,g]) ~ .-1, weights = c(wgt[,,g]), tau = qtl, data = data.frame(xrTC))
      betarTC[g,] = mod$coefficients
    }
    Zbg = aperm(array(x.randomTC %*% t(betarTC), c(nObs, G, m)), c(1,3,2))

    ## estimate TV parameters in the longitudinal data model
    yhg2 = yhg - c(Xbeta + Zbg)
    for (h in 1:m){
      mod = rq(c(yhg2[,h,]) ~ .-1, weights = c(wgt[,h,]), tau = qtl, data = data.frame(xrTV))
      betarTV[h,] = mod$coefficients
    }
    Wbetah = array(x.randomTV %*% t(betarTV), c(nObs, m, G))

    linear.predictor = c(Xbeta +  Zbg + Wbetah)
    scale = sum(c(wgt) * c(rho(x = yhg - linear.predictor, qtl = qtl))) / sum(wgt)


    # compute density
    # ******************
    fithg  = array(dal(yhg, linear.predictor, scale, qtl = qtl), c(nObs,m, G))

    # compute log-likelihood
    out = lkComputeTCTV(delta, Gamma, pg, fithg, m, G, sbj.obs, time.obs, n, T, observed)
    lk = out$lk; li = out$li; A = out$A

    if(verbose)
      if(iter/10 == floor(iter/10)) cat(sprintf("%11g", c(iter, m, G, lk, (lk -lk0))), "\n", sep = " | ")
  }
  if(verbose){
    if(iter/10 > floor(iter/10)) cat(sprintf("%11g", c(iter, m, G, lk, (lk -lk0))), "\n", sep = " | ")
    cat("------------|-------------|-------------|-------------|-------------|\n")
  }

  sigmaErr = sqrt(varAL(scale, qtl))
  npar = prTC*G+ ranIntTV*m + pf + (m-1)*(m+1) + (G-1) + 1

  aic = -2*lk + (2*npar)
  bic = -2*lk + npar *(log(n))

  # arrange output
  if(!is.null(betaf)) names(betaf) = namesFix
  colnames(betarTC) = namesRanTC
  colnames(betarTV) = namesRanTV

  rownames(betarTC) = paste("Comp", 1:nrow(betarTC), sep="")
  rownames(betarTV) = paste("St", 1:nrow(betarTV), sep="")
  names(pg) = paste("Comp", 1:nrow(betarTC), sep="")
  names(delta) = paste("St", 1:nrow(betarTV), sep="")
  rownames(Gamma) = paste("fromSt", 1:nrow(betarTV), sep="")
  colnames(Gamma) = paste("toSt", 1:nrow(betarTV), sep="")

  res = list()
  res$betaf = betaf
  res$betarTC = betarTC
  res$betarTV = betarTV
  res$delta = delta
  res$Gamma = Gamma
  res$pg = pg
  res$scale = scale
  res$sigma.e = sigmaErr
  res$lk = lk
  res$npar = npar
  res$AIC = aic
  res$BIC = bic
  res$qtl = qtl
  res$m = m
  res$G = G
  res$nsbjs = n
  res$nobs = nObs

  return(res)
}

