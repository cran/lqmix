lqmixTV = function(formula, randomTV, group, time, m, data, qtl=0.5, eps=10^-5, maxit=1000, se=TRUE, R=200, start=0, parInit = list(betaf=NULL,  betarTV=NULL, delta=NULL, Gamma=NULL, scale=NULL), verbose=TRUE, seed=NULL, parallel=FALSE, ncores=2, ...){

  # start = 0 -- deterministic start
  # start = 1 -- random start
  # start = 2 -- start from given values

  # ---- possible errors -----
  # ***************************

  if(start == 2 & is.null(unlist(parInit))) stop("No input parameters have been given with start = 2.")
  if(!inherits(formula, "formula") || length(formula) != 3) stop("\nFixed coefficient model must be a formula of the form \"y ~ x\".")
  if(!inherits(randomTV, "formula") || length(randomTV) != 2) stop("\nTV random coefficient model must be a formula of the form \"~ w\".")

  if (qtl <= 0 | qtl >= 1) stop("Quantile level out of range")
  if(m == 1) stop("The specified model corresponds to a linear quantile regression model with no random coefficients. Please, use the function lqr().")


  opt = ("opt" %in% names(list(...)))

  ranIntTV = length(grep("1", randomTV))>0

  mformt = as.character(time) # time variable
  mform2time = gsub(" ", "", mformt)
  mformTV = strsplit(as.character(randomTV)[2], "\\|")[[1]]
  mformTV = gsub(" ", "", mformTV)
  formula.rTV = reformulate(mformTV[1], intercept = ranIntTV) # random formula terms
  mform2sbj = group # id variable

  if (!(all(mform2sbj %in% names(data)))) stop("The specified grouping variable is not contained in the data frame.")
  if (!(all(mform2time %in% names(data)))) stop("The specified time variable is not contained in the data frame.")

  randomTerms = attr(terms(formula.rTV), "intercept") + length(attr(terms(formula.rTV),"term.labels"))>0
  if(!randomTerms) stop("The specifid model corresponds to a linear quantile regression model with no random coefficients. Please, use the function lqr().")

  # ---- initial settings ----
  # ***************************
  # remove incomplete data
  names = c(unique(unlist(lapply(c(formula, formula.rTV), all.vars))), mform2sbj, mform2time)
  asOneFormula = eval(parse(text = paste("~", paste(names, collapse = "+")))[[1]])
  data = model.frame(asOneFormula, data)

  # ordering the dataset wrt subjects
  origOrder = row.names(data)
  sbj.obs = data[,mform2sbj]
  origOrder.sbj = unique(sbj.obs)

  ord = order(sbj.obs)
  data = data[ord,]
  sbj.obs = as.factor(data[,mform2sbj])
  n = length(unique(sbj.obs))
  levels(sbj.obs) = 1:n
  sbj.obs = as.numeric(sbj.obs)
  Ti = table(sbj.obs)
  nObs = length(sbj.obs)

  # identify observed values
  time.obs = as.factor(data[,mform2time])
  T = length(unique(time.obs))
  levels(time.obs) = 1:T
  time.obs = as.numeric(time.obs)

  # check whether missingness corresponds to dropout
  timeidx = sort(unlist(sapply(Ti, function(xx){1:xx})))
  time.input = sort(time.obs)

  sbjtime.input = paste("id",sbj.obs, "t", time.obs, sep="")
  sbj = rep(1:n, each = T)
  time = rep(1:T, n)
  order.time = order(time)
  all = paste("id", sbj, "t", time, sep="")
  observed = (all %in% sbjtime.input) # ordered wrt units
  last.obs = cumsum(Ti)

  # order data wrt to time
  # n*T object
  time = time[order.time]
  sbj = sbj[order.time]
  observed = observed[order.time]

  # nObs objects
  order.time = order(time.obs)
  data = data[order.time,]
  time.obs = as.factor(data[,mform2time])
  sbj.obs = as.factor(data[,mform2sbj])
  revOrder = match(origOrder, row.names(data))
  revOrder.sbj = match(origOrder.sbj, unique(sbj.obs))
  levels(sbj.obs) = 1:n
  levels(time.obs) = 1:T
  time.obs = as.numeric(time.obs)


  # build response and design matrices
  xy = xyBuildTV(formula=formula, randomTV=formula.rTV, data=data)

  y.obs = xy$y.obs; x.fixed = xy$x.fixed; x.random = xy$x.random
  namesFix = xy$nameFix; namesRan = xy$namesRan
  fixed = xy$fixed; fixInt = xy$fixInt; ranInt = xy$ranInt; ranSlope = xy$ranSlope

  #check whether TV random coefficients are associated to TC covariates. If so, error
  if(ranSlope){
    TVcovar = function(cov, sbj){
      any(!tapply(cov, sbj.obs, function(x) all(x == x[1])))
    }
    # Select the correct columns (excluding intercept if ranInt == 1)
    xRand = if (ranInt == 1) x.random[, -1, drop = FALSE] else as.matrix(x.random)
    checkTVcovar = apply(xRand, 2, function(xx){TVcovar(xx, sbj.obs)})
    if(any(checkTVcovar)) stop(paste("TV random coefficients are only admitted for TC covariates."))
  }


  if(!opt)  {
  oo = suppressWarnings(lqmixTV.fit(y=y.obs, x.fixed=x.fixed, namesFix=namesFix, x.random=x.random, namesRan=namesRan,
                   sbj.obs=sbj.obs, time.obs=time.obs, observed=observed, m=m, qtl=qtl, n=n, T=T, Ti=Ti, nObs=nObs,
                   order.time=order.time, ranInt=ranInt, ranSlope=ranSlope, fixed=fixed, start=start, eps=eps,
                   maxit=maxit, parInit=parInit, verbose=verbose, seed=seed))

  }else oo = parInit
  see = se

  if(se){
    if(verbose & !opt) cat("Computing standard errors ...\n") else if(verbose & opt) cat("Computing standard errors for the optimal model...\n")

    # number of parameters
    if(fixed | fixInt) pf = ncol(x.fixed) else pf = 0
    pr = ncol(x.random)
    if(parallel==TRUE) cl = makeCluster(ncores) else cl = makeCluster(1)

    cc=0
    registerDoSNOW(cl)
    if(verbose){
      pb = txtProgressBar(max = R, style = 3)
      progress = function(x) setTxtProgressBar(pb,x)
      opts = list(progress = progress)
    }else opts = list()

    if(!is.null(seed)) set.seed(seed)
    ooo.se = foreach(cc = (1 : R), .options.snow = opts) %dorng% {
      tries = 0

      while(tries <= (R*10)){
        tries = tries + 1
        if (tries == R*10)  stop("Standard errors may not be computed.")

        xx.fixed = matrix(, n*T, pf)
        xx.random = matrix(, n*T, pr)
        yy = c(matrix(, n*T, 1))

        xx.fixed[observed,] = x.fixed
        xx.random[observed,] = x.random
        yy[observed] = y.obs

        # build new data matrices based on sampled units
        sample.unit = sample(1:n, n, replace=TRUE)
        whTmax = tapply(time.obs, sbj.obs, max)
        Ti.sample = Ti[sample.unit]
        T.sample = max(whTmax)
        nObs.sample = sum(Ti.sample)

        x.fix.sample = x.ran.sample = y.sample = c()
        for(t in 1:T.sample){
          x.fix.sample = rbind(x.fix.sample, matrix(matrix(xx.fixed[(time == t),], n, pf)[sample.unit,], n, pf))
          x.ran.sample = rbind(x.ran.sample, matrix(matrix(xx.random[(time == t),], n, pr)[sample.unit,], n, pr))
          y.sample = c(y.sample, yy[(time == t)][sample.unit])
        }

        observed.sample = !is.na(y.sample)
        if(fixed | fixInt){
          x.fix.sample = matrix(x.fix.sample[observed.sample,], nrow = nObs.sample)
          colnames(x.fix.sample) = namesFix
        }else x.fix.sample = NULL

        x.ran.sample = matrix(x.ran.sample[observed.sample,], nrow = nObs.sample)
        if(ranSlope) colnames(x.ran.sample) = namesRan

        y.sample = y.sample[observed.sample]
        sbj.obs.sample = sbj[observed.sample]
        time.obs.sample = time[observed.sample]
        order.time.sample = order(time[observed.sample][order(sbj.obs.sample)])

        oo.se = tryCatch(suppressWarnings(lqmixTV.fit(y=y.sample, x.fixed=x.fix.sample, namesFix=namesFix, x.random=x.ran.sample, namesRan=namesRan,
                         sbj.obs=sbj.obs.sample, time.obs=time.obs.sample, observed=observed.sample, m=m, qtl=qtl, n=n, T=T.sample, Ti=Ti.sample, nObs=nObs.sample,
                         order.time=order.time.sample, ranInt = ranInt, ranSlope=ranSlope, fixed=fixed, start=start, eps=eps,
                         maxit=maxit, parInit=oo, verbose=FALSE)), error = function(e){e})

        if(!is(oo.se, "error")) break

      }

      boot = list ()
      boot$betaf = oo.se$betaf
      ord = order(oo.se$betarTV[,1])
      boot$betar = oo.se$betarTV[ord,]
      boot$delta = oo.se$delta[ord]
      boot$Gamma = oo.se$Gamma[ord,ord]
      boot$scale = oo.se$scale

      return(boot)
    }
    if(verbose) close(pb)
    stopImplicitCluster()
    parallel::stopCluster(cl)
    rm(cl)

    vc = as.matrix(cov(t(sapply(ooo.se,unlist))))
    se = sqrt(diag(vc))
    varcov = list()
    if(fixed | fixInt){
      wh = grep("betaf", names(se))
      oo$se.betaf = se[wh]
      varcov$betaf = as.matrix(vc[wh, wh])
      names(oo$se.betaf) = colnames(varcov$betaf) = rownames(varcov$betaf) = namesFix
    }
    wh = grep("betar", names(se))
    oo$se.betarTV = matrix(se[wh], nrow = m)
    varcov$betarTV  = as.matrix(vc[wh, wh])

    colnames(varcov$betarTV) = rownames(varcov$betarTV) = c(sapply(namesRan, function(xx) paste(xx, paste("St", 1:m, sep=""), sep="_")))

    colnames(oo$se.betarTV) = namesRan
    rownames(oo$se.betarTV) = paste("St", 1:m, sep="")

    oo$se.delta = se[grep("delta", names(se))]
    names(oo$se.delta) = paste("St", 1:m, sep="")
    oo$se.Gamma = matrix(se[grep("Gamma", names(se))], m, m)
    rownames(oo$se.Gamma) = paste("fromSt", 1:m, sep="")
    colnames(oo$se.Gamma) = paste("toSt", 1:m, sep="")
    oo$se.scale = se[grep("scale", names(se))]

    oo$vcov = varcov
    cat("\n")
  }

  if(any(Ti != T)){
    if((any(time.input != timeidx))) oo$miss = "non-monotone" else oo$miss = "monotone"
  }else oo$miss = "none"

  oo$pf = oo$pr = NULL
  oo$model ="TV"

  oo$call = match.call()


  if(fixed | fixInt){
    oo$mmf = as.matrix(x.fixed[revOrder,])
    colnames(oo$mmf) = namesFix
  }
  oo$mmrTV = as.matrix(x.random[revOrder,])
  colnames(oo$mmrTV) = namesRan

  oo$y = y.obs[revOrder]

  oo$postTV = as.matrix(oo$postTV[revOrder.sbj,])
  rownames(oo$postTV) = origOrder.sbj

  class(oo) = "lqmix"

  return(oo)
}

