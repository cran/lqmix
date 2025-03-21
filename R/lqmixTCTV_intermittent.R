lqmixTCTV = function(formula, randomTC, randomTV, group, time, m, G, data, qtl=0.5, eps=10^-5, maxit=1000, se=TRUE, R=50, start=0, parInit = list(betaf=NULL,  betarTC=NULL, betarTV=NULL, pg=NULL, delta = NULL, Gamma = NULL, scale=NULL),verbose=TRUE, seed=NULL, parallel=parallel, ...){

  # start = 0 -- deterministic start
  # start = 1 -- random start
  # start = 2 -- start from given values

  # ---- possible errors ------
  # ***************************

  if(start == 2 & is.null(unlist(parInit))) stop("No input parameters have been given with start = 2.")

  if (qtl <= 0 | qtl >= 1) stop("Quantile level out of range")

  # ---- initial settings ----
  # ***************************
  if(!is.data.frame(data))  stop("`data' must be a data frame.")

  if(!inherits(formula, "formula") || length(formula) != 3) stop("\nFixed coefficient model must be a formula of the form \"y ~ x\".")
  if(!inherits(randomTV, "formula") || length(randomTV) != 2) stop("\nTV random coefficient model must be a formula of the form \"~ w\".")
  if (!inherits(randomTC, "formula") || length(randomTC) != 2) stop("\nTC random coefficient model must be a formula of the form \"~ z\".")

  if(m == 1 & G == 1) stop("The specified model corresponds to a linear quantile regression model with no random coefficients. Please, use the function lqr().")
  else if(m == 1) stop("The specified model corresponds to a linear quantile mixture with TC random coefficients only. Please, use the function lqmixTC().")
  else if(G == 1) stop("The specified model corresponds to a linear quantile mixture with TV random coefficients only. Please, use the function lqmixTV().")


  ranIntTC = length(grep("1", randomTC))>0
  ranIntTV = length(grep("1", randomTV))>0

  mformt = as.character(time)
  mform2time = gsub(" ", "", mformt)

  mformTC = strsplit(as.character(randomTC)[2], "\\|")[[1]]
  mformTC = gsub(" ", "", mformTC)
  formula.rTC = reformulate(mformTC[1], intercept = ranIntTC) # random formula terms
  mform2sbjTC = group # id variable

  mformTV = strsplit(as.character(randomTV)[2], "\\|")[[1]]
  mformTV = gsub(" ", "", mformTV)
  formula.rTV = reformulate(mformTV[1], intercept = ranIntTV) # random formula terms
  mform2sbjTV = group # id variable

  if (!(all(mform2sbjTC %in% names(data)))) stop("The specified clustering variable is not contained in the data frame.")
  if (!(all(mform2time %in% names(data)))) stop("The specified time variable is not contained in the data frame.")

  randomTermsTC = attr(terms(formula.rTC), "intercept") + length(attr(terms(formula.rTC),"term.labels"))>0
  randomTermsTV = attr(terms(formula.rTV), "intercept") + length(attr(terms(formula.rTV),"term.labels"))>0

  if(!(randomTermsTC & randomTermsTV)) stop("The specifid model corresponds to a linear quantile regression model with no random coefficients. Please, use the function lqr().")
  else if(!(randomTermsTC) & randomTermsTV) stop("The specifid model corresponds to a linear quantile regression model with TV random coefficients only. Please, use the function lqmixTV().")
  else if(randomTermsTC & !(randomTermsTV)) stop("The specifid model corresponds to a linear quantile regression model with TC random coefficients only. Please, use the function lqmixTC().")


  # ---- initial settings ----
  # ***************************
  # remove incomplete data
  names = c(unique(unlist(lapply(c(formula, formula.rTC, formula.rTV), all.vars))), mform2sbjTC, mform2time)
  asOneFormula = eval(parse(text = paste("~", paste(names, collapse = "+")))[[1]])
  data = model.frame(asOneFormula, data)

  # response
  namesY = formula[[2]]

  # define fixed and random model frames and matrices
  # fixed model frame
  mff = model.frame(formula, data)
  # random model frame TC
  mfrTC = model.frame(formula.rTC, data)
  # random model matrix TC
  mmrTC = model.matrix(formula.rTC, mfrTC)
  # random model frame TV
  mfrTV = model.frame(formula.rTV, data)
  # random model matrix TV
  mmrTV = model.matrix(formula.rTV, mfrTV)

  # intercept derived from the formula
  fixInt = attr(terms(mff), "intercept") == 1
  ranIntTC = attr(terms(mfrTC), "intercept") == 1
  ranIntTV = attr(terms(mfrTV), "intercept") == 1

  # terms specified as both TC and TV -- error
  randomTermsTC = c(ranIntTC*1, attr(terms(formula.rTC),"term.labels"))
  randomTermsTV = c(ranIntTV*1, attr(terms(formula.rTV),"term.labels"))

  wh = intersect(randomTermsTC, randomTermsTV)
  if(length(wh)>0) stop(paste("One or more terms are specified as both TC and TV random coefficients:", paste(wh, collapse=", ")))

  # if random intercept, remove it from fixed formula
  if((ranIntTC == 1 | ranIntTV == 1) & fixInt == 1) fixInt = 0

  # identify the type of intercept: 0 -> fixed, 1 -> random TC, 2 -> TV, 999 -> no intercept in the model
  ranInt = 999  # Default value
  if (ranIntTC == 1 & fixInt == 0) ranInt = 1 else if (ranIntTV == 1 & fixInt == 0) ranInt = 2 else if (ranIntTC == 0 & ranIntTV == 0 & fixInt == 1) ranInt = 0


  # ranInt = ifelse(ranIntTC == 1 & fixInt == 0, 1,
  #                 ifelse(ranIntTV == 1 & fixInt == 0, 2,
  #                 ifelse(ranIntTC == 0 & ranIntTV == 0 & fixInt == 1, 0, 999)))

  # variable names
  termsFix = attr(terms(formula),"term.labels")
  termsRanTC = attr(terms(formula.rTC),"term.labels")
  termsRanTV = attr(terms(formula.rTV),"term.labels")

  # slopes in common between random and fixed model formula
  termsFix = setdiff(termsFix, termsRanTC)

  # wh = which(termsFix %in% c(termsRanTC, termsRanTV))
  # if(length(wh)>0) termsFix = termsFix[-wh]

  # any fixed slope? + reformulate the fixed model formula
  fixed = FALSE

  if(length(termsFix)>0){
    formula = reformulate(termsFix, response=as.character(formula[2]))
    fixed = TRUE
    mmf = model.matrix(formula, mff)
    if(!fixInt) mmf = mmf[, -1, drop = FALSE]
    # if(! as.logical(fixInt)){
    #   mmf = as.matrix(mmf[,-1])
    #   namesFix = attr(mmf0, "dimnames")[[2]][-1]
    #   colnames(mmf) = namesFix
    #   fixed = TRUE
    # }
    # fixed = TRUE
  }else if(fixInt){
    formula = reformulate("1", response = as.character(formula[2]))
    mmf = model.matrix(formula, mff)
    fixed = TRUE
  } else formula = mmf0 = mmf = NULL
  namesFix = colnames(mmf)

  # any random slope?
  ranSlopeTC = length(termsRanTC) > 0
  # ranSlopeTC = FALSE
  # if(length(termsRanTC>0)) ranSlopeTC = 1

  ranSlopeTV = length(termsRanTV) > 0
  # ranSlopeTV = FALSE
  # if(length(termsRanTV>0)) ranSlopeTV = 1

  #if(! as.logical(fixInt)) namesFix = colnames(mmf0)[-1] else namesFix = colnames(mmf0)
  namesRanTC = colnames(mmrTC)
  namesRanTV = colnames(mmrTV)

  # ordering the dataset
  sbj.obs = data[,mform2sbjTV]
  ord = order(sbj.obs)
  data = data[ord,]
  sbj.obs = as.factor(data[,mform2sbjTV])
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
  order.time = order(time)
  time = time[order.time]
  sbj = sbj[order.time]
  observed = observed[order.time]

  # nObs objects
  order.time = order(time.obs)
  data = data[order.time,]
  time.obs = time.obs[order.time]
  sbj.obs = sbj.obs[order.time]


  # prepare covariates and response
  # *******************************
  # random covariates
  mfrTC = model.frame(formula.rTC, data)
  x.randomTC = model.matrix(formula.rTC, mfrTC)

  mfrTV = model.frame(formula.rTV, data)
  x.randomTV = model.matrix(formula.rTV, mfrTV)

  if(ranSlopeTV){
    TVcovar = function(cov, sbj){
      any(!tapply(cov, sbj.obs, function(x) all(x == x[1])))
    }
    xRand = if (ranInt == 1) x.randomTV[, -1, drop = FALSE] else as.matrix(x.randomTV)
    checkTVcovar = apply(xRand, 2, function(xx){TVcovar(xx, sbj.obs)})
    if(any(checkTVcovar)) stop(paste("TV random coefficients are only admitted for TC covariates."))
  }

  # fixed covariates
  if(!is.null(formula)){
    mff = model.frame(formula, data)
    x.fixed = model.matrix(formula, mff)
    if(!fixInt) x.fixed = x.fixed[, -1, drop = FALSE]
  }else x.fixed = NULL

  # response variable
  y.obs = model.response(mff)

  if(!("opt" %in% names(list(...)))) {
    oo = lqmixTCTV.fit(y=y.obs, x.fixed=x.fixed, namesFix=namesFix,
                       x.randomTC=x.randomTC, namesRanTC=namesRanTC,
                       x.randomTV=x.randomTV, namesRanTV=namesRanTV,
                       sbj.obs=sbj.obs, time.obs=time.obs, observed=observed, m=m, G=G, qtl=qtl, n=n, T=T, Ti=Ti, nObs=nObs,
                       order.time=order.time, ranInt=ranInt, ranSlopeTC=ranSlopeTC,ranSlopeTV=ranSlopeTV, fixed=fixed, start=start, eps=eps,
                       maxit=maxit, parInit=parInit,verbose=verbose,seed=seed)
    } else oo = parInit
    see = se

    if(se){
      if(verbose) cat("Computing standard errors ...\n")

      # number of parameters
      if(fixed | fixInt) pf = ncol(x.fixed) else pf = 0
      prTC = ncol(x.randomTC)
      prTV = ncol(x.randomTV)

      if(parallel==TRUE) cl = makeCluster(2) else cl = makeCluster(1)
      #registerDoParallel(cl)
      cc=0
      registerDoSNOW(cl)
      pb = txtProgressBar(max = R, style = 3)
      progress = function(x) setTxtProgressBar(pb,x)
      opts = list(progress = progress)

      ooo.se = foreach(cc = (1 : R), .options.snow = opts) %dopar%{
      #ooo.se = foreach(cc = (1 : R)) %dopar%{

        stopImplicitCluster()
        tries = 0
        while(tries <= (R*10)){
          tries = tries + 1

          if (tries == R*10)  stop("Standard errors may not be computed.")

          # build the complete data matrices of size (n*T)*(?)
          xx.fixed = matrix(, n*T, pf)
          xx.randomTC = matrix(, n*T, prTC)
          xx.randomTV = matrix(, n*T, prTV)
          yy = c(matrix(, n*T, 1))

          xx.fixed[observed,] = x.fixed
          xx.randomTC[observed,] = x.randomTC
          xx.randomTV[observed,] = x.randomTV
          yy[observed] = y.obs

          # build new data matrices based on sampled units
          if(!is.null(seed)) set.seed(seed*tries*cc) else set.seed(tries*cc)
          sample.unit = sample(1:n, n, replace=TRUE)
          whTmax = tapply(time.obs, sbj.obs, max)
          Ti.sample = Ti[sample.unit]
          T.sample = max(whTmax)
          nObs.sample = sum(Ti.sample)

          time = rep(1:T.sample, each = n)
          x.fix.sample = x.ranTC.sample = x.ranTV.sample = y.sample = c()

          for(t in 1:T.sample){
            x.fix.sample = rbind(x.fix.sample, matrix(matrix(xx.fixed[(time == t),], n, pf)[sample.unit,], n, pf))
            x.ranTC.sample = rbind(x.ranTC.sample, matrix(matrix(xx.randomTC[(time == t),], n, prTC)[sample.unit,], n, prTC))
            x.ranTV.sample = rbind(x.ranTV.sample, matrix(matrix(xx.randomTV[(time == t),], n, prTV)[sample.unit,], n, prTV))
            y.sample = c(y.sample, yy[(time == t)][sample.unit])
          }

          observed.sample = !is.na(y.sample)

          if(fixed | fixInt){
            x.fix.sample = matrix(x.fix.sample[observed.sample,], nrow = nObs.sample)
            colnames(x.fix.sample) = namesFix
          }else x.fix.sample = NULL

          x.ranTC.sample = matrix(x.ranTC.sample[observed.sample,], nrow = nObs.sample)
          colnames(x.ranTC.sample) = namesRanTC

          x.ranTV.sample = matrix(x.ranTV.sample[observed.sample,], nrow = nObs.sample)
          colnames(x.ranTV.sample) = namesRanTV

          y.sample = y.sample[observed.sample]
          sbj.obs.sample = sbj[observed.sample]
          time.obs.sample = time[observed.sample]
          order.time.sample = order(time.obs.sample[order(sbj.obs.sample)])

          oo.se = tryCatch(lqmixTCTV.fit(y=y.sample, x.fixed=x.fix.sample, namesFix=namesFix, x.randomTC=x.ranTC.sample, namesRanTC=namesRanTC,
                                         x.randomTV=x.ranTV.sample, namesRanTV=namesRanTV,
                                         sbj.obs=sbj.obs.sample, time.obs=time.obs.sample, observed=observed.sample, m=m, G=G, qtl=qtl, n=n, T=T.sample, Ti=Ti.sample,
                                         nObs=nObs.sample,
                                         order.time=order.time.sample, ranInt=ranInt, ranSlopeTC=ranSlopeTC, ranSlopeTV=ranSlopeTV, fixed=fixed,
                                         start=2, eps=eps,
                                         maxit=maxit, parInit=oo, verbose=FALSE), error = function(e){e})
          if(!is(oo.se, "error")) break

        }

        boot = list ()
        boot$betaf = oo.se$betaf
        ord = order(oo.se$betarTC[,1])
        boot$betarTC = oo.se$betarTC[ord, ]
        boot$pg = oo.se$pg[ord]

        ord = order(oo.se$betarTV[,1])
        boot$betarTV = oo.se$betarTV[ord,]
        boot$delta = oo.se$delta[ord]
        boot$Gamma = oo.se$Gamma[ord, ord]

        boot$scale = oo.se$scale

        return(boot)
      }
      close(pb)
      stopImplicitCluster()
      parallel::stopCluster(cl)
      rm(cl)

      se = apply(sapply(ooo.se,unlist), 1, sd)

      if(fixed | fixInt){
        oo$se.betaf = se[grep("betaf", names(se))]
        names(oo$se.betaf) = namesFix
      }

      oo$se.betarTC = matrix(se[grep("betarTC", names(se))], nrow = G)
      colnames(oo$se.betarTC) = namesRanTC
      rownames(oo$se.betarTC) = paste("Comp", 1:G, sep="")

      oo$se.betarTV = matrix(se[grep("betarTV", names(se))], nrow = m)
      colnames(oo$se.betarTV) = namesRanTV
      rownames(oo$se.betarTV) = paste("St", 1:m, sep="")

      oo$se.delta = se[grep("delta", names(se))]
      names(oo$se.delta) = paste("St", 1:m, sep="")
      oo$se.Gamma = matrix(se[grep("Gamma", names(se))], m, m)
      rownames(oo$se.Gamma) = paste("fromSt", 1:m, sep="")
      colnames(oo$se.Gamma) = paste("toSt", 1:m, sep="")
      oo$se.pg = se[grep("pg", names(se))]
      names(oo$se.pg) = paste("Comp", 1:G, sep="")
      oo$se.scale = se[grep("scale", names(se))]
      cat("\n")
    }

    if(any(Ti != T)){
      if((any(time.input != timeidx))) oo$miss = "non-monotone" else oo$miss = "monotone"
      #
      # {
      #
      #   if(verbose){
      #     message("\nData affected by non-monotone missingness: parameter estimates may be biased.")
      #    }
      # }
    }else oo$miss = "none"

    oo$prTC = oo$prTC = NULL
    oo$prTV = oo$prTC = NULL
    oo$pf = oo$pr = NULL
    oo$model ="TCTV"

    oo$call = match.call()
    class(oo) = "lqmix"

    #message(mess)
    return(oo)
}

