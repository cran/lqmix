lqmixTV = function(formula, randomTV, group, time, m, data, qtl=0.5, eps=10^-5, maxit=1000, se=TRUE, R=50, start=0, parInit = list(betaf=NULL,  betarTV=NULL, delta=NULL, Gamma=NULL, scale=NULL), verbose=TRUE, seed=NULL, parallel=FALSE, ...){

  # start = 0 -- deterministic start
  # start = 1 -- random start
  # start = 2 -- start from given values

  # ---- possible errors -----
  # ***************************
  if(se==TRUE & is.null(R)){
    mess = "\n The number of bootstrap samples for computing standard errors has not been specified. The default value R=50 has been used."
  }else mess=NULL

  if(start == 2 & is.null(unlist(parInit))) stop("No input parameters have been given with start = 2.")
  if(!is.data.frame(data))  stop("`data' must be a data frame.")
  if(!inherits(formula, "formula") || length(formula) != 3) stop("\nFixed coefficient model must be a formula of the form \"y ~ x\".")
  if(!inherits(randomTV, "formula") || length(randomTV) != 2) stop("\nTV random coefficient model must be a formula of the form \"~ w\".")
  if (qtl <= 0 | qtl >= 1) stop("Quantile level out of range")
  if(m == 1) stop("The specified model corresponds to a linear quantile regression model with no random coefficients. Please, use the function lqr().")


  mformt = as.character(time) # time variable
  mform2time = gsub(" ", "", mformt)
  mformTV = strsplit(as.character(randomTV)[2], "\\|")[[1]]
  mformTV = gsub(" ", "", mformTV)
  formula.rTV = reformulate(mformTV[1]) # random formula terms
  mform2sbj = group # id variable

  if (!(all(mform2sbj %in% names(data)))) stop("The specified clustering variable is not contained in the data frame.")
  if (!(all(mform2time %in% names(data)))) stop("The specified time variable is not contained in the data frame.")

  randomTerms = attr(terms(formula.rTV), "intercept") + length(attr(terms(formula.rTV),"term.labels"))>0
  if(!randomTerms) stop("The specifid model corresponds to a linear quantile regression model with no random coefficients. Please, use the function lqr().")

  # ---- initial settings ----
  # ***************************
  # printing options
  if(verbose & length(as.logical(list(...)))==0) cat("Model TV - qtl =", qtl, "\n")

  # remove incomplete data
  names = c(unique(unlist(lapply(c(formula, formula.rTV), all.vars))), mform2sbj, mform2time)
  asOneFormula = eval(parse(text = paste("~", paste(names, collapse = "+")))[[1]])
  data = model.frame(asOneFormula, data)

  # response
  namesY = as.character(formula)[[2]]

  # define fixed and random model frames and matrices
  # fixed model frame
  mff = model.frame(formula, data)
  # random model frame
  mfr = model.frame(formula.rTV, data)
  # random model matrix
  mmr = model.matrix(formula.rTV, mfr)


  # intercept derived from the formula
  fixInt = attr(terms(mff), "intercept") == 1
  ranInt = attr(terms(mfr), "intercept") == 1
  # if random intercept, remove intercept from fixed formula
  if(ranInt == 1 & fixInt == 1) fixInt = 0

  # identify the type of intercept: 0 -> fixed, 1 -> random TC, 999 -> no intercept in the model
  ranInt = ifelse(ranInt == 1 & fixInt == 0, 1, ifelse(ranInt == 0 & fixInt == 1, 0, 999))

  # variable names
  termsFix = attr(terms(formula),"term.labels")
  termsRan = attr(terms(formula.rTV),"term.labels")
  # ATT: intercepts do not appear in the list

  # slopes in common between random and fixed model formula
  wh = which(termsFix %in% termsRan)
  if(length(wh)>0) termsFix = termsFix[-wh]

  # any fixed slope? + reformulate the fixed model formula
  fixed = FALSE
  if(length(termsFix)>0){
    formula = reformulate(termsFix, response=as.character(formula[2]), intercept = as.logical(fixInt))
    # fixed model matrix
    mmf = model.matrix(formula, mff)
    fixed = TRUE
  }else if(fixInt){
    formula = reformulate("1", response = as.character(formula[2]))
    mmf = model.matrix(formula, mff)
  } else formula = mmf = NULL

  # any random slope?
  ranSlope = FALSE
  if(length(termsRan>0)) ranSlope = 1

  # namesFix and namesRan
  namesFix = colnames(mmf)
  namesRan = colnames(mmr)

  # ordering the dataset
  sbj.obs = data[,mform2sbj]
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
  levels(sbj.obs) = 1:T
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
  # ********************************
  # random covariates
  mfr = model.frame(formula.rTV, data)
  x.random = model.matrix(formula.rTV, mfr)

  #check whether TV random coefficients are associated to TC covariates. Else error
  if(ranSlope){
    TVcovar = function(cov, sbj){
      tapply(cov, sbj, function(xx){length(table(xx))>1})
    }
    if(ranInt == 1) checkTVcovar = colSums(apply(as.matrix(x.random[,-1], nrow = nObs), 2, function(xx){TVcovar(xx, sbj.obs)}))>0
      else checkTVcovar = colSums(apply(as.matrix(x.random, nrow = nObs), 2, function(xx){TVcovar(xx, sbj.obs)}))>0
    if(any(checkTVcovar)) stop(paste("TV random coefficients are only admitted for TC covariates."))
  }

  # fixed covariates
  if(!is.null(formula)){
    mff = model.frame(formula, data)
    x.fixed = model.matrix(formula, mff)
  }else x.fixed = NULL

  # response variable
  #y.obs = data[,as.character(namesY)]
  y.obs = model.response(mff)
  if(!("opt" %in% names(list(...))))  {
  oo = suppressWarnings(lqmixTV.fit(y=y.obs, x.fixed=x.fixed, namesFix=namesFix, x.random=x.random, namesRan=namesRan,
                   sbj.obs=sbj.obs, time.obs=time.obs, observed=observed, m=m, qtl=qtl, n=n, T=T, Ti=Ti, nObs=nObs,
                   order.time=order.time, ranInt=ranInt, ranSlope=ranSlope, fixed=fixed, start=start, eps=eps,
                   maxit=maxit, parInit=parInit, verbose=verbose, seed=seed))
  }else oo = parInit
  see = se
  if(se){
    if(verbose) cat("Computing standard errors ... ")

    # number of parameters
    if(fixed | fixInt) pf = ncol(x.fixed) else pf = 0
    pr = ncol(x.random)
    if(parallel==TRUE) cl = makeCluster(2) else cl = makeCluster(1)
    registerDoParallel(cl)

    #if(parallel==TRUE) registerDoParallel(2) else registerDoParallel(1)
    count=0
    ooo.se = foreach(count = (1 : R)) %dopar%{

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
        set.seed(tries*count)
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

      # if(verbose) cat(count, " ... ")
      # if ("search" %in% names(list(...))) cat(count, " ... ")

      boot = list ()
      boot$betaf = oo.se$betaf
      ord = order(oo.se$betarTV[,1])
      boot$betar = oo.se$betarTV[ord,]
      boot$delta = oo.se$delta[ord]
      boot$Gamma = oo.se$Gamma[ord,ord]
      boot$scale = oo.se$scale

      return(boot)
    }
    stopImplicitCluster()
    parallel::stopCluster(cl)
    rm(cl)
    se = apply(sapply(ooo.se,unlist), 1, sd)


    if(fixed | fixInt){
      oo$se.betaf = se[grep("betaf", names(se))]
      names(oo$se.betaf) = namesFix
    }

    oo$se.betarTV = matrix(se[grep("betar", names(se))], nrow = m)
    colnames(oo$se.betarTV) = namesRan
    rownames(oo$se.betarTV) = paste("St", 1:m, sep="")

    oo$se.delta = se[grep("delta", names(se))]
    names(oo$se.delta) = paste("St", 1:m, sep="")
    oo$se.Gamma = matrix(se[grep("Gamma", names(se))], m, m)
    rownames(oo$se.Gamma) = paste("fromSt", 1:m, sep="")
    colnames(oo$se.Gamma) = paste("toSt", 1:m, sep="")
    oo$se.scale = se[grep("scale", names(se))]

    cat("\n")
  }

  if(any(Ti != T)){
    if((any(time.input != timeidx))){
      oo$miss = "non-monotone"
      if(verbose){
       message("\nData affected by non-monotone missingness: parameter estimates may be biased.")
      }
    }else oo$miss = "monotone"
  }else oo$miss = "none"

  oo$pf = oo$pr = NULL
  oo$model ="TV"

  oo$call = match.call()
  class(oo) = "lqmix"


 # set.seed(NULL)
  return(oo)
  message(mess)

}

