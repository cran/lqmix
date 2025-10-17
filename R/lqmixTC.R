lqmixTC  = function(formula, randomTC, group, time, G, data, qtl=0.5, eps=10^-5, maxit=1000, se=TRUE, R=200, start=0, parInit = list(betaf=NULL, betarTC=NULL, pg=NULL, scale=NULL), verbose=TRUE, seed=NULL, parallel=FALSE, ncores=2, ...){

  # start = 0 -- deterministic start
  # start = 1 -- random start
  # start = 2 -- start from given values


  # ---- possible errors ------
  # ***************************

  if(start == 2 & is.null(unlist(parInit))) stop("No input parameters have been given with start = 2.")
  if(!inherits(formula, "formula") || length(formula) != 3) stop("\nFixed coefficient model must be a formula of the form \"y ~ x\".")
  if (!inherits(randomTC, "formula") || length(randomTC) != 2) stop("\nTC random coefficient model must be a formula of the form \"~ z\".")
  if (qtl <= 0 | qtl >= 1) stop("Quantile level out of range")

  if(G == 1) stop("The specified model corresponds to a linear quantile regression model with no random coefficients. Please, use the function lqr().")

  opt = ("opt" %in% names(list(...)))

  ranIntTC = length(grep("1", randomTC))>0

  mformt = as.character(time) # time variable
  mform2time = gsub(" ", "", mformt)
  mformTC = strsplit(as.character(randomTC)[2], "\\|")[[1]]
  mformTC = gsub(" ", "", mformTC)
  formula.rTC = reformulate(mformTC[1], intercept = ranIntTC) # random formula terms
  mform2sbj = as.character(group) # id variable

  if (!(all(mform2sbj %in% names(data)))) stop("The specified grouping variable is not contained in the data frame.")
  if (!(all(mform2time %in% names(data)))) stop("The specified time variable is not contained in the data frame.")

  randomTerms = attr(terms(formula.rTC), "intercept") + length(attr(terms(formula.rTC),"term.labels"))>0
  if(!randomTerms) stop("The specifid model corresponds to a linear quantile regression model with no random coefficients. Please, use the function lqr().")

  # ---- initial settings ----
  # ***************************
  # remove incomplete data
  names = c(unique(unlist(lapply(c(formula, formula.rTC), all.vars))), mform2sbj, mform2time)
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
  # n*T objects
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
  xy = xyBuildTC(formula=formula, randomTC=formula.rTC, data=data)

  y.obs = xy$y.obs; x.fixed = xy$x.fixed; x.random = xy$x.random
  namesFix = xy$nameFix; namesRan = xy$namesRan
  fixed = xy$fixed; fixInt = xy$fixInt; ranInt = xy$ranInt; ranSlope = xy$ranSlope

  if(!opt) {
  oo = lqmixTC.fit(y=y.obs, x.fixed=x.fixed, namesFix=namesFix, x.random=x.random, namesRan=namesRan,
                   id=sbj.obs, G=G, qtl=qtl, n=n, T=T, Ti=Ti, nObs=nObs,
                   order.time=order.time, ranInt=ranInt, ranSlope=ranSlope, fixed=fixed, start=start, eps=eps,
                   maxit=maxit, parInit=parInit,verbose=verbose,seed=seed)
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

        # build the complete data matrices of size (n*T)*(?)
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
        time = rep(1:T.sample, each = n)

        x.fix.sample = x.ran.sample = y.sample = c()
        for(t in 1:T){
          x.fix.sample = rbind(x.fix.sample, matrix(matrix(xx.fixed[(time == t),], n, pf)[sample.unit,], n, pf))
          x.ran.sample = rbind(x.ran.sample, matrix(matrix(xx.random[(time == t),], n, pr)[sample.unit,], n, pr))
          y.sample = c(y.sample, yy[(time == t)][sample.unit])
        }

        observed.sample = !is.na(y.sample)
        nObs.sample = sum(observed.sample)

        if(fixed | fixInt){
          x.fix.sample = matrix(x.fix.sample[observed.sample,], nrow = nObs.sample)
          colnames(x.fix.sample) = namesFix
        }else x.fix.sample = NULL

        x.ran.sample = matrix(x.ran.sample[observed.sample,], nrow = nObs.sample)
        colnames(x.ran.sample) = namesRan

        y.sample = y.sample[observed.sample]
        sbj.obs.sample = sbj[observed.sample]
        order.time.sample = order(time[observed.sample][order(sbj.obs.sample)])

        oo.se = tryCatch(lqmixTC.fit(y=y.sample, x.fixed=x.fix.sample, namesFix=namesFix, x.random=x.ran.sample, namesRan=namesRan,
                                   id=sbj.obs.sample, G=G, qtl=qtl,
                                   n=n, T=T.sample, Ti=Ti.sample, nObs=nObs.sample,
                                   order.time=order.time.sample,ranInt=ranInt,ranSlope=ranSlope, fixed=fixed, start=2, eps=eps,
                                   maxit=maxit, parInit=oo, verbose=FALSE), error=function(e){e})
        if(!is(oo.se, "error")) break
      }

      boot = list ()
      boot$betaf = oo.se$betaf
      ord = order(oo.se$betar[,1])
      boot$betar = oo.se$betar[ord,]
      boot$pg = oo.se$pg[ord]
      boot$scale = oo.se$scale
      return(boot)
    }

    if(verbose) close(pb)
    stopImplicitCluster()
    parallel::stopCluster(cl)
    rm(cl)

    vc = cov(t(sapply(ooo.se,unlist)))
    se = sqrt(diag(vc))
    varcov = list()
    if(fixed | fixInt){
      wh = grep("betaf", names(se))
      oo$se.betaf = se[wh]
      varcov$betaf = as.matrix(vc[wh, wh])
      names(oo$se.betaf) = colnames(varcov$betaf) = rownames(varcov$betaf) = namesFix
    }
    wh = grep("betar", names(se))
    oo$se.betarTC = matrix(se[wh], nrow = G)
    varcov$betarTC  = as.matrix(vc[wh, wh])
    colnames(varcov$betarTC) = rownames(varcov$betarTC) = c(sapply(namesRan, function(xx) paste(xx, paste("Comp", 1:G, sep=""), sep="_")))

    colnames(oo$se.betarTC) = namesRan
    rownames(oo$se.betarTC) = paste("Comp", 1:G, sep="")

    oo$se.pg = se[grep("pg", names(se))]
    names(oo$se.pg) = paste("Comp", 1:G, sep="")
    oo$se.scale = se[grep("scale", names(se))]
    cat("\n")
    oo$vcov = varcov
  }

  if(any(Ti != T)){
    if((any(time.input != timeidx))){
      oo$miss = "non-monotone"
    }else oo$miss = "monotone"
  }else oo$miss = "none"

  oo$pf = oo$pr = NULL
  oo$model = "TC"

  oo$call = match.call()
  if(fixed | fixInt){
    oo$mmf = as.matrix(x.fixed[revOrder,])
    colnames(oo$mmf) = namesFix
  }
  oo$mmrTC = as.matrix(x.random[revOrder,])
  colnames(oo$mmrTC) = namesRan

  oo$y = y.obs[revOrder]

  oo$postTC = as.matrix(oo$postTC[revOrder.sbj,])
  rownames(oo$postTC) = origOrder.sbj


  class(oo) = "lqmix"
  return(oo)
}


