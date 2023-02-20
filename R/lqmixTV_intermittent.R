#' Linear Quantile Mixture with Time-Varying (TV) Random Coefficients
#'
#' Estimate a finite mixture of linear quantile regression models with Time-Varying (TV), discrete, random coefficients
#'
#' @param formula an object of class formula of the form y ~ x1 + x2 + ... + xp for fixed coefficients.
#' @param randomTV a one-sided formula of the form ~ w1 + w2 + ... + wv | group. w1,.,wv denote the variables associated to the TV random coefficients (1 for the intercept), while group is the indicator for the grouping factor, i.e. the factor identifying the unit longitudinal measurements refer to. Only TC variables are allowed.
#' @param time a string indicating the time variable.
#' @param m number of hidden states associated to the TV random coefficients
#' @param data a data frame containing the variables named in formula, randomTV, and time
#' @param qtl quantile to be estimated
#' @param eps tolerance level for (relative) convergence of the EM algorithm
#' @param maxit maximum number of iterations for the EM algorithm
#' @param se standard error computation
#' @param R number of bootstrap samples for computing standard errors
#' @param start type of starting values (0 = deterministic, 1 = random, 2 = initial values in input)
#' @param parInit list of initial model parameters when \code{start=2}
#' @param verbose if set to FALSE, no printed output is given during the function execution
#' @param ... further arguments to be passed to of from methods
#'
#' @details
#' The function computes ML estimates for the parameters of a linear quantile mixture model, based on TV random coefficients.
#' That is, a linear quantile regression model based on a hidden Markov specification.
#' Estimates are derived by maximizing the (log-)likelihood of a Laplace regression where the location parameter is modeled as a function of fixed and TV random coefficients, as proposed by Farcomeni (2012).
#' The function requires data in long-format and two additional columns indicating the unit identifier and the time occasion.
#' Two formulas specify the model, namely \code{formula},  and \code{formulaTV}:
#' \code{formula} is associated to fixed coefficients; \code{formulaTC} is associated to TV random coefficients.
#' In this latter, only TC variables are allowed.
#'
#' The function admits the presence of missing data, both in terms of drop-outs (monotone missing data)
#' and intermittent missing, under a missing-at-random assumption.
#' Note that, due to the presence of TV random coefficients, intermittent missingness may cause biased inference.
#'
#' If \code{se=TRUE}, standard errors based on a block bootstrap procedure are computed.
#'
#' @import
#' quantreg
#' stats
#' methods
#'
#' @importFrom Rdpack reprompt
#'
#' @return Return an object of \code{\link{class}} \code{lqmix}. This is a list containing the following elements:
#' \item{betaf}{a vector containing fixed regression coefficients}
#' \item{betarTV}{a matrix containing the TV random coefficients}
#' \item{Init}{the initial probability vector of the the hidden Markov chain associated to TV random coefficients}
#' \item{Trans}{the transition probability matrix of the hidden Markov chain associated to TV random coefficients}
#' \item{scale}{the scale parameter}
#' \item{sigma.e}{the standard deviation of error terms}
#' \item{lk}{the log-likelihood}
#' \item{npar}{the total number of model parameters}
#' \item{aic}{the AIC value}
#' \item{bic}{the BIC value}
#' \item{qtl}{the estimated quantile}
#' \item{nsbjs}{the number of subjects}
#' \item{nobs}{the total number of observations}
#' \item{se.betaf}{the standard errors for fixed regression coefficients}
#' \item{se.betarTV}{the standard errors for TV random coefficients}
#' \item{se.scale}{the standard error for the scale parameter}
#' \item{miss}{the missingness type}
#' \item{model}{the estimated model}
#' \item{call}{the matched call}
#'
#' @references{
#'   \insertRef{ref:lqmixTV}{lqmix}
#' }
#' @references{
#'   \insertRef{ref:LM}{lqmix}
#' }
#' @references{
#'   \insertRef{ref:hmm}{lqmix}
#' }
#'
#' @examples
#' outTV = lqmixTV(formula=meas~trt+time+trt:time,randomTV=~1|id,
#'                 time="time",m=2,data=pain,R=10)
#' @export


lqmixTV = function(formula, randomTV, time, m, data, qtl=0.5, eps=10^-5, maxit=1000, se=TRUE, R=50, start=0, parInit = list(betaf=NULL,  betarTV=NULL, delta=NULL, Gamma=NULL, scale=NULL), verbose=TRUE, ...){


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
  if(!inherits(randomTV, "formula") || length(randomTV) != 2) stop("\nTV random coefficient model must be a formula of the form \"~ w|id\".")
  if (qtl <= 0 | qtl >= 1) stop("Quantile level out of range")
  if(m == 1) stop("The specified model corresponds to a linear quantile regression model with no random coefficients. Please, use the function lqr().")


  mformt = as.character(time)
  mform2time = gsub(" ", "", mformt)
  mformTV = strsplit(as.character(randomTV)[2], "\\|")[[1]]
  mformTV = gsub(" ", "", mformTV)
  formula.rTV = reformulate(mformTV[1]) # random formula terms
  mform2sbj = mformTV[2] # id variable

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
  y.obs = data[,as.character(namesY)]
  if(!("opt" %in% names(list(...))))  {
  oo = suppressWarnings(lqmixTV.fit(y=y.obs, x.fixed=x.fixed, namesFix=namesFix, x.random=x.random, namesRan=namesRan,
                   sbj.obs=sbj.obs, time.obs=time.obs, observed=observed, m=m, qtl=qtl, n=n, T=T, Ti=Ti, nObs=nObs,
                   order.time=order.time, ranInt=ranInt, ranSlope=ranSlope, fixed=fixed, start=start, eps=eps,
                   maxit=maxit, parInit=parInit, verbose = verbose))
  }else oo = parInit
  see = se
  if(se){
    if(verbose) cat("Computing standard errors: ")

    # number of parameters
    if(fixed | fixInt) pf = ncol(x.fixed) else pf = 0
    pr = ncol(x.random)

    boot.se = c()
    count = tries = 0
    done = FALSE
    while(count < R & tries <= (R*10)){
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

      if(!is(oo.se, "error")){
        count = count + 1
        if(verbose) cat(count, " ... ")
        if ("search" %in% names(list(...))) cat(count, " ... ")

        boot = list ()
        boot$betaf = oo.se$betaf
        boot$betar = oo.se$betarTV
        boot$delta = oo.se$delta
        boot$Gamma = oo.se$Gamma
        boot$scale = oo.se$scale

        boot.se = rbind(boot.se, unlist(lapply(boot, function(x) c(t(x)))))

      }
    }
    se = apply(boot.se, 2, sd)

    if(fixed | fixInt){
      oo$se.betaf = se[grep("betaf", names(se))]
      names(oo$se.betaf) = namesFix
    }

    oo$se.betarTV = matrix(se[grep("betar", names(se))], nrow = m, byrow = TRUE)
    colnames(oo$se.betarTV) = namesRan
    rownames(oo$se.betarTV) = paste("St", 1:m, sep="")

    oo$se.delta = se[grep("delta", names(se))]
    names(oo$se.delta) = paste("St", 1:m, sep="")
    oo$se.Gamma = matrix(se[grep("Gamma", names(se))], m, m, byrow = T)
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

  oo$call <- match.call()
  class(oo) <- "lqmix"


  set.seed(NULL)
  return(oo)
  message(mess)

}

