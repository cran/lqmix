#' Linear Quantile Mixture with Time-Constant (TC) Random Coefficients
#'
#'  Estimate a finite mixture of linear quantile regression models with Time-Constant (TC), discrete, random coefficients.
#'
#' @param formula an object of class formula of the form y ~ x1 + x2 + ... + xp for fixed coefficients
#' @param randomTC a one-sided formula of the form ~ z1 + z2 + ... + zr | group. z1,.,zr denote the variables associated to the TC random coefficients (1 for the intercept), while group is the indicator for the grouping factor, i.e. the factor identifying the unit longitudinal measurements refer to
#' @param time a string indicating the time variable
#' @param G number of mixture components defining the TC random coefficients
#' @param data a data frame containing the variables named in formula, randomTC, and time
#' @param qtl quantile to be estimated
#' @param eps tolerance level for (relative) convergence of the EM algorithm
#' @param maxit maximum number of iterations for the EM algorithm
#' @param se standard error computation
#' @param R number of bootstrap samples for computing standard errors
#' @param start type of starting values (0 = deterministic, 1 = random, 2 = initial values in input)
#' @param parInit list of initial model parameters when \code{start = 2}
#' @param verbose if set to FALSE, no printed output is given during the function execution
#' @param ... further arguments to be passed to of from methods
#'
#' @details
#' The function computes ML estimates for the parameters of a linear quantile mixture model, based on TC random coefficients.
#' That is, a linear quantile regression model based on a finite mixture specification.
#' Estimates are derived by maximizing the (log-)likelihood of a Laplace regression where the location parameter is modeled as a function of fixed and TC random coefficients, as proposed by Alfo' et. al (2017).
#'
#' The function requires data in long-format and two additional columns indicating the unit identifier and the time occasion.
#' Two formulas specify the model, namely \code{formula} and \code{formulaTC}:
#' \code{formula} is associated to fixed coefficients; \code{formulaTC} is associated to TC random coefficients.
#'
#' The function admits the presence of missing data, both in terms of drop-outs (monotone missing data)
#' and intermittent missing, under a missing-at-random assumption.
#'
#' If \code{se=TRUE}, standard errors based on a block bootstrap procedure are computed.
#'
#' @import
#' quantreg
#' stats
#' methods
#'
#'@importFrom Rdpack reprompt
#'
#' @return Return an object of \code{\link{class}} \code{lqmix}. This is a list containing the following elements:
#' \item{betaf}{a vector containing fixed regression coefficients}
#' \item{betarTC}{a matrix containing the TC random coefficients}
#' \item{pg}{the prior probabilities of the finite mixture associated to TC random coefficients}
#' \item{scale}{the scale parameter}
#' \item{sigma.e}{the standard deviation of error terms}
#' \item{lk}{the log-likelihood}
#' \item{npar}{the total number of model parameters}
#' \item{aic}{the AIC value}
#' \item{bic}{the BIC value}
#' \item{qtl}{the estimated quantile}
#' \item{G}{the number of mixture components}
#' \item{nsbjs}{the number of subjects}
#' \item{nobs}{the total number of observations}
#' \item{se.betaf}{the standard errors for fixed regression coefficients}
#' \item{se.betarTC}{the standard errors for TC random coefficients}
#' \item{se.pg}{the standard errors for prior probabilities of the finite mixture associated to TC random coefficients}
#' \item{se.scale}{the standard error for the scale parameter}
#' \item{miss}{the missingness type}
#' \item{model}{the estimated model}
#' \item{call}{the matched call}
#'
#' @references{
#'   \insertRef{ref:lqmixTC}{lqmix}
#' }
#' @references{
#'   \insertRef{ref:npml1}{lqmix}
#' }
#' @references{
#'   \insertRef{ref:npml2}{lqmix}
#' }
#' @examples
#' outTC = lqmixTC(formula=meas~trt+time+trt:time,randomTC=~1|id,
#'                 time="time",G=2,data=pain,se=TRUE,R=10)
#' @export


lqmixTC = function(formula, randomTC, time, G, data, qtl=0.5, eps=10^-5, maxit=1000, se=TRUE, R=50, start=0, parInit = list(betaf=NULL, betarTC=NULL, pg=NULL, scale=NULL), verbose=TRUE, ...){

  # start = 0 -- deterministic start
  # start = 1 -- random start
  # start = 2 -- start from given values


  # ---- possible errors ------
  # ***************************
  if(se==TRUE & is.null(R)){
    mess = "\n The number of bootstrap samples for computing standard errors has not been specified. The default value R=50 has been used."
  }else mess=NULL

  if(start == 2 & is.null(unlist(parInit))) stop("No input parameters have been given with start = 2.")

  if(!is.data.frame(data))  stop("`data' must be a data frame.")
  if(!inherits(formula, "formula") || length(formula) != 3) stop("\nFixed coefficient model must be a formula of the form \"y ~ x\".")
  if (!inherits(randomTC, "formula") || length(randomTC) != 2) stop("\nTC random coefficient model must be a formula of the form \"~ z|id\".")
  if (qtl <= 0 | qtl >= 1) stop("Quantile level out of range")

  if(G == 1) stop("The specified model corresponds to a linear quantile regression model with no random coefficients. Please, use the function lqr().")
  mformt = as.character(time)
  mform2time = gsub(" ", "", mformt)
  mformTC = strsplit(as.character(randomTC)[2], "\\|")[[1]]
  mformTC = gsub(" ", "", mformTC)
  formula.rTC = reformulate(mformTC[1]) # random formula terms
  mform2sbj = mformTC[2] # id variable

  if (!(all(mform2sbj %in% names(data)))) stop("The specified clustering variable is not contained in the data frame.")
  if (!(all(mform2time %in% names(data)))) stop("The specified time variable is not contained in the data frame.")

  randomTerms = attr(terms(formula.rTC), "intercept") + length(attr(terms(formula.rTC),"term.labels"))>0
  if(!randomTerms) stop("The specifid model corresponds to a linear quantile regression model with no random coefficients. Please, use the function lqr().")

  # ---- initial settings ----
  # ***************************
  # printing options
  if(verbose & length(as.logical(list(...)))==0) cat("Model TC - qtl =", qtl, "\n")

  # remove incomplete data
  names = c(unique(unlist(lapply(c(formula, formula.rTC), all.vars))), mform2sbj, mform2time)
  asOneFormula = eval(parse(text = paste("~", paste(names, collapse = "+")))[[1]])
  data = model.frame(asOneFormula, data)

  # response
  namesY = formula[[2]]

  # define fixed and random model frames and matrices
  # fixed model frame
  mff = model.frame(formula, data)
  # random model frame
  mfr = model.frame(formula.rTC, data)
  # random model matrix
  mmr = model.matrix(formula.rTC, mfr)


  # intercept derived from the formula
  fixInt = attr(terms(mff), "intercept") == 1
  ranInt = attr(terms(mfr), "intercept") == 1
  # if random intercept, remove intercept from fixed formula
  if(ranInt == 1 & fixInt == 1) fixInt = 0

  # identify the type of intercept: 0 -> fixed, 1 -> random TC, 999 -> no intercept in the model
  ranInt = ifelse(ranInt == 1 & fixInt == 0, 1, ifelse(ranInt == 0 & fixInt == 1, 0, 999))

  # variable names
  termsFix = attr(terms(formula),"term.labels")
  termsRan = attr(terms(formula.rTC),"term.labels")
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


  # ordering the dataset wrt subjects
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
  # n*T objects
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
  mfr = model.frame(formula.rTC, data)
  x.random = model.matrix(formula.rTC, mfr)

  # fixed covariates
  if(!is.null(formula)){
    mff = model.frame(formula, data)
    x.fixed = model.matrix(formula, mff)
  }else x.fixed = NULL

  # response variable
  y.obs = data[,as.character(namesY)]

  if(!("opt" %in% names(list(...)))) {
  oo = lqmixTC.fit(y=y.obs, x.fixed=x.fixed, namesFix=namesFix, x.random=x.random, namesRan=namesRan,
                   id=sbj.obs, G=G, qtl=qtl, n=n, T=T, Ti=Ti, nObs=nObs,
                   order.time=order.time, ranInt=ranInt, ranSlope=ranSlope, fixed=fixed, start=start, eps=eps,
                   maxit=maxit, parInit=parInit,verbose=verbose)
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

      # build the complete data matrices of size (n*T)*(?)
      xx.fixed = matrix(, n*T, pf)
      xx.random = matrix(, n*T, pr)
      yy = c(matrix(, n*T, 1))

      xx.fixed[observed,] = x.fixed
      xx.random[observed,] = x.random
      yy[observed] = y.obs

      # build new data matrices based on sampled units
      set.seed(tries)
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
      if(!is(oo.se, "error")){
        count = count + 1
        if(verbose) cat(count, " ... ")
        if ("search" %in% names(list(...))) cat(count, " ... ")

        boot = list ()
        boot$betaf = oo.se$betaf
        boot$betar = oo.se$betar
        boot$pg = oo.se$pg
        boot$scale = oo.se$scale
        boot$sigma.b = oo.se$sigma.b

        boot.se = rbind(boot.se, unlist(lapply(boot, function(x) c(t(x)))))

      }
    }

    se = apply(boot.se, 2, sd)

    if(fixed | fixInt){
      oo$se.betaf = se[grep("betaf", names(se))]
      names(oo$se.betaf) = namesFix
    }
    oo$se.betarTC = matrix(se[grep("betar", names(se))], nrow = G, byrow = TRUE)
    colnames(oo$se.betarTC) = namesRan
    rownames(oo$se.betarTC) = paste("Comp", 1:G, sep="")

    oo$se.pg = se[grep("pg", names(se))]
    names(oo$se.pg) = paste("Comp", 1:G, sep="")
    oo$se.scale = se[grep("scale", names(se))]
    cat("\n")
  }


  if(any(Ti != T)){
    if((any(time.input != timeidx))){
      oo$miss = "non-monotone"
    }else oo$miss = "monotone"
  }else oo$miss = "none"

  oo$pf = oo$pr = NULL
  oo$model = "TC"

  oo$call <- match.call()
  class(oo) <- "lqmix"

  set.seed(NULL)
  return(oo)

  message(mess)

}


