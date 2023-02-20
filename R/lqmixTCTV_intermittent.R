#' Linear Quantile Mixture with Time-Constant (TC) and Time-Varuying (TV) Random Coefficients
#'
#' Estimate a finite mixture of linear quantile regression models with both Time-Constant (TC) and Time-Varying (TV), discrete, random coefficients.
#'
#' @param formula an object of class formula of the form y ~ x1 + x2 + ... + xp for fixed coefficients
#' @param randomTC a one-sided formula of the form ~ z1 + z2 + ... + zr | group. z1,.,zr denote the variables associated to TC random coefficients (1 for the intercept), while group is the indicator for the grouping factor, i.e. the factor identifying the unit longitudinal measurements refer to
#' @param randomTV a one-sided formula of the form ~ w1 + w2 + ... + wl | group. w1,.,wl denote the variables associated to TV random coefficients (1 for the intercept), while group is the indicator for the grouping factor, i.e. the factor identifying the unit longitudinal measurements refer to. Only TC variables are allowed
#' @param time a string indicating the time variable
#' @param m number of hidden states associated the TV random coefficients
#' @param G number of mixture components associated the TC random coefficients
#' @param data a data frame containing the variables named in formula, randomTC, randomTV, and time
#' @param qtl quantile to be estimated
#' @param eps tolerance level for (relative) convergence of the EM algorithm
#' @param maxit maximum number of iterations for the EM algorithm
#' @param se standard error computation
#' @param R number of bootstrap samples for computing standard errors
#' @param start type of starting values (0 = deterministic, 1 = random, 2 = initial values in input)
#' @param start type of starting values (0 = deterministic, 1 = random, 2 = initial values in input)
#' @param parInit list of initial model parameters when \code{start=2}
#' @param verbose if set to FALSE, no printed output is given during the function execution
#' @param ... further arguments to be passed to of from methods
#'
#' @details
#' The function computes ML estimates for the parameters of a linear quantile mixture model, based on TC and TV random coefficients.
#' That is, a linear quantile regression model based on a mixed hidden Markov specification.
#' Estimates are derived by maximizing the (log-)likelihood of a Laplace regression where the location parameter is modeled as a function of fixed, TC, and TV random coefficients, as proposed by Marino et. al (2018).
#'
#' The function requires data in long-format and two additional columns indicating the unit identifier and the time occasion.
#' Three formulas specify the model, namely \code{formula}, \code{formulaTC}, and \code{formulaTV}:
#' \code{formula} is associated to fixed coefficients; \code{formulaTC} is associated to TC random coefficients; \code{formulaTC} is associated to TV random coefficients.
#' In this latter, only TC variables are allowed.
#'
#' The function admits the presence of missing data, both in terms of drop-outs (monotone missing data)
#' and intermittent missing, under a missing-at-random assumption.
#' Note that, due to the presence of TV random coefficients, intermittent missingness may cause biased inference.
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
#' \item{betarTV}{a matrix containing the TV random coefficients}
#' \item{pg}{the prior probabilities of the finite mixture associated to TC random coefficients}
#' \item{delta}{the initial probability vector of the hidden Markov chain associated to TV random coefficients}
#' \item{Gamma}{the transition probability matrix of the hidden Markov chain associated to TV random coefficients}
#' \item{scale}{the scale parameter}
#' \item{sigma.e}{the standard deviation of error terms}
#' \item{lk}{the log-likelihood}
#' \item{npar}{the total number of model parameters}
#' \item{aic}{the AIC value}
#' \item{bic}{the BIC value}
#' \item{qtl}{the estimated quantile}
#' \item{m}{the number of hidden states}
#' \item{G}{the number of mixture components}
#' \item{nsbjs}{the number of subjects}
#' \item{nobs}{the total number of observations}
#' \item{se.betaf}{the standard errors for fixed regression coefficients}
#' \item{se.betarTC}{the standard errors for TC random coefficients}
#' \item{se.betarTV}{standard errors for TV random coefficients}
#' \item{se.Init}{the standard errors for the initial probabilities of the hidden Markov chain associated to TV random coefficients}
#' \item{se.Trans}{the standard errors for the transition probabilities of the hidden Markov chain associated to TV random coefficients}
#' \item{se.Mprob}{the standard errors for the prior probabilities of the finite mixture associated to TC random coefficients}
#' \item{se.scale}{the standard error for the scale parameter}
#' \item{miss}{the missingness type}
#' \item{model}{the estimated model}
#' \item{call}{the matched call}
#'
#' @references{
#'   \insertRef{ref:lqmixTCTV}{lqmix}
#' }
#' @references{
#'   \insertRef{ref:mhmm1}{lqmix}
#' }
#' @references{
#'   \insertRef{ref:mhmm2}{lqmix}
#' }
#'
#' @examples
#' outTCTV = lqmixTCTV(formula=meas~trt+time+trt:time,randomTC=~time|id,
#'                     randomTV=~1|id,time="time",m=2,G=2,data=pain,R=10)
#' @export
#'

lqmixTCTV = function(formula, randomTC, randomTV, time, m, G, data, qtl=0.5, eps=10^-5, maxit=1000, se=TRUE, R=50, start=0, parInit = list(betaf=NULL,  betarTC=NULL, betarTV=NULL, pg=NULL, delta = NULL, Gamma = NULL, scale=NULL),verbose=TRUE, ...){

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
  if (!inherits(randomTV, "formula") || length(randomTC) != 2) stop("\nTV random coefficient model must be a formula of the form \"~ w|id\".")
  if (qtl <= 0 | qtl >= 1) stop("Quantile level out of range")



  # ---- initial settings ----
  # ***************************
  if(!is.data.frame(data))  stop("`data' must be a data frame.")

  if(!inherits(formula, "formula") || length(formula) != 3) stop("\nFixed coefficient model must be a formula of the form \"y ~ x\".")
  if(length(unlist(sapply(c("~", "\\|"), grep, gsub(" ", "", as.character(randomTV))))) != 2) stop("\nTV random coefficient model must be a formula of the form \"~ w|id\".")
  if(length(unlist(sapply(c("~", "\\|"), grep, gsub(" ", "", as.character(randomTC))))) != 2) stop("\nTC random coefficient model must be a formula of the form \"~ z|id\".")

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
  mform2sbjTC = mformTC[2] # id variable

  mformTV = strsplit(as.character(randomTV)[2], "\\|")[[1]]
  mformTV = gsub(" ", "", mformTV)
  formula.rTV = reformulate(mformTV[1], intercept = ranIntTV) # random formula terms
  mform2sbjTV = mformTV[2] # id variable

  if(mform2sbjTV != mform2sbjTC) stop("The specified clustering variable differs in the TC and TV random coefficient model formulas.")

  if (!(all(mform2sbjTC %in% names(data)))) stop("The specified clustering variable is not contained in the data frame.")
  if (!(all(mform2time %in% names(data)))) stop("The specified time variable is not contained in the data frame.")

  randomTermsTC = attr(terms(formula.rTC), "intercept") + length(attr(terms(formula.rTC),"term.labels"))>0
  randomTermsTV = attr(terms(formula.rTV), "intercept") + length(attr(terms(formula.rTV),"term.labels"))>0
  if(!(randomTermsTC & randomTermsTV)) stop("The specifid model corresponds to a linear quantile regression model with no random coefficients. Please, use the function lqr().")
  else if(!(randomTermsTC) & randomTermsTV) stop("The specifid model corresponds to a linear quantile regression model with TV random coefficients only. Please, use the function lqmixTV().")
  else if(randomTermsTC & !(randomTermsTV)) stop("The specifid model corresponds to a linear quantile regression model with TC random coefficients only. Please, use the function lqmixTC().")


  # ---- initial settings ----
  # ***************************
  # printing options
  if(verbose  & length(as.logical(list(...)))==0) cat("Model TCTV - qtl =", qtl, "\n")

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

  # if random intercept, remove intercept from fixed formula
  if((ranIntTC == 1 | ranIntTV == 1) & fixInt == 1) fixInt = 0

  # identify the type of intercept: 0 -> fixed, 1 -> random TC, 2 -> TV, 999 -> no intercept in the model
  ranInt = ifelse(ranIntTC == 1 & fixInt == 0, 1,
                  ifelse(ranIntTV == 1 & fixInt == 0, 2,
                  ifelse(ranIntTC == 0 & ranIntTV == 0 & fixInt == 1, 0, 999)))

  # variable names
  termsFix = attr(terms(formula),"term.labels")
  termsRanTC = attr(terms(formula.rTC),"term.labels")
  termsRanTV = attr(terms(formula.rTV),"term.labels")
  # ATT: intercepts do not appear in the list

  # slopes in common between random and fixed model formula
  wh = which(termsFix %in% c(termsRanTC, termsRanTV))
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
  ranSlopeTC = FALSE
  if(length(termsRanTC>0)) ranSlopeTC = 1

  ranSlopeTV = FALSE
  if(length(termsRanTV>0)) ranSlopeTV = 1

  # namesFix and namesRan
  namesFix = colnames(mmf)
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
  # *******************************
  # random covariates
  mfrTC = model.frame(formula.rTC, data)
  x.randomTC = model.matrix(formula.rTC, mfrTC)

  mfrTV = model.frame(formula.rTV, data)
  x.randomTV = model.matrix(formula.rTV, mfrTV)

  # fixed covariates
  if(!is.null(formula)){
    mff = model.frame(formula, data)
    x.fixed = model.matrix(formula, mff)
  }else x.fixed = NULL

  # response variable
  y.obs = data[,as.character(namesY)]

  if(ranSlopeTV){
    TVcovar = function(cov, sbj){
      tapply(cov, sbj, function(xx){length(table(xx))>1})
    }
    if(ranInt == 1) checkTVcovar = colSums(apply(as.matrix(x.randomTV[,-1], nrow = nObs), 2, function(xx){TVcovar(xx, sbj.obs)}))>0
    else checkTVcovar = colSums(apply(as.matrix(x.randomTV, nrow = nObs), 2, function(xx){TVcovar(xx, sbj.obs)}))>0
    if(any(checkTVcovar)) stop(paste("TV random coefficients are only admitted for TC covariates."))
  }
  if(!("opt" %in% names(list(...)))) {
    oo = lqmixTCTV.fit(y=y.obs, x.fixed=x.fixed, namesFix=namesFix,
                       x.randomTC=x.randomTC, namesRanTC=namesRanTC,
                       x.randomTV=x.randomTV, namesRanTV=namesRanTV,
                       sbj.obs=sbj.obs, time.obs=time.obs, observed=observed, m=m, G=G, qtl=qtl, n=n, T=T, Ti=Ti, nObs=nObs,
                       order.time=order.time, ranInt=ranInt, ranSlopeTC=ranSlopeTC,ranSlopeTV=ranSlopeTV, fixed=fixed, start=start, eps=eps,
                       maxit=maxit, parInit=parInit,verbose=verbose)
    } else oo = parInit
    see = se
    if(se){

      if(verbose) cat("Computing standard errors: ")

      # number of parameters
      if(fixed | fixInt) pf = ncol(x.fixed) else pf = 0
      prTC = ncol(x.randomTC)
      prTV = ncol(x.randomTV)

      boot.se = c()
      count = tries = 0
      while(count < R & tries <= (R*10)){

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
        set.seed(tries)
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

        if(!is(oo.se, "error")){
          count = count + 1
          if(verbose) cat(count, " ... ")
          if ("search" %in% names(list(...))) cat(count, " ... ")
          boot = list ()
          boot$betaf = oo.se$betaf
          boot$betarTC = oo.se$betarTC
          boot$betarTV = oo.se$betarTV
          boot$delta = oo.se$delta
          boot$Gamma = oo.se$Gamma
          boot$pg = oo.se$pg
          boot$scale = oo.se$scale

          boot.se = rbind(boot.se, unlist(lapply(boot, function(x) c(t(x)))))

        }
      }

      se = apply(boot.se, 2, sd)

      if(fixed | fixInt){
        oo$se.betaf = se[grep("betaf", names(se))]
        names(oo$se.betaf) = namesFix
      }

      oo$se.betarTC = matrix(se[grep("betarTC", names(se))], nrow = G, byrow = TRUE)
      colnames(oo$se.betarTC) = namesRanTC
      rownames(oo$se.betarTC) = paste("Comp", 1:G, sep="")

      oo$se.betarTV = matrix(se[grep("betarTV", names(se))], nrow = m, byrow = TRUE)
      colnames(oo$se.betarTV) = namesRanTV
      rownames(oo$se.betarTV) = paste("St", 1:m, sep="")

      oo$se.delta = se[grep("delta", names(se))]
      names(oo$se.delta) = paste("St", 1:m, sep="")
      oo$se.Gamma = matrix(se[grep("Gamma", names(se))], m, m, byrow = T)
      rownames(oo$se.Gamma) = paste("fromSt", 1:m, sep="")
      colnames(oo$se.Gamma) = paste("toSt", 1:m, sep="")
      oo$se.pg = se[grep("pg", names(se))]
      names(oo$se.pg) = paste("Comp", 1:G, sep="")
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

    oo$prTC = oo$prTC = NULL
    oo$prTV = oo$prTC = NULL
    oo$pf = oo$pr = NULL
    oo$model ="TCTV"

    oo$call = match.call()
    class(oo) = "lqmix"
    set.seed(NULL)
    return(oo)

  message(mess)

}

