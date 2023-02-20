#' Linear Quantile Regression
#'
#' Estimate a linear quantile regression model with no random coefficients
#'
#' @param formula an object of class formula of the form y ~ x1 + x2 + ... + xp for fixed coefficients
#' @param data a data frame containing the variables named in formula and time
#' @param qtl quantile to be estimated
#' @param se standard error computation
#' @param R number of bootstrap sample for computing standard errors
#' @param verbose if set to FALSE, no printed output is given during the function execution
#' @param ... further arguments to be passed to of from methods
#'
#' @details
#' The function computes ML estimates for the parameters of a linear quantile regression model for independent observations
#' Estimates are derived by maximizing the (log-)likelihood of a Laplace regression where the location parameter is modeled as a function
#' of fixed coefficients only.
#'
#' If \code{se=TRUE}, standard errors based on a bootstrap procedure are computed.
#'
#' @import
#' quantreg
#' stats
#' methods
#'
#' @importFrom Rdpack reprompt
#'
#' @return Return an object of \code{\link{class}} \code{lqr}. This is a list containing the following elements:
#' \item{betaf}{a vector containing fixed regression coefficients}
#' \item{scale}{the scale parameter}
#' \item{sigma.e}{the standard deviation of error terms}
#' \item{lk}{the log-likelihood}
#' \item{npar}{the total number of model parameters}
#' \item{aic}{the AIC value}
#' \item{bic}{the BIC value}
#' \item{qtl}{the estimated quantile}
#' \item{nobs}{the total number of observations}
#' \item{se.betaf}{the standard errors for fixed regression coefficients}
#' \item{se.scale}{the standard error for the scale parameter}
#' \item{model}{the estimated model}
#' \item{call}{the matched call}
#'
#' @references{
#'   \insertRef{ref:lqr}{lqmix}
#' }
#'
#' @examples
#' out0 = lqr(formula=meas~trt+time+trt:time,data=pain,se=TRUE,R=10)
#' @export

lqr = function(formula, data, qtl=0.5, se=TRUE, R=50, verbose=TRUE, ...){

  # ---- possible errors -----
  # ***************************
  if(se==TRUE & is.null(R)){
    mess = "\n The number of bootstrap samples for computing standard errors has not been specified. The default value R=50 has been used."
  }else mess=NULL

  if(!is.data.frame(data))  stop("`data' must be a data frame.")
  if(!inherits(formula, "formula") || length(formula) != 3) stop("\nFixed coefficient model must be a formula of the form \"y ~ x\".")

  if (qtl <= 0 | qtl >= 1) stop("Quantile level out of range")

  if(verbose & length(as.logical(list(...)))==0) cat("Model homogeneous", "- qtl =", qtl,"\n")

  # remove incomplete data
  names = c(unique(unlist(lapply(c(formula), all.vars))))
  asOneFormula = eval(parse(text = paste("~", paste(names, collapse = "+")))[[1]])
  data = model.frame(asOneFormula, data)

  # response
  namesY = as.character(formula)[[2]]

  # fixed model frame
  mff = model.frame(formula, data)
  mmf = model.matrix(formula, mff)

  # variable names
  termsFix = attr(terms(formula),"term.labels")
  # namesFix
  namesFix = colnames(mmf)

  # fixed covariates
  x.fixed = model.matrix(formula, mff)
  # response variable
  y.obs = data[,as.character(namesY)]
  nObs = dim(x.fixed)[1]

  oo = lqr.fit(y=y.obs,x.fixed=x.fixed,namesFix=namesFix,qtl=qtl,nObs=nObs,verbose=verbose)

  see = se
  if(se){
    if(verbose) cat("Computing standard errors: ")
    boot.se = c()
    count = tries = 0
    done = FALSE

    while(count < R & tries <= (R*10)){
      tries = tries + 1
      if (tries == R*10)  stop("Standard errors may not be computed.")

      # build the complete data matrices of size (n*T)*(?)
      pf = oo$pf

      # build new data matrices based on sampled units
      set.seed(tries)
      sample.unit = sample(1:nObs, nObs, replace=TRUE)

      x.fix.sample = x.fixed[sample.unit,]
      y.sample = y.obs[sample.unit]

      colnames(x.fix.sample) = namesFix

      oo.se = tryCatch(suppressWarnings(lqr.fit(y=y.sample, x.fixed=x.fix.sample,
                               namesFix=namesFix, qtl=qtl, nObs=nObs,verbose=FALSE)),error=function(e){e})
      if(!is(oo.se, "error")){
        count = count + 1
        if(verbose==T) cat(count, " ... ")
        if (length(as.logical(list(...)))>0) cat(count, " ... ")

        boot = list ()
        boot$betaf = oo.se$betaf
        boot$scale = oo.se$scale

        boot.se = rbind(boot.se, unlist(lapply(boot, function(x) c(t(x)))))

      }
    }

    se = apply(boot.se, 2, sd)

    oo$se.betaf = se[grep("betaf", names(se))]
    names(oo$se.betaf) = namesFix
    oo$se.scale = se[grep("scale", names(se))]


  }

  oo$pf = oo$pr = NULL
  oo$model = "Homogeneous"
  oo$call = match.call()
  class(oo) <- "lqr"

  return(oo)


}

