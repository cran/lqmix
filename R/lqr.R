#' Linear Quantile Regression
#'
#' Estimate a linear quantile regression model for independent data (no random coefficients).
#'
#' @param formula an object of class \code{formula}: a symbolic description of the model to be fitted
#' @param data a data frame containing the variables named in \code{formula}
#' @param qtl quantile to be estimated
#' @param se standard error computation
#' @param R number of bootstrap samples for computing standard errors
#' @param verbose if set to FALSE, no printed output is given during the function execution
#' @param seed an integer value for random numbers generation, used for bootstrap standard errors
#' @param parallel if set to TRUE, a parallelized code is use for standard error computation (if se=TRUE)
#' @param ncores number of cores used for computing bootstrap standard errors (if required)
#' @param ... not used
#'
#' @details
#' The function computes ML estimates for the parameters of a linear quantile regression model for independent observations.
#' Estimates are derived by maximizing the (log-)likelihood of a Laplace regression, where the location parameter is modeled as a function
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
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom utils setTxtProgressBar
#' @importFrom utils txtProgressBar
#' @importFrom doSNOW registerDoSNOW
#' @importFrom foreach foreach
#' @importFrom doParallel stopImplicitCluster
#' @importFrom doRNG %dorng%
#'
#' @return Return an object of \code{class} \code{lqr}. This is a list containing the following elements:
#' \item{betaf}{a vector containing fixed regression coefficients}
#' \item{scale}{the scale parameter}
#' \item{sigma.e}{the standard deviation of error terms}
#' \item{lk}{the log-likelihood}
#' \item{npar}{the total number of model parameters}
#' \item{AIC}{the AIC value}
#' \item{BIC}{the BIC value}
#' \item{qtl}{the estimated quantile}
#' \item{nobs}{the total number of observations}
#' \item{se.betaf}{the standard errors for the regression coefficients}
#' \item{se.scale}{the standard error for the scale parameter}
#' \item{model}{the estimated model}
#' \item{mmf}{the model matrix associated to the regression coefficients}
#' \item{y}{the model response}
#' \item{call}{the matched call}
#' \item{formula}{the model formula}
#'
#' @references{
#'   \insertRef{ref:KoeBas}{lqmix}
#' }
#'
#' @examples
#' out0 = lqr(formula=meas~trt+time+trt:time,data=pain,se=TRUE,R=10)
#' @export

lqr = function(formula, data, qtl=0.5, se=TRUE, R=100, verbose=TRUE, seed=NULL, parallel=FALSE, ncores=2, ...){

  # ---- possible errors -----
  # ***************************
  if(is.null(data)) stop("No input dataset has been given.")
  if(!is.data.frame(data))  stop("`data' must be a data frame.")
  if(!inherits(formula, "formula") || length(formula) != 3) stop("\nFixed coefficient model must be a formula of the form \"y ~ x\".")

  if (qtl <= 0 | qtl >= 1) stop("Quantile level out of range")

  opt = ("opt" %in% names(list(...)))

  # remove incomplete data
  names = c(unique(unlist(lapply(c(formula), all.vars))))
  asOneFormula = eval(parse(text = paste("~", paste(names, collapse = "+")))[[1]])
  data = model.frame(asOneFormula, data)

  xy = xyBuildHOM(formula=formula, data=data)

  y.obs = xy$y.obs; x.fixed = xy$x.fixed
  namesFix = xy$nameFix
  nObs = dim(x.fixed)[1]

  oo = lqr.fit(y=y.obs,x.fixed=x.fixed,namesFix=namesFix,qtl=qtl,nObs=nObs,verbose=verbose)

  see = se
  if(se){
    if(verbose & !opt) cat("Computing standard errors ...\n") else if(verbose & opt) cat("Computing standard errors for the optimal model...\n")
    if(parallel==TRUE) cl = makeCluster(ncores) else cl = makeCluster(1)
    cc = 0
    registerDoSNOW(cl)
    if(verbose){
      pb = txtProgressBar(max = R, style = 3)
      progress = function(x) setTxtProgressBar(pb,x)
      opts = list(progress = progress)
    }else opts = list()

    ooo.se = foreach(cc = (1 : R), .options.snow = opts) %dorng% {
      tries = 0

      while(tries <= (R*10)){
        tries = tries + 1

        if (tries == R*10)  stop("Standard errors may not be computed.")

        # build the complete data matrices of size (n*T)*(?)
        pf = oo$pf

        # build new data matrices based on sampled units
        if(!is.null(seed)) set.seed(seed*tries*cc)# else set.seed(tries*cc)

        sample.unit = sample(1:nObs, nObs, replace=TRUE)

        x.fix.sample = matrix(x.fixed[sample.unit,], nObs)
        y.sample = y.obs[sample.unit]

        colnames(x.fix.sample) = namesFix

        oo.se = tryCatch(suppressWarnings(lqr.fit(y=y.sample, x.fixed=x.fix.sample,
                                 namesFix=namesFix, qtl=qtl, nObs=nObs,verbose=FALSE)),error=function(e){e})
        if(!is(oo.se, "error")) break
      }


      boot = list ()
      boot$betaf = oo.se$betaf
      boot$scale = oo.se$scale

      return(boot)
    }

    if(verbose) close(pb)
    stopImplicitCluster()
    parallel::stopCluster(cl)
    rm(cl)

    varcov = cov(t(sapply(ooo.se,unlist)))
    se = sqrt(diag(varcov))

    wh = grep("betaf", names(se))
    oo$se.betaf = se[wh]
    varcov = as.matrix(varcov[wh, wh])
    names(oo$se.betaf) = colnames(varcov) = rownames(varcov) = namesFix

    oo$se.scale = se[grep("scale", names(se))]
    oo$vcov = varcov
  }

  if(!inherits(oo, "error")){
    oo$pf = oo$pr = NULL
    oo$model = "Homogeneous"

    oo$call <- match.call()
    oo$formula <- formula
    oo$mmf = x.fixed
    oo$y = y.obs
    oo$formula = formula
    class(oo) <- "lqr"

  }
  return(oo)
}

