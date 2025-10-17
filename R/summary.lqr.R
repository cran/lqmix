#' Summary of an \code{lqr} object
#'
#' Summary method for the \code{\link{class}} \code{\link{lqr}}.
#'
#' @param object an \code{lqr} object
#' @param ... not used
#'
#' @return Return an object of \code{class} \code{summary.lqr}.
#' This is a list of summary statistics for the fitted linear quantile regression model given in \code{object}, with the following elements:
#'
#' \item{fix}{a matrix with estimates, standard errors, Z statistics, and p-values for the regression coefficients}
#' \item{scale}{the scale parameter}
#' \item{sigma.e}{the standard deviation of error terms}
#' \item{lk}{the log-likelihood}
#' \item{npar}{the total number of model parameters}
#' \item{aic}{the AIC value}
#' \item{bic}{the BIC value}
#' \item{qtl}{the estimated quantile}
#' \item{nobs}{the total number of observations}
#' \item{model}{the estimated model}
#' \item{call}{the matched call}
#'
#' @export

summary.lqr <- function(object, ...){

  if(any(!is.null(c(object$se.betaf)))){
    names = c("Estimate", "St.Error", "z.value", "P(>|z|)")

    est = c(object$betaf)
    sef = c(object$se.betaf)
    zvalf = c(object$betaf/sef)
    pvalf = c(1.96*pnorm(-abs(zvalf)))

    tabf = cbind(Estimate = est,
                 St.Err = sef,
                 t.value = zvalf,
                 p.value = pvalf)
    colnames(tabf) = names

    lk = object$lk
    nobs = object$nobs
    model = object$mod
    scale = object$scale
    sigma.e = object$sigma.e
    npar = object$npar
    aic = object$aic
    bic = object$bic
    qtl = object$qtl
    nobs = object$nobs

    res = list()
    res$call = match.call()

    res$fix = tabf
    res$scale = scale
    res$sigma.e = sigma.e
    res$lk = lk
    res$npar = npar
    res$aic = aic
    res$bic = bic
    res$qtl = qtl
    res$nobs = nobs
    res$model = model
    if(!is.null(object$call)) res$call = match.call()

    class(res) = "summary.lqr"
    return(res)
  }else{
    print(object)
    message("Model inference not allowed: standard errors have not been computed.")
  }
}
