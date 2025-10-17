#' Predictions from an \code{lqr} object
#'
#' Returns the predicted values for an object of \code{\link{class}} \code{\link{lqr}}.
#'
#' @param object an \code{lqr} object
#' @param newdata an optional data frame in which to look for variables with which to predict. If omitted, the fitted values are produced
#' @param \dots not used
#'
#' @details
#' The function computes predictions for an object of \code{\link{class}} \code{\link{lqr}}.
#'
#' @return A vector of predictions.
#'
#' @export

predict.lqr = function (object, newdata = NULL, ...){


  qtl = object$qtl
  formula = object$formula
  betaf = object$betaf

  if(!is.null(newdata)){

    if (!inherits(newdata, "data.frame"))
      stop("'newdata' must be a data frame")
    if (!all(all.vars(formula) %in% names(newdata)))
      stop("newdata must have all terms in the formula from main call")

    nObs = dim(newdata)[1]
    rown = rownames(newdata)
    # build design matrices
    xy = xyBuildHOM(formula=formula, data=newdata)
    y.obs = xy$y.obs; x.fixed = xy$x.fixed

  }else{
    nObs = object$nobs
    x.fixed = object$mmf
    rown = rownames(object$mmf)
  }
  linear.predictor = c(x.fixed %*% betaf)

  names(linear.predictor) = rown
  return(linear.predictor)
}
