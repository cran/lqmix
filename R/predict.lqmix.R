#' Predictions from an \code{lqmix} object
#'
#' Returns the predicted values for an object of class \code{\link{lqmix}}.
#'
#' @param object an \code{lqmix} object
#' @param newdata an optional data frame in which to look for variables with which to predict. If omitted, the fitted values are produced
#' @param \dots not used
#'
#' @details
#' The function computes predictions for an object of \code{\link{class}} \code{\link{lqmix}}.
#' If the fitted model is based on TC, discrete, random coefficients only, a matrix of size \code{nsbjs x G} is given as output; if the fitted model is based on TV, discrete, random coefficients only, a matrix of size \code{nobs x m}; if the fitted model is based on both TC and TV, discrete, random coefficients, an array of size \code{nobs x G x m} is returned.
#
#' @return A matrix or an array of predictions, based on the estimated model.
#'
#' @export

predict.lqmix = function (object, newdata = NULL, ...){


  qtl = object$qtl
  formula = object$formula
  randomTC = object$randomTC
  randomTV = object$randomTV

  group = object$group
  time = object$time

  mod = object$mod
  betaf = object$betaf
  betarTC = object$betarTC
  betarTV = object$betarTV


  G = object$G
  m = object$m

  if(!is.null(newdata)){

    if (!inherits(newdata, "data.frame"))
      stop("'newdata' must be a data frame")
    if (!all(all.vars(formula) %in% names(newdata)))
      stop("newdata must have all terms in 'fixed' formula from main call")
    if (!all(all.vars(randomTC) %in% names(newdata)))
      stop("newdata must have all terms in 'randomTC' formula from main call")
    if (!all(all.vars(randomTV) %in% names(newdata)))
      stop("newdata must have all terms in 'randomTV' formula from main call")
    if(!object$group %in% names(newdata))
      stop("newdata must contain the grouping variable")
    if(!object$time %in% names(newdata))
      stop("newdata must contain the time variable")

    nObs = dim(newdata)[1]
    rown = rownames(newdata)

    # build design matrices
    if(mod == "TC"){
      ranIntTC = length(grep("1", randomTC))>0
      mformTC = strsplit(as.character(randomTC)[2], "\\|")[[1]]
      mformTC = gsub(" ", "", mformTC)
      formula.rTC = reformulate(mformTC[1], intercept = ranIntTC) # random formula terms

      xy = xyBuildTC(formula=formula, randomTC=formula.rTC, data=newdata)
      y.obs = xy$y.obs; x.fixed = xy$x.fixed; x.randomTC = xy$x.random


      if(!is.null(betaf)) Xbeta = matrix(x.fixed %*% betaf, nObs, G) else Xbeta=0
      Zbg = x.randomTC %*% t(betarTC)
      linear.predictor = Xbeta +  Zbg

    }else if(mod == "TV"){

      ranIntTV = length(grep("1", randomTV))>0
      mformTV = strsplit(as.character(randomTV)[2], "\\|")[[1]]
      mformTV = gsub(" ", "", mformTV)
      formula.rTV = reformulate(mformTV[1], intercept = ranIntTV) # random formula terms

      xy = xyBuildTV(formula=formula, randomTV=formula.rTV, data=newdata)
      y.obs = xy$y.obs; x.fixed = xy$x.fixed; x.randomTV = xy$x.random

      if(!is.null(betaf)) Xbeta = matrix(x.fixed %*% betaf, nObs, m) else Xbeta=0
      Wbetah = x.randomTV %*% t(betarTV)
      linear.predictor = Xbeta +  Wbetah


    }else if(mod == "TCTV"){
      ranIntTC = length(grep("1", randomTC))>0
      mformTC = strsplit(as.character(randomTC)[2], "\\|")[[1]]
      mformTC = gsub(" ", "", mformTC)
      formula.rTC = reformulate(mformTC[1], intercept = ranIntTC) # random formula terms

      ranIntTV = length(grep("1", randomTV))>0
      mformTV = strsplit(as.character(randomTV)[2], "\\|")[[1]]
      mformTV = gsub(" ", "", mformTV)
      formula.rTV = reformulate(mformTV[1], intercept = ranIntTV) # random formula terms

      xy = xyBuildTCTV(formula=formula, randomTC=formula.rTC, randomTV=formula.rTV, data=newdata)
      y.obs = xy$y.obs; x.fixed = xy$x.fixed; x.randomTC = xy$x.randomTC; x.randomTV = xy$x.randomTV


      if(!is.null(betaf)) Xbeta = array(x.fixed %*% betaf, c(nObs, m, G))  else Xbeta = 0
      Zbg = aperm(array(x.randomTC %*% t(betarTC), c(nObs, G, m)), c(1,3,2))
      Wbetah = array(x.randomTV %*% t(betarTV), c(nObs, m, G))
      linear.predictor = Xbeta + Zbg + Wbetah
      dimnames(linear.predictor)[2:3] = list(names(object$pg), names(object$delta))
    }


  }else{
    nObs = object$nobs
    x.fixed = object$mmf
    x.randomTC = object$mmrTC
    x.randomTV = object$mmrTV
    rown = rownames(object$mmf)

    if(mod == "TC"){

      if(!is.null(betaf)) Xbeta = matrix(x.fixed %*% betaf, nObs, G) else Xbeta=0
      Zbg = x.randomTC %*% t(betarTC)
      linear.predictor = Xbeta +  Zbg

    }else if(mod == "TV"){
      if(!is.null(betaf)) Xbeta = matrix(x.fixed %*% betaf, nObs, m) else Xbeta=0
      Wbetah = x.randomTV %*% t(betarTV)
      linear.predictor = Xbeta +  Wbetah


    }else if(mod == "TCTV"){

      if(!is.null(betaf)) Xbeta = array(x.fixed %*% betaf, c(nObs, m, G))  else Xbeta = 0
      Zbg = aperm(array(x.randomTC %*% t(betarTC), c(nObs, G, m)), c(1,3,2))
      Wbetah = array(x.randomTV %*% t(betarTV), c(nObs, m, G))
      linear.predictor = Xbeta + Zbg + Wbetah
      dimnames(linear.predictor)[2:3] = list(names(object$pg), names(object$delta))
    }

  }
  rownames(linear.predictor) = rown
  return(linear.predictor)
}
