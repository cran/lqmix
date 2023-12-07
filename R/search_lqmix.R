#' Search the Global Maximum of a Linear Quantile Mixture
#'
#' Search the global maximum of the log-likelihood function for a finite mixture of linear quantile regression models with TC and/or TV, discrete, random coefficients, for varying number of components and/or states
#'
#' @param formula an object of \code{\link{class}} \code{formula}: a symbolic description of the model to be fitted
#' @param randomTC a one-sided formula of the form \code{~z1+z2+...+zr}, where \code{z1,..., zr} denote the variables associated to TC random coefficients (1 for the intercept)
#' @param randomTV a one-sided formula of the form \code{~w1+w2+...+wl}, where \code{w1,..., wl} denote the variables associated to TV random coefficients (1 for the intercept). Note that only TC variables are allowed
#' @param group a string indicating the grouping variable, i.e., the factor identifying the unit longitudinal measurements refer to
#' @param time a string indicating the time variable
#' @param Gv vector of possible number of mixture components associated to TC random coefficients (if present)
#' @param mv vector of possible number of states associated to the TV random coefficients (if present)
#' @param method method to use for selecting the optimal model. Possible values are \code{"lk"}, \code{"aic"}, or \code{"bic"}
#' @param data a data frame containing the variables named in \code{formula}, \code{randomTC}, \code{randomTV}, and \code{time}
#' @param qtl quantile to be estimated
#' @param eps tolerance level for (relative) convergence of the EM algorithm
#' @param maxit maximum number of iterations for the EM algorithm
#' @param se standard error computation for the optimal model
#' @param R number of bootstrap samples for computing standard errors
#' @param nran number of repetitions of each random initialization
#' @param seed an integer value for random numbers generation
#' @param verbose if set to FALSE, no printed output is given during the function execution
#' @param parallel if set to TRUE, a parallelized code is use for standard error computation (if se=TRUE)
#'
#' @details
#' The function allows to identify the optimal model specification in terms of number of mixture components and/or hidden states associated to
#' TC and/or TV random coefficients, respectively.
#' This is done by considering a multi-start strategy based on both deterministic and random starting points.
#' The number or random tries is proportional to the number of mixture components and/or hidden states associated to the random coefficients in the model.
#'
#' If \code{method="lk"}, the optimal model selected by the function is that providing the highest log-likelihood value;
#' if \code{method="AIC"}, (\code{method="BIC"}, respectively), the optimal model selected by the function is that providing the lowest AIC (BIC, respectively) value.
#'
#' If \code{se=TRUE}, standard errors based on a block bootstrap procedure are computed for the identified optimal model.
#'
#' @import
#' quantreg
#' stats
#' methods
#' doParallel
#' foreach
#'
#'
#' @return Return an object of \code{\link{class}} \code{search_lqmix}. This is a list containing the following elements:
#' \item{optimal}{the identified optimal model}
#' \item{allmodels}{the output of each estimated model}
#' \item{lkv}{the vector of likelihood values for each estimated model}
#' \item{aicv}{the vector of AIC values for each estimated model}
#' \item{bicv}{the vector of BIC values for each estimated model}
#' \item{qtl}{the estimated quantile}
#' \item{mv}{the vector of possible number of states associated to TV random coefficients (if present)}
#' \item{Gv}{the vector of possible number of mixture components associated to TC random coefficients (if present)}
#' \item{method}{the method used to select the optimal model}
#' \item{call}{the matched call}
#'
#' @examples
#' sTC = search_lqmix(formula=meas~trt+time+trt:time,
#'                    randomTC=~1,group="id",time="time",Gv=1:3,method="bic",data=pain,se=FALSE)
#'\donttest{
#' sTV = search_lqmix(formula=meas~trt+time+trt:time,
#' randomTV=~1,group="id",time="time",mv=1:3,method="bic",data=pain,se=FALSE)
#'
#' sTCTV = search_lqmix(formula=meas~trt+time+trt:time,
#' randomTC=~time,randomTV=~1,group="id",time="time",mv=1:3,Gv=1:3,method="bic",data=pain,se=FALSE)
#'}
#' @export

search_lqmix = function(formula,randomTC=NULL,randomTV=NULL,group,time,Gv=NULL,mv=NULL,data,method="bic",nran=0,qtl=0.5,eps=10^-5,maxit=1000,se=TRUE,R=50,verbose=TRUE,seed=NULL,parallel=FALSE){

  # possible errors
  if(is.null(Gv) & is.null(mv)) stop("No values for both Gv and mv are provided.")
  if(!(method %in% c("lk","bic","aic"))) stop("The method specified for selecting the optimal model is not supported.")


  if(verbose){
    cat("Search the optimal linear quantile mixture model", "\n")
    cat("*************************************************", "\n")
  }

  # n. of bootstrap samples for se computation = 50 by default
  if(se & is.null(R)) R = 50

  if(is.null(randomTV) & !is.null(randomTC)){ # model = "TC"
    # possible errors
    if(is.null(Gv) & !is.null(randomTC)) stop("The argument Gv has not been specified.")
    if(!is.null(Gv) & is.null(randomTC)) stop("The argument randomTC has not been specified.")
    if(!is.null(mv)) stop("A value for the argument mv has been provided. Specify the argument randomTV as well.")


    model = "TC"
    Gv = unique(sort(Gv))
    out = vector("list",1)

    lkv = bicv = aicv = c()
    out[[1]] = vector("list", max(Gv))
    for(g in Gv){
      if(g == 1){
        if(verbose){ cat("Model homogeneous", "- qtl =", qtl,"\n")
          if(nran>0)cat("Random start: 0 ... \n")
        }
        out[[1]][[g]] = try(lqr(formula=formula,data=data,verbose=verbose,qtl=qtl,se=F,search=T))
      }else if(nran>0){
        if(verbose) cat("Model TC - qtl =", qtl, "\n")
        oo = try(lqmixTC(formula=formula,randomTC=randomTC,group=group,time=time,G=g,eps=10^-2,data=data,qtl=qtl,verbose=F,maxit=maxit,se=F))
        if(verbose) cat("Random start: 0 ... ")
        for(h in 1:(nran*(g-1))){
          if(verbose) cat(h, " ... ")
          if(!is.null(seed)) seed=seed+h
          ooh = try(lqmixTC(formula=formula,randomTC=randomTC,group=group,time=time,G=g,eps=10^-2,data=data,qtl=qtl,start=1,verbose=F,maxit=maxit,se=F,seed=seed))
          if(ooh$lk > oo$lk) oo = ooh
        }
        if(verbose)cat("\n")
        out[[1]][[g]] = lqmixTC(formula=formula,randomTC=randomTC,group=group,time=time,G=g,data=data,qtl=qtl,start=2,parInit=oo,verbose=verbose,maxit=maxit,eps=eps,se=F,search=T)
      }else{
        if(verbose) cat("Model TC - qtl =", qtl, "\n")
        out[[1]][[g]] = lqmixTC(formula=formula,randomTC=randomTC,group=group,time=time,G=g,data=data,qtl=qtl,start=0,verbose=verbose,maxit=maxit,eps=eps,se=F,search=T)
      }

      if(!inherits(out[[1]][[g]], "try-error")){
        lkv = c(lkv, out[[1]][[g]]$lk)
        bicv = c(bicv, out[[1]][[g]]$BIC)
        aicv = c(aicv, out[[1]][[g]]$AIC)
      }else{
        lkv = c(lkv, NA)
        bicv = c(bicv, NA)
        aicv = c(aicv, NA)
      }

    }
  }
  else if(is.null(randomTC) & !is.null(randomTV)){ # model = "TV"
    # possible errors
    if(is.null(mv) & !is.null(randomTV)) stop("The argument mv has not been specified.")
    if(!is.null(mv) & is.null(randomTV)) stop("The argument randomTV has not been specified.")
    if(!is.null(Gv)) stop("A value for the argument Gv has been provided. Specify the argument randomTC as well.")


    model = "TV"
    mv = unique(sort(mv))
    out = vector("list",max(mv))
    lkv = bicv = aicv = c()

    for(m in mv){
      out[[m]] = vector("list",1)
      if(m == 1){
        if(verbose){
          cat("Model homogeneous", "- qtl =", qtl,"\n")
          if(nran>0)cat("Random start: 0 ... \n")
        }
        out[[m]][[1]] = try(lqr(formula=formula,data=data,verbose=verbose,qtl=qtl,se=F,search=T))
      }else if(nran>0){
        if(verbose) cat("Model TV - qtl =", qtl, "\n")
        oo = try(lqmixTV(formula=formula,randomTV=randomTV,group=group,time=time,m=m,eps=10^-2,data=data,qtl=qtl,verbose=FALSE,se=F))
        if(verbose) cat("Random start: 0 ... ")
        for(h in 1:(nran*(m-1))){
          if(verbose) cat(h, " ... ")
          if(!is.null(seed)) seed=seed+h
          ooh = try(lqmixTV(formula=formula,randomTV=randomTV,group=group,time=time,m=m,eps=10^-2,data=data,qtl=qtl,start=1,verbose=FALSE,se=F,seed=seed))
          if(ooh$lk > oo$lk) oo = ooh
        }
        if(verbose)cat("\n")
        out[[m]][[1]] = try(lqmixTV(formula=formula,randomTV=randomTV,group=group,time=time,m=m,data=data,qtl=qtl,start=2,parInit=oo,verbose=verbose,maxit=maxit,eps=eps,se=F,search=T))
      }else{
        if(verbose) cat("Model TV - qtl =", qtl, "\n")
        out[[m]][[1]] = try(lqmixTV(formula=formula,randomTV=randomTV,group=group,time=time,m=m,data=data,qtl=qtl,start=0,verbose=verbose,maxit=maxit,eps=eps,se=F,search=T))
      }

      if(!inherits(out[[m]][[1]], "try-error")){

        lkv = c(lkv, out[[m]][[1]]$lk)
        bicv = c(bicv, out[[m]][[1]]$BIC)
        aicv = c(aicv, out[[m]][[1]]$AIC)
      }else{
        lkv = c(lkv, NA)
        bicv = c(bicv, NA)
        aicv = c(aicv, NA)
      }
    }

  }else{

    if(any(is.null(Gv), is.null(mv)) & !any(is.null(randomTC), is.null(randomTV))) stop("The argument(s) Gv and/or mv has/have not been specified.")
    if(!any(is.null(Gv), is.null(mv)) & any(is.null(randomTC), is.null(randomTV))) stop("The argument(s) randomTC and/or randomTV has/have not been specified.")

    model = "TCTV"
    mv = sort(unique(mv))
    Gv = sort(unique(Gv))

    out = vector("list",max(mv))
    lkv = bicv = aicv = c()

    for(m in mv){
      out[[m]] =  vector("list",max(Gv))
      for(g in sort(Gv)){
        if(m == 1 & g == 1){
          if(verbose){
            cat("Model homogeneous", "- qtl =", qtl,"\n")
            if(nran>0) cat("Random start: 0 ... \n")
          }
          out[[m]][[g]] = try(lqr(formula=formula,data=data,verbose=verbose,qtl=qtl,se=F,search=T)) # model Homogeneous
        }else if(m > 1 & g == 1){ # model = TV
          if(verbose) cat("Model TV", "- qtl =", qtl,"\n")
          if(nran>0){
            oo = try(lqmixTV(formula=formula,randomTV=randomTV,group=group,time=time,m=m,eps=10^-2,data=data,qtl=qtl,verbose=FALSE,se=F))
            if(verbose) cat("Random start: 0 ... ")
            for(h in 1:(nran*(m-1))){
              if(verbose) cat(h, " ... ")
              if(!is.null(seed)) seed=seed+h
              ooh = try(lqmixTV(formula=formula,randomTV=randomTV,group=group,time=time,m=m,eps=10^-2,data=data,qtl=qtl,start=1,verbose=FALSE,se=F,seed=seed))
              if(ooh$lk > oo$lk) oo = ooh
            }
            if(verbose)cat("\n")
            out[[m]][[g]] = try(lqmixTV(formula=formula,randomTV=randomTV,group=group,time=time,m=m,data=data,qtl=qtl,start=2,parInit=oo,verbose=verbose,maxit=maxit,eps=eps,se=F,search=T))
          }else out[[m]][[g]] = try(lqmixTV(formula=formula,randomTV=randomTV,group=group,time=time,m=m,data=data,qtl=qtl,,start=0,verbose=verbose,maxit=maxit,eps=eps,se=F,search=T))

        }else if(m == 1 & g > 1){ # model = TC
          if(verbose) cat("Model TC - qtl =", qtl, "\n")
          if(nran>0){
            oo = try(lqmixTC(formula=formula,randomTC=randomTC,group=group,time=time,G=g,eps=10^-2,data=data,qtl=qtl,verbose=FALSE,se=F))
            if(verbose) cat("Random start: 0 ... ")
            for(h in 1:(nran*(g-1))){
              if(verbose) cat(h, " ... ")
              if(!is.null(seed)) seed=seed+h
              ooh = try(lqmixTC(formula=formula,randomTC=randomTC,group=group,time=time,G=g,eps=10^-2,data=data,qtl=qtl,start=1,verbose=FALSE,se=F, seed=seed))
              if(ooh$lk > oo$lk) oo = ooh
            }
            if(verbose) cat("\n")
            out[[m]][[g]] = lqmixTC(formula=formula,randomTC=randomTC,group=group,time=time,G=g,data=data,qtl=qtl,start=2,parInit=oo,verbose=verbose,maxit=maxit,eps=eps,se=F,search=T)
          }else out[[m]][[g]] = lqmixTC(formula=formula,randomTC=randomTC,group=group,time=time,G=g,data=data,qtl=qtl,start=0,verbose=verbose,maxit=maxit,eps=eps,se=F,search=T)

        }else{ # model = TCTV
          if(verbose) cat("Model TCTV - qtl =", qtl, "\n")
          if(nran>0){
            oo = try(lqmixTCTV(formula=formula,randomTC=randomTC,randomTV=randomTV,group=group,time=time,m=m,G=g,eps=10^-2,data=data,qtl=qtl,verbose=F,se=F))
            if(verbose) cat("Random start: 0 ... ")
            for(h in 1:(nran*(m-1)*(g-1))){
              if(verbose)cat(h, " ... ")
              if(!is.null(seed)) seed=seed+h
               ooh = try(lqmixTCTV(formula=formula,randomTC=randomTC,randomTV=randomTV,group=group,time=time,m=m,G=g,eps=10^-2,data=data,qtl=qtl,start=1,verbose=F,se=F,seed=seed))
              if(ooh$lk > oo$lk) oo = ooh
            }
            if(verbose)cat("\n")
            out[[m]][[g]] = try(lqmixTCTV(formula=formula,randomTC=randomTC,randomTV=randomTV,group=group,time=time,m=m,G=g,data=data,qtl=qtl,start=2,parInit=oo,verbose=verbose,maxit=maxit,eps=eps,se=F,search=TRUE))
          }else out[[m]][[g]] = try(lqmixTCTV(formula=formula,randomTC=randomTC,randomTV=randomTV,group=group,time=time,m=m,G=g,data=data,qtl=qtl,start=0,verbose=verbose,maxit=maxit,eps=eps,se=F,search=TRUE))
        }

        if(!inherits(out[[m]][[g]], "try-error")){
          lkv = c(lkv, out[[m]][[g]]$lk)
          bicv = c(bicv, out[[m]][[g]]$BIC)
          aicv = c(aicv, out[[m]][[g]]$AIC)
        }else{
          lkv = c(lkv, NA)
          bicv = c(bicv, NA)
          aicv = c(aicv, NA)
        }

      }
    }

  }

  if(model == "TCTV") names(bicv) = names(aicv) = names(lkv) = paste(paste("m=",rep(mv, each = length(Gv)), sep=""), paste("G=",rep(Gv, length(mv)), sep=""), sep="-")
  else if(model == "TC") names(bicv) = names(aicv) = names(lkv) = paste("m=1-", paste("G=",Gv, sep=""), sep="")
  else names(bicv) = names(aicv) = names(lkv) = paste(paste("m=",mv, sep=""),"G=1", sep="-")

  # selection of the optimal model
  if(method == "lk") wh = which.max(lkv)
  else if(method == "aic") wh = which.min(aicv)
  else wh = which.min(bicv)


  if(model == "TCTV"){
    whTC = as.numeric(strsplit(strsplit(names(bicv)[wh], "-")[[1]][2],"=")[[1]][2])
    whTV = as.numeric(strsplit(strsplit(names(bicv)[wh], "-")[[1]][1],"=")[[1]][2])

    if(se){
      if(verbose) cat("Computing standard errors for the optimal model ... ")

      if(parallel==TRUE) registerDoParallel(2) else registerDoParallel(1)

      if(whTC == 1 & whTV == 1){
        if(verbose) optimal.model = try(lqr(formula=formula,data=data,verbose=F,qtl=qtl,se=T,R=R,search=T,opt=TRUE)) # model Homogeneous
        else optimal.model = try(lqr(formula=formula,data=data,verbose=F,qtl=qtl,se=T,R=R,opt=TRUE)) # model Homogeneous
      }else if(whTC == 1 & whTV > 1){
        modInit = out[[whTV]][[1]]
        if(verbose) optimal.model = try(lqmixTV(formula=formula,randomTV=randomTV,group=group,time=time,m=whTV,data=data,qtl=qtl,start=2,parInit=modInit,verbose=F,maxit=maxit,eps=eps,se=T,R=R,search=T,opt=TRUE,parallel=parallel))
        else optimal.model = try(lqmixTV(formula=formula,randomTV=randomTV,group=group,time=time,m=whTV,data=data,qtl=qtl,start=2,parInit=modInit,verbose=F,maxit=maxit,eps=eps,se=T,R=R,opt=TRUE,parallel=parallel))
       }else if(whTC > 1 & whTV == 1){
         modInit = out[[1]][[whTC]]
        if(verbose) optimal.model = try(lqmixTC(formula=formula,randomTC=randomTC,group=group,time=time,G=whTC,data=data,qtl=qtl,start=2,parInit=modInit,verbose=F,maxit=maxit,eps=eps,se=T,R=R,search=T,opt=TRUE,parallel=parallel))
        else optimal.model = try(lqmixTC(formula=formula,randomTC=randomTC,group=group,time=time,G=whTC,data=data,qtl=qtl,start=2,parInit=modInit,verbose=F,maxit=maxit,eps=eps,se=T,R=R,opt=TRUE,parallel=parallel))
      }else{
         modInit = out[[whTV]][[whTC]]

         if(verbose) optimal.model = try(lqmixTCTV(formula=formula,randomTC=randomTC,randomTV=randomTV,group=group,time=time,m=whTV,G=whTC,data=data,qtl=qtl,start=2,parInit=modInit,verbose=F,eps=eps,se=T,R=R,search=T,opt=TRUE,parallel=parallel))
         else optimal.model = try(lqmixTCTV(formula=formula,randomTC=randomTC,randomTV=randomTV,group=group,time=time,m=whTV,G=whTC,data=data,qtl=qtl,start=2,parInit=modInit,verbose=F,eps=eps,se=T,R=R,maxit=1,opt=TRUE,parallel=parallel))
      }
    }else optimal.model = out[[whTV]][[whTC]]

    optimal.model$call = NULL

  }else if(model == "TC"){ # TC

    whTC = as.numeric(strsplit(strsplit(names(bicv)[wh], "-")[[1]][2],"=")[[1]][2])

    if(se){
      if(verbose) cat("Computing standard errors for the optimal model ... ")
      if(whTC == 1){
        if(verbose) optimal.model = try(lqr(formula=formula,data=data,qtl=qtl,se=T,R=R,verbose=F,search=T,opt=TRUE)) # model Homogeneous
        else optimal.model = try(lqr(formula=formula,data=data,qtl=qtl,se=T,R=R,verbose=F,opt=TRUE)) # model Homogeneous
      }else{
        modInit = out[[1]][[whTC]]
        if(verbose) optimal.model = try(lqmixTC(formula=formula,randomTC=randomTC,group=group,time=time,G=whTC,data=data,qtl=qtl,start=2,parInit=modInit,verbose=F,maxit=maxit,eps=eps,se=T,R=R,search=T,opt=TRUE,parallel=parallel))
        else optimal.model = try(lqmixTC(formula=formula,randomTC=randomTC,group=group,time=time,G=whTC,data=data,qtl=qtl,start=2,parInit=modInit,verbose=F,maxit=maxit,eps=eps,se=T,R=R,opt=TRUE,parallel=parallel))
      }
    }else optimal.model = out[[1]][[whTC]]
    optimal.model$call = NULL

  }else if(model == "TV"){ #TV
    whTV = as.numeric(strsplit(strsplit(names(bicv)[wh], "-")[[1]][1],"=")[[1]][2])

    if(se){
      if(verbose) cat("Computing standard errors for the optimal model ... ")
      if(whTV == 1){
        if(verbose) optimal.model = try(lqr(formula=formula,data=data,verbose=F,qtl=qtl,se=T,R=R,search=T,opt=TRUE)) # model Homogeneous
        else optimal.model = try(lqr(formula=formula,data=data,verbose=F,qtl=qtl,se=T,R=R,opt=TRUE)) # model Homogeneous
      }else{
        modInit = out[[whTV]][[1]]
        if(verbose) optimal.model = try(lqmixTV(formula=formula,randomTV=randomTV,group=group,time=time,m=whTV,data=data,qtl=qtl,start=2,parInit=modInit,verbose=F,maxit=-1,eps=eps,se=T,R=R,search=T,opt=TRUE,parallel=parallel))
        else optimal.model = try(lqmixTV(formula=formula,randomTV=randomTV,group=group,time=time,m=whTV,data=data,qtl=qtl,start=2,parInit=modInit,verbose=F,maxit=-1,eps=eps,se=T,R=R,opt=TRUE,parallel=parallel))
      }
    }else optimal.model = out[[whTV]][[1]]
    optimal.model$call = NULL
  }
  if(model == "TCTV") o = list(optimal=optimal.model, allmodels=out,lkv=lkv,aicv=aicv,bicv=bicv,qtl=qtl,mv=mv,Gv=Gv,method=method)
  else if(model == "TC") o = list(optimal=optimal.model, allmodels=out,lkv=lkv,aicv=aicv,bicv=bicv,qtl=qtl,Gv=Gv,method=method)
  else o = list(optimal=optimal.model, allmodels=out,lkv=lkv,aicv=aicv,bicv=bicv,qtl=qtl,mv=mv,method=method)


  o$call <- match.call()
  class(o) <- "search_lqmix"
  return(o)


}
