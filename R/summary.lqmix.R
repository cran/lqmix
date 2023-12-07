#' Summary of an \code{lqmix} Object
#'
#' Summary method for the \code{\link{class}} \code{lqmix}
#'
#'
#' @param object an \code{lqmix} object
#' @param ... not used
#'
#' @return Return an object of \code{\link{class}} \code{summary.lqmix}.
#' This is a list of summary statistics for the fitted linear quantile mixture model given in \code{object}, with the following elements:
#' \item{fix}{a matrix with estimates, standard errors, Z statistics, and p-values for the fixed regression coefficients}
#' \item{ranTC}{a matrix with estimates, standard errors, Z statistics, and p-values for the TC random coefficients (if present)}
#' \item{ranTV}{a matrix with estimates, standard errors, Z statistics, and p-values for the TV random coefficients (if present)}
#' \item{pg}{a matrix with estimates and standard errors for the prior probabilities of the finite mixture associated to TC random coefficients (if present)}
#' \item{delta}{a matrix with estimates and standard errors for the initial probabilities of the hidden Markov chain associated to TV random coefficients (if present)}
#' \item{Gamma}{a matrix with estimates and standard errors for the transition probabilities of the hidden Markov chain associated to TV random coefficients (if present)}
#' \item{scale}{the scale parameter}
#' \item{sigma.e}{the standard deviation of error terms}
#' \item{logLik}{the log-likelihood at convergence of the EM algorithm}
#' \item{npar}{the total number of model parameters}
#' \item{AIC}{the AIC value}
#' \item{BIC}{the BIC value}
#' \item{qtl}{the estimated quantile}
#' \item{G}{the number of mixture components associated to TC random coefficients (if present)}
#' \item{m}{the number of hidden states associated to TV random coefficients (if present)}
#' \item{nsbj}{the number of subjects}
#' \item{nobs}{the total number of observations}
#' \item{miss}{the missingness type}
#' \item{model}{the estimated model}
#' \item{call}{the matched call}
#' @export


summary.lqmix = function(object, ...){

  if(any(!is.null(c(object$se.betaf, object$se.betarTC, object$se.betarTV)))){
    names = c("Estimate", "St.Error", "z.value", "P(>|z|)")

    if(!is.null(object$betaf)){
      est = c(object$betaf)
      sef = c(object$se.betaf)
      zvalf = c(object$betaf/sef)
      pvalf = c(1.96*pnorm(-abs(zvalf)))

      tabf = cbind(Estimate = est,
                   St.Err = sef,
                   t.value = zvalf,
                   p.value = pvalf)
      colnames(tabf) = names

    }else tabf = NULL

    if(!is.null(object$betarTC) & !is.null(object$betarTV)){ #TCTV

      # Time-Constant random coefficients
      est = c(object$betarTC)
      seTC = c(object$se.betarTC)
      zvalTC = c(object$betarTC/seTC)
      pvalTC = c(1.96*pnorm(-abs(zvalTC)))

      tabTC = cbind(Estimate = est,
                   St.Err = seTC,
                   z.value = zvalTC,
                   p.value = pvalTC)
      colnames(tabTC) = names

      rownames(tabTC) = c(sapply(colnames(object$betarTC), function(xx) paste(xx, paste("Comp", 1:nrow(object$betarTC), sep=""), sep="_")))

      # Time-Varying random coefficients
      est = c(object$betarTV)
      seTV = c(object$se.betarTV)
      zvalTV = c(object$betarTV/seTV)
      pvalTV = c(1.96*pnorm(-abs(zvalTV)))

      tabTV = cbind(Estimate = est,
                    St.Err = seTV,
                    z.value = zvalTV,
                    p.value = pvalTV)
      colnames(tabTV) = names
      rownames(tabTV) = c(sapply(colnames(object$betarTV), function(xx) paste(xx, paste("St", 1:nrow(object$betarTV), sep=""), sep="_")))

      tabMprob = cbind(Estimate = object$pg, St.Err = object$se.pg)

      tabInit = cbind(Estimate = object$delta, St.Err = object$se.delta)

      tmp1 = c(t(object$Gamma))
      names(tmp1) = paste(rep(rownames(object$Gamma), each=object$m), rep(colnames(object$Gamma), object$m), sep="")
      tmp2 = c(t(object$se.Gamma))
      names(tmp1) = paste(rep(rownames(object$Gamma), each=object$m), rep(colnames(object$Gamma), object$m), sep="")
      tabTrans = cbind(Estimate = tmp1, St.Err = tmp2)

    }else if (!is.null(object$betarTC) & is.null(object$betarTV)){ # TC


      # Time-Constant random coefficients
      est = c(object$betarTC)
      seTC = c(object$se.betarTC)
      zvalTC = c(object$betarTC/seTC)
      pvalTC = c(1.96*pnorm(-abs(zvalTC)))

      tabTC = cbind(Estimate = est,
                    St.Err = seTC,
                    z.value = zvalTC,
                    p.value = pvalTC)

      colnames(tabTC) = names
      rownames(tabTC) = c(sapply(colnames(object$betarTC), function(xx) paste(xx, paste("Comp", 1:nrow(object$betarTC), sep=""), sep="_")))

      tabTV = NULL

      tabMprob = cbind(Estimate = object$pg, St.Err = object$se.pg)

      tabInit = NULL
      tabTrans = NULL

    }else{ # TV

      # Time-Varying random coefficients
      estTV = c(object$betarTV)
      seTV = c(object$se.betarTV)
      zvalTV = c(object$betarTV/seTV)
      pvalTV = c(1.96*pnorm(-abs(zvalTV)))

      tabTV = cbind(Estimate = estTV,
                    St.Err = seTV,
                    z.value = zvalTV,
                    p.value = pvalTV)

      colnames(tabTV) = names
      rownames(tabTV) = c(sapply(colnames(object$betarTV), function(xx) paste(xx, paste("St", 1:nrow(object$betarTV), sep=""), sep="_")))

      tabTC = NULL
      tabMprob = NULL
      tabInit = cbind(Estimate = object$delta, St.Err = object$se.delta)

      tmp1 = c(t(object$Gamma))
      names(tmp1) = paste(rep(rownames(object$Gamma), each=object$m), rep(colnames(object$Gamma), object$m), sep="")
      tmp2 = c(t(object$se.Gamma))
      names(tmp1) = paste(rep(rownames(object$Gamma), each=object$m), rep(colnames(object$Gamma), object$m), sep="")
      tabTrans = cbind(Estimate = tmp1, St.Err = tmp2)

    }

    lk = object$lk
    nobs = object$nobs
    nsbjs = object$nsbjs
    mod = object$mod
    miss = object$miss
    scale = object$scale
    sigma.e = object$sigma.e
    qtl = object$qtl
    aic = object$aic
    bic = object$bic
    if(object$model == "TC") G = object$G else G = NULL
    m = object$m
    G = object$G
    npar = object$npar

    res = list()

    if(!is.null(tabf)) res$fix = tabf
    if(!is.null(tabTC)) res$ranTC = tabTC
    if(!is.null(tabTV)) res$ranTV = tabTV

    if(!is.null(tabMprob)) res$pg = tabMprob
    if(!is.null(tabInit)) res$delta = tabInit
    if(!is.null(tabTrans)) res$Gamma = tabTrans

    res$scale = scale
    res$sigma.e = sigma.e
    res$lk = lk
    res$npar = npar
    res$aic = aic
    res$bic = bic
    res$qtl = qtl
    res$G = G
    res$m = m
    res$nsbjs = nsbjs
    res$nobs = nobs


    res$miss = miss
    res$model = mod
    if(!is.null(object$call)) res$call = match.call()

    class(res) = "summary.lqmix"
    return(res)

  }else{
    print(object)
    cat("\n")
    message("Model inference not allowed: standard errors have not been computed.")
  }
}
