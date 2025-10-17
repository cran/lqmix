#' Summary of a \code{search_lqmix} object
#'
#' Summary method for the \code{\link{class}} \code{\link{search_lqmix}}.
#'
#'
#' @param object a \code{search_lqmix} object
#' @param ... not used
#'
#' @return Return an object of \code{\link{class}} \code{summary.search_lqmix}.
#' This is a list of summary statistics for the optimal linear quantile mixture model given in \code{object$optimal}, with the following elements:
#' \item{fix}{a matrix with estimates, standard errors, Z statistics, and p-values for the fixed regression coefficients for the optimal fitted model}
#' \item{ranTC}{a matrix with estimates, standard errors, Z statistics, and p-values for the TC random coefficients, if present for the optimal fitted model}
#' \item{ranTV}{a matrix with estimates, standard errors, Z statistics, and p-values for the TV random coefficients, if present for the optimal fitted model}
#' \item{pg}{a matrix with estimates and standard errors for the prior probabilities of the finite mixture associated to TC random coefficients, if present for the optimal fitted model}
#' \item{delta}{a matrix with estimates and standard errors for the initial probabilities of the hidden Markov chain associated to TV random coefficients, if present for the optimal fitted model}
#' \item{Gamma}{a matrix with estimates and standard errors for the transition probabilities of the hidden Markov chain associated to TV random coefficients, if present for the optimal fitted model}
#' \item{scale}{the scale parameter for the optimal model}
#' \item{sigma.e}{the standard deviation of error terms for the optimal model}
#' \item{logLik}{the log-likelihood at convergence of the EM algorithm for the optimal model}
#' \item{npar}{the total number of model parameters for the optimal model}
#' \item{AIC}{the AIC value for the optimal model}
#' \item{BIC}{the BIC value for the optimal model}
#' \item{qtl}{the estimated quantile}
#' \item{G}{the number of mixture components associated to TC random coefficients, if present for the optimal fitted model}
#' \item{m}{the number of hidden states associated to TV random coefficients, if present for the optimal fitted model}
#' \item{nsbj}{the number of subjects}
#' \item{nobs}{the total number of observations}
#' \item{miss}{the missingness type}
#' \item{model}{the identified optimal model}
#' \item{call}{the matched call}
#'
#' @export


summary.search_lqmix = function(object, ...){

  opt = object$optimal

  if(any(!is.null(c(opt$se.betaf, opt$se.betarTC, opt$se.betarTV)))){
    names = c("Estimate", "St.Error", "z.value", "P(>|z|)")

    if(!is.null(opt$betaf)){
      est = c(opt$betaf)
      sef = c(opt$se.betaf)
      zvalf = c(opt$betaf/sef)
      pvalf = c(1.96*pnorm(-abs(zvalf)))

      tabf = cbind(Estimate = est,
                   St.Error = sef,
                   t.value = zvalf,
                   p.value = pvalf)
      colnames(tabf) = names

    }else tabf = NULL

    if(!is.null(opt$betarTC) & !is.null(opt$betarTV)){ #TCTV

      # Time-Constant random coefficients
      est = c(opt$betarTC)
      seTC = c(opt$se.betarTC)
      zvalTC = c(opt$betarTC/seTC)
      pvalTC = c(1.96*pnorm(-abs(zvalTC)))

      tabTC = cbind(Estimate = est,
                   St.Error = seTC,
                   z.value = zvalTC,
                   p.value = pvalTC)
      colnames(tabTC) = names

      rownames(tabTC) = c(sapply(colnames(opt$betarTC), function(xx) paste(xx, paste("Comp", 1:nrow(opt$betarTC), sep=""), sep="_")))

      # Time-Varying random coefficients
      est = c(opt$betarTV)
      seTV = c(opt$se.betarTV)
      zvalTV = c(opt$betarTV/seTV)
      pvalTV = c(1.96*pnorm(-abs(zvalTV)))

      tabTV = cbind(Estimate = est,
                    St.Error = seTV,
                    z.value = zvalTV,
                    p.value = pvalTV)
      colnames(tabTV) = names
      rownames(tabTV) = c(sapply(colnames(opt$betarTV), function(xx) paste(xx, paste("St", 1:nrow(opt$betarTV), sep=""), sep="_")))

      tabMprob = cbind(Estimate = opt$pg, St.Error = opt$se.pg)

      tabInit = cbind(Estimate = opt$delta, St.Error = opt$se.delta)

      tmp1 = c(t(opt$Gamma))
      names(tmp1) = paste(rep(rownames(opt$Gamma), each=opt$m), rep(colnames(opt$Gamma), opt$m), sep="")
      tmp2 = c(t(opt$se.Gamma))
      names(tmp1) = paste(rep(rownames(opt$Gamma), each=opt$m), rep(colnames(opt$Gamma), opt$m), sep="")
      tabTrans = cbind(Estimate = tmp1, St.Error = tmp2)

    }else if (!is.null(opt$betarTC) & is.null(opt$betarTV)){ # TC


      # Time-Constant random coefficients
      est = c(opt$betarTC)
      seTC = c(opt$se.betarTC)
      zvalTC = c(opt$betarTC/seTC)
      pvalTC = c(1.96*pnorm(-abs(zvalTC)))

      tabTC = cbind(Estimate = est,
                    St.Error = seTC,
                    z.value = zvalTC,
                    p.value = pvalTC)

      colnames(tabTC) = names
      rownames(tabTC) = c(sapply(colnames(opt$betarTC), function(xx) paste(xx, paste("Comp", 1:nrow(opt$betarTC), sep=""), sep="_")))

      tabTV = NULL

      tabMprob = cbind(Estimate = opt$pg, St.Error = opt$se.pg)

      tabInit = NULL
      tabTrans = NULL

    }else{ # TV

      # Time-Varying random coefficients
      estTV = c(opt$betarTV)
      seTV = c(opt$se.betarTV)
      zvalTV = c(opt$betarTV/seTV)
      pvalTV = c(1.96*pnorm(-abs(zvalTV)))

      tabTV = cbind(Estimate = estTV,
                    St.Error = seTV,
                    z.value = zvalTV,
                    p.value = pvalTV)

      colnames(tabTV) = names
      rownames(tabTV) = c(sapply(colnames(opt$betarTV), function(xx) paste(xx, paste("St", 1:nrow(opt$betarTV), sep=""), sep="_")))

      tabTC = NULL
      tabMprob = NULL
      tabInit = cbind(Estimate = opt$delta, St.Error = opt$se.delta)

      tmp1 = c(t(opt$Gamma))
      names(tmp1) = paste(rep(rownames(opt$Gamma), each=opt$m), rep(colnames(opt$Gamma), opt$m), sep="")
      tmp2 = c(t(opt$se.Gamma))
      names(tmp1) = paste(rep(rownames(opt$Gamma), each=opt$m), rep(colnames(opt$Gamma), opt$m), sep="")
      tabTrans = cbind(Estimate = tmp1, St.Error = tmp2)

    }

    lk = opt$lk
    nobs = opt$nobs
    nsbjs = opt$nsbjs
    mod = opt$mod
    miss = opt$miss
    scale = opt$scale
    sigma.e = opt$sigma.e
    qtl = opt$qtl
    aic = opt$aic
    bic = opt$bic
    if(opt$model == "TC") G = opt$G else G = NULL
    m = opt$m
    G = opt$G
    npar = opt$npar

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
    if(!is.null(opt$call)) res$call = match.call()

    class(res) = "summary.search_lqmix"
    return(res)

  }else{
    print(opt)
     cat("\n")
    message("Model inference not allowed: standard errors have not been computed.")
  }
}
