#' Linear Quantile Mixture with TC and/or TV, discrete, random coefficients
#'
#' Estimate a finite mixture of linear quantile regression models with TC and/or TV discrete random coefficients, for a given number of components and/or states.
#'
#' @param formula an object of class \code{formula}: a symbolic description of the model to be fitted
#' @param randomTC a one-sided formula of the form \code{~z1+z2+...+zd}, where \code{z1,..., zd} denote the variables associated to TC random coefficients (1 for the intercept)
#' @param randomTV a one-sided formula of the form \code{~w1+w2+...+wl}, where \code{w1,..., wl} denote the variables associated to TV random coefficients (1 for the intercept). Note that only TC variables are allowed
#' @param group a string indicating the grouping variable, i.e., the factor identifying the unit longitudinal measurements refer to
#' @param time a string indicating the time variable
#' @param G number of mixture components associated to TC random coefficients
#' @param m number of states associated to TV random coefficients
#' @param data a data frame containing the variables named in \code{formula}, \code{randomTC}, \code{randomTV}, \code{group}, and \code{time}
#' @param qtl quantile to be estimated
#' @param eps tolerance level for convergence of the EM algorithm
#' @param maxit maximum number of iterations for the EM algorithm
#' @param se if set to TRUE, standard errors are computed
#' @param R number of bootstrap samples for computing standard errors
#' @param start type of starting values (0 = deterministic, 1 = random, 2 = initial values in input)
#' @param parInit list of initial model parameters when \code{start=2}
#' @param seed an integer value for random numbers generation, used for random parameter initialization and bootstrap standard errors
#' @param verbose if set to FALSE, no printed output is given during the function execution
#' @param parallel if set to TRUE, a parallelized code is use for standard error computation (if \code{se = TRUE})
#' @param ncores number of cores used for computing bootstrap standard errors (if required)
#'
#' @details
#'
#'
#' The function computes ML estimates for a linear quantile mixture model with TC and/or TV random coefficients.
#' Estimates are derived by maximizing the (log-)likelihood of an Asymmetric Laplace regression, where the location parameter is modeled as a function of fixed coefficients, together with TC, TV, or TC and TV discrete random coefficients, as proposed by Alfo' et. al (2017), Farcomeni (2012), and Marino et. al (2018), respectively.
#'
#' The function requires data in long-format and two additional columns indicating the group identifier and the time occasion.
#' The model is specified by means of the arguments \code{formula}, \code{formulaTC}, and \code{formulaTV}:
#' \code{formula} is associated to fixed coefficients; \code{formulaTC} is associated to TC random coefficients; \code{formulaTV} is associated to TV random coefficients.
#' In this latter, only TC variables (predictors) are allowed.
#'
#' The function allows for missing data, including dropouts (monotone missing data) and intermittent missingness, under a missing-at-random assumption.
#' Note that, when TV random coefficients are considered, intermittent missingness may cause biased inference.
#'
#' If \code{se=TRUE}, standard errors based on a block bootstrap procedure are computed.
#'
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
#' @return Return an object of \code{class} \code{lqmix}. This is a list containing the following elements:
#' \item{betaf}{a vector of fixed regression coefficients}
#' \item{betarTC}{a matrix of TC random coefficients, if present in the model}
#' \item{betarTV}{a matrix of TV random coefficients, if present in the model}
#' \item{pg}{the prior probabilities of the finite mixture associated to TC random coefficients, if present in the model}
#' \item{delta}{the initial probability vector of the hidden Markov chain associated to TV random coefficients, if present in the model}
#' \item{Gamma}{the transition probability matrix of the hidden Markov chain associated to TV random coefficients, if present in the model}
#' \item{scale}{the scale parameter}
#' \item{sigma.e}{the standard deviation of error terms}
#' \item{lk}{the log-likelihood at convergence of the EM algorithm}
#' \item{npar}{the total number of model parameters}
#' \item{aic}{the AIC value}
#' \item{bic}{the BIC value}
#' \item{qtl}{the estimated quantile}
#' \item{G}{the number of mixture components associated to TC random coefficients (if present)}
#' \item{m}{the number of hidden states associated to TV random coefficients (if present)}
#' \item{nsbjs}{the number of subjects (units)}
#' \item{nobs}{the total number of observations}
#' \item{se.betaf}{the standard errors for fixed regression coefficients}
#' \item{se.betarTC}{the standard errors for TC random coefficients, if present in the model}
#' \item{se.betarTV}{the standard errors for TV random coefficients, if present in the model}
#' \item{se.Mprob}{the standard errors for the prior probabilities of the finite mixture associated to TC random coefficients, if present in the model}
#' \item{se.Init}{the standard errors for the initial probabilities of the hidden Markov chain associated to TV random coefficients, if present in the model}
#' \item{se.Trans}{the standard errors for the transition probabilities of the hidden Markov chain associated to TV random coefficients, if present in the model}
#' \item{se.scale}{the standard error for the scale parameter}
#' \item{vcov}{the bootstrap variance-covariance matrices of the regression coefficients}
#' \item{postTC}{a matrix of size \code{nsbjs x G} with the estimated posterior probabilities for the finite mixture components associated to TC random coefficients, if present in the model}
#' \item{postTV}{a matrix of size \code{nobs x m} with the estimated posterior probabilities for the hidden states associated to TV random coefficients, if present in the model}
#' \item{miss}{the missingness type}
#' \item{model}{the estimated model}
#' \item{call}{the matched call}
#' \item{mmf}{the model matrix associated with fixed regression coefficients}
#' \item{mmrTC}{the model matrix associated with TC random coefficients, if present in the model}
#' \item{mmrTV}{the model matrix associated with TV random coefficients, if present in the model}
#' \item{y}{the model response}
#' \item{formula}{the fixed model formula}
#' \item{randomTC}{the model formula for the TC random coefficients, if present in the model}
#' \item{randomTV}{the model formula for the TV random coefficients, if present in the model}
#' \item{group}{the grouping variable}
#' \item{time}{the time variable}
#'
#' @references{
#'   \insertRef{ref:lqmixTC}{lqmix}
#' }
#' @references{
#'   \insertRef{ref:lqmixTV}{lqmix}
#' }
#'
#' @references{
#'   \insertRef{ref:lqmixTCTV}{lqmix}
#' }
#'
#' @examples
#' outTC = lqmix(formula=meas~trt+time+trt:time,randomTC=~1,
#'               group="id",time="time",G=2,data=pain,se=TRUE,R=10)
#'
#'\donttest{
#' outTV = lqmix(formula=meas~trt+time+trt:time,randomTV=~1,
#'               group="id",time="time",m=2,data=pain,R=10)
#'
#' outTCTV = lqmix(formula=meas~trt+time+trt:time,randomTC=~time,
#'                     randomTV=~1,group="id",time="time",m=2,G=2,data=pain,R=10)
#'}
#' @export

lqmix = function(formula, randomTC=NULL, randomTV=NULL, group, time, G=NULL, m=NULL, data, qtl=0.5, eps=10^-5, maxit=1000, se=TRUE, R=200, start=0, parInit=list(betaf=NULL, betarTC=NULL, betarTV=NULL, pg=NULL, delta=NULL, Gamma=NULL, scale=NULL), verbose=TRUE, seed=NULL, parallel=FALSE, ncores=2){

  # handling errors

  if(is.null(data)) stop("No input dataset has been given.")
  if(!is.data.frame(data))  stop("`data' must be a data frame.")

  if(start == 2 & is.null(unlist(parInit))) stop("No input parameters have been given with start = 2.")
  if (qtl <= 0 | qtl >= 1) stop("Quantile level out of range")

 # if(is.null(G) & is.null(m)) stop("No values for both G and m are provided.")

  if(!is.null(randomTC) & is.null(randomTV)) { # TC model
    if (is.null(G)) stop("Argument 'G' must be specified to fit a linear quantile mixture with TC random coefficients.")
    if (!is.null(m)) stop("Argument 'm' is provided, but no 'randomTV' formula is specified. To include TV random coefficients, define 'randomTV' as well.")

    model = "TC"
    oo = lqmixTC(formula=formula,randomTC=randomTC,group=group,time=time,G=G,eps=eps,data=data,qtl=qtl,verbose=verbose,maxit=maxit,se=se,R=R,start=start,parInit=parInit,seed=seed,parallel=parallel,ncores=ncores)

  }else if (is.null(randomTC) && !is.null(randomTV)) { # TV model
      if(is.null(m)) stop("Argument 'm' must be specified to fit a linear quantile mixture with TV random coefficients.")
      if(!is.null(G)) stop("Argument 'G' is provided, but no 'randomTC' formula is specified. To include TC random ciefficients, define 'randomTC' as well.")

    model = "TV"
    oo = lqmixTV(formula=formula,randomTV=randomTV,group=group,time=time,m=m,eps=eps,data=data,qtl=qtl,verbose=verbose,maxit=maxit,se=se,R=R,start=start,parInit=parInit,seed=seed,parallel=parallel,ncores=ncores)

  }else if (!is.null(randomTC) && !is.null(randomTV)) { # TCTV model
    if (is.null(G) || is.null(m)) stop("Both 'G' and 'm' must be specified to fit a linear quantile mixture with TC and TV random coefficients.")

    model = "TCTV"
    oo = lqmixTCTV(formula=formula,randomTC=randomTC,randomTV=randomTV,group=group,time=time,m=m,G=G,eps=eps,data=data,qtl=qtl,verbose=verbose,maxit=maxit,se=se,R=R,start=start,parInit=parInit,seed=seed,parallel=parallel,ncores=ncores)

  }else{
    stop("No random-coefficients formulas ('randomTC' or 'randomTV') have been specified. The model corresponds to a standard linear quantile regression without random coefficients. Please use the function lqr().")
  }

  if(!inherits(oo, "error")){
    oo$call <- match.call()
    oo$formula <- formula
    oo$randomTC <- randomTC
    oo$randomTV <- randomTV
    oo$group <- group
    oo$time <- time
    class(oo) <- "lqmix"
  }
  return(oo)

}
