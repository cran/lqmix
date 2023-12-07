#' Linear Quantile Mixture with TC and/or TV, discrete, random coefficients
#'
#' Estimate a finite mixture of linear quantile regression models with TC and/or TV, discrete, random coefficients, for a given number of components and/or states
#'
#' @param formula an object of class \code{\link{formula}}: a symbolic description of the model to be fitted
#' @param randomTC a one-sided formula of the form \code{~z1+z2+...+zr}, where \code{z1,..., zr} denote the variables associated to TC random coefficients (1 for the intercept)
#' @param randomTV a one-sided formula of the form \code{~w1+w2+...+wl}, where \code{w1,..., wl} denote the variables associated to TV random coefficients (1 for the intercept). Note that only TC variables are allowed
#' @param group a string indicating the grouping variable, i.e., the factor identifying the unit longitudinal measurements refer to
#' @param time a string indicating the time variable
#' @param G number of mixture components associated to TC random coefficients
#' @param m number of states associated to the TV random coefficients
#' @param data a data frame containing the variables named in \code{formula}, \code{randomTC}, \code{randomTV}, and \code{time}
#' @param qtl quantile to be estimated
#' @param eps tolerance level for (relative) convergence of the EM algorithm
#' @param maxit maximum number of iterations for the EM algorithm
#' @param se standard error computation for the optimal model
#' @param R number of bootstrap samples for computing standard errors
#' @param start type of starting values (0 = deterministic, 1 = random, 2 = initial values in input)
#' @param parInit list of initial model parameters when \code{start=2}. For a list of
#' @param seed an integer value for random numbers generation
#' @param verbose if set to FALSE, no printed output is given during the function execution
#' @param parallel if set to TRUE, a parallelized code is use for standard error computation (if se=TRUE)
#'
#' @details
#'
#'
#' The function computes ML estimates for the parameters of a linear quantile mixture model, based on TC and/or TV random coefficients.
#' Estimates are derived by maximizing the (log-)likelihood of a Laplace regression where the location parameter is modeled as a function of fixed coefficients, together with TC and/or TV discrete random coefficients, as proposed by Alfo' et. al (2017), Farcomeni (2012), and Marino et. al (2018), respectively.
#'
#' The function requires data in long-format and two additional columns indicating the group identifier and the time occasion.
#' The model is specified by means of the arguments \code{formula}, \code{formulaTC}, and \code{formulaTV}:
#' \code{formula} is associated to fixed coefficients; \code{formulaTC} is associated to TC random coefficients; \code{formulaTV} is associated to TV random coefficients.
#' In this latter, only TC variables (predictors) are allowed.
#'
#' The function admits the presence of missing data, both in terms of drop-outs (monotone missing data) and intermittent missing, under a missing-at-random assumption.
#' Note that, when TV random coefficients are considered, intermittent missingness may cause biased inference.
#'
#' If \code{se=TRUE}, standard errors based on a block bootstrap procedure are computed.
#'
#'
#' @import
#' quantreg
#' stats
#' methods
#' doParallel
#' foreach
#'
#'
#' @importFrom Rdpack reprompt
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#'
#' @return Return an object of \code{\link{class}} \code{lqmix}. This is a list containing the following elements:
#' \item{betaf}{a vector containing fixed regression coefficients}
#' \item{betarTC}{a matrix containing the TC random coefficients, if present in the model}
#' \item{betarTV}{a matrix containing the TV random coefficients, if present in the model}
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
#' \item{nsbjs}{the number of subjects}
#' \item{nobs}{the total number of observations}
#' \item{se.betaf}{the standard errors for fixed regression coefficients}
#' \item{se.betarTC}{the standard errors for TC random coefficients (if present)}
#' \item{se.betarTV}{the standard errors for TV random coefficients (if present)}
#' \item{se.Mprob}{the standard errors for the prior probabilities of the finite mixture associated to TC random coefficients (if present)}
#' \item{se.Init}{the standard errors for the initial probabilities of the hidden Markov chain associated to TV random coefficients(if present)}
#' \item{se.Trans}{the standard errors for the transition probabilities of the hidden Markov chain associated to TV random coefficients (if present)}
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

lqmix = function(formula,randomTC=NULL,randomTV=NULL,group,time,G=NULL,m=NULL,data,qtl=0.5,eps=10^-5,maxit=1000,se=TRUE,R=50,start=0,parInit=list(betaf=NULL,  betarTC=NULL, betarTV=NULL, pg=NULL, delta = NULL, Gamma = NULL, scale=NULL),verbose=TRUE,seed=NULL,parallel=FALSE){

  # general possible errors
  if(se==TRUE & is.null(R)){
    R = 50
    mess = "\n The number of bootstrap samples for computing standard errors has not been specified. The default value R=50 has been used."
  }else mess=NULL

  if(start == 2 & is.null(unlist(parInit))) stop("No input parameters have been given with start = 2.")

  if (qtl <= 0 | qtl >= 1) stop("Quantile level out of range")

  if(is.null(G) & is.null(m)) stop("No values for both G and m are provided. The specified model corresponds to a linear quantile regression model with no random coefficients. Please, use the function lqr().")
  if(!is.data.frame(data))  stop("`data' must be a data frame.")


  if(is.null(randomTV) & !is.null(randomTC)){ # model = "TC"

    # possible errors
    if(is.null(G) & !is.null(randomTC)) stop("The argument G has not been specified.")
    if(!is.null(G) & is.null(randomTC)) stop("The argument randomTC has not been specified.")
    if(!is.null(m)) stop("A value for the argument m has been provided. Specify the argument randomTV as well.")

    model = "TC"

    oo = lqmixTC(formula=formula,randomTC=randomTC,group=group,time=time,G=G,eps=eps,data=data,qtl=qtl,verbose=verbose,maxit=maxit,se=se,R=R,start=start,seed=seed,parallel=parallel)

  }

  else if(is.null(randomTC) & !is.null(randomTV)){ # model = "TV"

    # possible errors
    if(is.null(m) & !is.null(randomTV)) stop("The argument m has not been specified.")
    if(!is.null(m) & is.null(randomTV)) stop("The argument randomTV has not been specified.")
    if(!is.null(G)) stop("A value for the argument Gv has been provided. Specify the argument randomTC as well.")


    model = "TV"
    oo = lqmixTV(formula=formula,randomTV=randomTV,group=group,time=time,m=m,eps=eps,data=data,qtl=qtl,verbose=verbose,se=se,R=R,start=start,seed=seed,parallel=parallel)


  }else{ # model = TCTV

    if(any(is.null(G), is.null(m)) & !any(is.null(randomTC), is.null(randomTV))) stop("The argument(s) G and/or m has/have not been specified.")
    if(!any(is.null(G), is.null(m)) & any(is.null(randomTC), is.null(randomTV))) stop("The argument(s) randomTC and/or randomTV has/have not been specified.")

    model = "TCTV"
    oo = lqmixTCTV(formula=formula,randomTC=randomTC,randomTV=randomTV,group=group,time=time,m=m,G=G,eps=eps,data=data,qtl=qtl,verbose=verbose,se=se,R=R,start=start,seed=seed,parallel=parallel)
  }

  if(!inherits(oo, "error")){
    oo$call <- match.call()
    class(oo) <- "lqmix"
  }
  return(oo)


}
