#' Density of the Asymmetric Laplace distribution
#'
#' Compute the density for the three parameter Asymmetric Laplace Distribution
#'
#' @param y vector of quantiles
#' @param mu location parameter
#' @param sigma scale parameter
#' @param qtl skewness parameter
#' @param log logical; if TRUE, probabilities are log-transformed
#'
#' @return Return the density for the asymmetric Laplace distribution
#'
#' @details
#' The function computes the density of the Asymmetric Laplace distribution, with location \eqn{\mu}, scale \eqn{\sigma > 0}
#' and skewness \code{qtl = q} in (0,1), as discussed by Koenker and Machado (1999) and Yu and Moyeed (2001), according to the following expression
#'
#' \deqn{f(y | \mu, \sigma, q) = \frac{q(1-q)}{\sigma} \exp(-\rho_{q} (\frac{y-\mu}{\sigma}))}
#'
#' @references{
#'   \insertRef{ref:dal1}{lqmix}
#' }
#' @references{
#'   \insertRef{ref:dal2}{lqmix}
#' }
#'
#' @export


dal = function (y, mu = 0, sigma = 1, qtl = 0.5, log = FALSE){

  eps = .Machine$double.eps^(2/3)
  if (qtl > 1 | qtl < 0) stop("Parameter 'qtl' must be in [0,1]")
  if (qtl == 0) qtl = eps
  if (qtl == 1) qtl = 1 - eps
  if (sigma < 0) warning("Scale parameter 'sigma' is negative")
  ind = ifelse(y < mu, 1, 0)
  val = log(qtl) + log(1 - qtl) - log(sigma) -((y - mu)/sigma * (qtl - ind))
  if (!log) exp(val) else val

}
