% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lqr.R
\name{lqr}
\alias{lqr}
\title{Linear Quantile Regression}
\usage{
lqr(formula, data, qtl = 0.5, se = TRUE, R = 50, verbose = TRUE, ...)
}
\arguments{
\item{formula}{an object of class \code{\link{formula}}: a symbolic description of the model to be fitted}

\item{data}{a data frame containing the variables named in \code{formula} and \code{time}}

\item{qtl}{quantile to be estimated}

\item{se}{standard error computation}

\item{R}{number of bootstrap samples for computing standard errors}

\item{verbose}{if set to FALSE, no printed output is given during the function execution}

\item{...}{further arguments to be passed to or from methods}
}
\value{
Return an object of \code{\link{class}} \code{lqr}. This is a list containing the following elements:
\item{betaf}{a vector containing fixed regression coefficients}
\item{scale}{the scale parameter}
\item{sigma.e}{the standard deviation of error terms}
\item{lk}{the log-likelihood}
\item{npar}{the total number of model parameters}
\item{AIC}{the AIC value}
\item{BIC}{the BIC value}
\item{qtl}{the estimated quantile}
\item{nobs}{the total number of observations}
\item{se.betaf}{the standard errors for fixed regression coefficients}
\item{se.scale}{the standard error for the scale parameter}
\item{model}{the estimated model}
\item{call}{the matched call}
}
\description{
Estimate a linear quantile regression model with no random coefficients
}
\details{
The function computes ML estimates for the parameters of a linear quantile regression model for independent observations.
Estimates are derived by maximizing the (log-)likelihood of a Laplace regression, where the location parameter is modeled as a function
of fixed coefficients only.

If \code{se=TRUE}, standard errors based on a bootstrap procedure are computed.
}
\examples{
out0 = lqr(formula=meas~trt+time+trt:time,data=pain,se=TRUE,R=10)
}
\references{
{
  \insertRef{ref:lqr}{lqmix}
}
}
