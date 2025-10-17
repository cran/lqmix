\name{lqmix-package}
\alias{lqmix-package}
\docType{package}
\title{Overview of the package \code{lqmix}}

\description{
The \code{lqmix} package allows for the estimation of finite mixtures of linear quantile regression models based on Time-Constant (TC) and/or Time-Varying (TV), discrete, random coefficients for the analysis of longitudinal data.
}
\details{
\code{lqmix} is an \code{R} package devoted to the estimation of a class of linear quantile regression models for longitudinal data, in the presence of Time-Constant (TC) and/or Time-Varying (TV), unit-specific, random coefficients, having unspecific distribution.
The parameters of this distribution, together with all the others characterizing the model, are estimated in a maximum likelihood framework, via an extended Expectation-Maximization (EM) algorithm. This approach leads to the estimation of discrete distributions for the random coefficients, which give rise to a likelihood function similar to that of standard finite mixture models (in the case of TC random coefficients only), hidden Markov models (in the case of TV random coefficients only), or mixed hidden Markov models with discrete effects (in the case of both TC and TV random coefficients).

Parameters' standard errors are estimated via a block-bootstrap procedure, while model selection is performed by either maximizing the log-likelihood function, or minimizing the Akaike Information Criterion, or the Bayesian Information Criterion.

Missing data are allowed and treated under a Missing at Random assumption.
}

\author{
Maria Francesca Marino [aut,cre], Marco Alfo' [aut], Nicola Salvati [aut], and Maria Giovanna Ranalli [aut]

Maintainer: Maria Francesca Marino <mariafrancesca.marino@unifi.it>
}

\references{
  \insertRef{ref:lqmixTC}{lqmix}

  \insertRef{ref:npml1}{lqmix}

  \insertRef{ref:npml2}{lqmix}

  \insertRef{ref:lqmixTV}{lqmix}

  \insertRef{ref:LM}{lqmix}

  \insertRef{ref:hmm}{lqmix}

  \insertRef{ref:lqmixTCTV}{lqmix}

  \insertRef{ref:mhmm1}{lqmix}

  \insertRef{ref:mhmm2}{lqmix}

}
