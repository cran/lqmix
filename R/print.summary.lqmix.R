#' Print the summary of an \code{lqmix} object
#'
#' Print the summary of an object of \code{\link{class}} \code{\link{lqmix}}.
#'
#' @param x a summary of an \code{lqmix} object
#' @param digits a non-null value for digits specifying the minimum number of significant digits to be printed
#' @param ... not used
#'
#' @return Return a summary of an \code{\link{lqmix}} object.
#'
#' @export


print.summary.lqmix = function(x, digits = max (3, getOption("digits") -3), ...){

  if(!is.null(x$call)){
    if(x$model == "TC"){
      cat(paste("Model: TC random coefficients with G=", x$G, " at qtl=", x$q, sep=""))
      cat("\n")
      cat("********************************************************", "\n")
    }else if(x$model == "TV"){
      cat(paste("Model: TV random coefficients with m=", x$m, " at qtl=", x$q, sep=""))
      cat("\n")
      cat("********************************************************", "\n")
    }else{
      cat(paste("Model: TC and TV random coefficients with m=", x$m, " and G=", x$G, " at qtl=", x$q, sep=""))
      cat("\n")
      cat("***********************************************************************", "\n")
    }
  }else{
    if(x$model == "TC"){
      cat(paste("Opt model: TC random coefficients with G=", x$G, " at qtl=", x$q, sep=""))
      cat("\n")
      cat("**********************************************************", "\n")
    }else if(x$model == "TV"){
      cat(paste("Opt model: TV random coefficients with m=", x$m, " at qtl=", x$q, sep=""))
      cat("\n")
      cat("**********************************************************", "\n")
    }else{
      cat(paste("Opt model: TC and TV random coefficients with m=", x$m, " and G=", x$G, " at qtl=", x$q, sep=""))
      cat("\n")
      cat("***************************************************************************", "\n")
    }
  }

 cat("\n---- Observed process ----\n")

  if(!is.null(x$fix)){
    cat("\nFixed Coefficients:\n")
    printCoefmat(round(x$fix, digits), P.values=TRUE, has.Pvalue=TRUE,
                 signif.stars = TRUE, signif.legend = F)
  }

  if(!is.null(x$ranTC) & !is.null(x$ranTV)){ #TCTV

    cat("\nTime-Constant Random Coefficients:\n")
    printCoefmat(round(x$ranTC, digits), P.values=TRUE, has.Pvalue=TRUE,
                 signif.stars = TRUE, signif.legend = F)

    cat("\nTime-Varying Random Coefficients:\n")
    printCoefmat(round(x$ranTV, digits), P.values=TRUE, has.Pvalue=TRUE,
                 signif.stars = TRUE, signif.legend = TRUE)

    cat("\nResidual scale parameter:", round(x$scale, digits), "- Residual standard deviation:", round(x$sigma.e, digits), "\n")

    cat("\n---- Latent process ----\n")

    cat("\nMixture probabilities:\n")
    printCoefmat(round(x$pg, digits), P.values=FALSE, has.Pvalue=FALSE,
                 signif.stars = FALSE, signif.legend = FALSE)

    cat("\nInitial probabilities:\n")
    printCoefmat(round(x$delta, digits), P.values=FALSE, has.Pvalue=FALSE,
                 signif.stars = FALSE, signif.legend = FALSE)

    cat("\nTransition probabilities:\n")
    printCoefmat(round(x$Gamma, digits), P.values=FALSE, has.Pvalue=FALSE,
                 signif.stars = FALSE, signif.legend = FALSE)

  }else if(!is.null(x$ranTC) & is.null(x$ranTV)){ #TC

    cat("\nTime-Constant Random Coefficients:\n")
    printCoefmat(round(x$ranTC, digits), P.values=TRUE, has.Pvalue=TRUE,
                 signif.stars = TRUE, signif.legend = TRUE)

    cat("\nResidual scale parameter:", round(x$scale, digits), "- Residual standard deviation:", round(x$sigma.e, digits), "\n")

    cat("\n---- Latent process ----\n")
    cat("\nMixture probabilities:\n")
    printCoefmat(round(x$pg, digits), P.values=FALSE, has.Pvalue=FALSE,
                 signif.stars = FALSE, signif.legend = FALSE)

  }else{ #TV

    cat("\nTime-Varying Random Coefficients:\n")
    printCoefmat(round(x$ranTV, digits), P.values=TRUE, has.Pvalue=TRUE,
                 signif.stars = TRUE, signif.legend = TRUE)

    cat("\nResidual scale parameter:", round(x$scale, digits), "- Residual standard deviation:", round(x$sigma.e, digits), "\n")

    cat("\n---- Latent process ----\n")

    cat("\nInitial probabilities:\n")
    printCoefmat(round(x$delta, digits), P.values=FALSE, has.Pvalue=FALSE,
                 signif.stars = FALSE, signif.legend = FALSE)

    cat("\nTransition probabilities:\n")
    printCoefmat(round(x$Gamma, digits), P.values=FALSE, has.Pvalue=FALSE,
                 signif.stars = FALSE, signif.legend = FALSE)

  }

  cat("\nLog-likelihood at convergence:", round(x$lk, digits))
  cat("\nNumber of observations:", x$nobs, "- Number of subjects:", x$nsbjs, "\n")

  if(x$miss == "non-monotone" & x$mod != "TC") message("Data affected by non-monotone missingness: parameter estimates may be biased.")

  invisible(x)
}
