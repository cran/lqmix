#' Print an \code{lqmix} object
#'
#' Print an object of \code{\link{class}} \code{lqmix}
#'
#' @param x an \code{lqmix} object
#' @param digits a non-null value for digits specifying the minimum number of significant digits to be printed
#' @param ... not used
#'
#' @return Return an \code{lqmix} object
#'
#' @export


print.lqmix = function(x, digits = max(3, getOption("digits") -3), ...){

  if(!is.null(x$call)){
    if(x$model == "TC"){
      cat(paste("Model: TC random coefficients with G=", x$G, " at qtl=", x$qtl, sep=""))
      cat("\n")
      cat("******************************************************", "\n")
    }else if(x$model == "TV"){
      cat(paste("Model: TV random coefficients with m=", x$m, " at qtl=", x$qtl, sep=""))
      cat("\n")
      cat("******************************************************", "\n")
    }else{
      cat(paste("Model: TC and TV random coefficients with m=", x$m, " and G=", x$G, " at qtl=", x$qtl, sep=""))
      cat("\n")
      cat("******************************************************************", "\n")
    }
  }else{
    if(x$model == "TC"){
      cat(paste("Opt model: TC random coefficients with G=", x$G, " at qtl=", x$qtl, sep=""))
      cat("\n")
      cat("********************************************************", "\n")
    }else if(x$model == "TV"){
      cat(paste("Opt model: TV random coefficients with m=", x$m, " at qtl=", x$qtl, sep=""))
      cat("\n")
      cat("**********************************************************", "\n")
    }else{
      cat(paste("Opt model: TC and TV random coefficients with m=", x$m, " and G=", x$G, " at qtl=", x$qtl, sep=""))
      cat("\n")
      cat("*********************************************************************", "\n")
    }
  }

  cat("\n---- Observed process ----\n")
  if(!is.null(x$betaf)){
    cat("\nFixed Coefficients:\n")
    print(round(x$betaf, digits))

  }

  if(!is.null(x$betarTC) & !is.null(x$betarTV)){ #TCTV
    cat("\nTime-Constant Random Coefficients:\n")
    rownames(x$betarTC) = paste("Comp", 1:nrow(x$betarTC), sep="")
    print(round(x$betarTC, digits))
    cat("\nTime-Varying Random Coefficients:\n")
    rownames(x$betarTV) = paste("St", 1:nrow(x$betarTV), sep="")
    print(round(x$betarTV, digits))

    cat("\nResidual scale parameter:", round(x$scale, digits), "- Residual standard deviation:", round(x$sigma.e, digits), "\n")

    cat("\n---- Latent process ----\n")
    cat("\nMixture probabilities:\n")
    names(x$pg) = paste("Comp", 1:nrow(x$betarTC), sep="")
    print(round(x$pg, digits))

    cat("\nInitial probabilities:\n")
    names(x$delta) = paste("St", 1:nrow(x$betarTV), sep="")
    print(round(x$delta, digits))

    cat("\nTransition probabilities:\n")
    rownames(x$Gamma) = paste("fromSt", 1:nrow(x$betarTV), sep="")
    colnames(x$Gamma) = paste("toSt", 1:nrow(x$betarTV), sep="")
    print(round(x$Gamma, digits))

  }else if (!is.null(x$betarTC) & is.null(x$betarTV)){ # TC

    cat("\nTime-Constant Random Coefficients:\n")
    rownames(x$betarTC) = paste("Comp", 1:nrow(x$betarTC), sep="")
    print(round(x$betarTC,digits))

    cat("\nResidual scale parameter:", round(x$scale, digits), "- Residual standard deviation:", round(x$sigma.e, digits), "\n")

    cat("\n---- Latent process ----\n")
    cat("\nMixture probabilities:\n")
    names(x$pg) = paste("Comp", 1:nrow(x$betarTC), sep="")
    print(round(x$pg, digits))

  }else{ # TV
    cat("\nTime-Varying Random Coefficients:\n")
    rownames(x$betarTV) = paste("St", 1:nrow(x$betarTV), sep="")
    print(round(x$betarTV,digits))

    cat("\nResidual scale parameter:", round(x$scale, digits), "- Residual standard deviation:", round(x$sigma.e, digits), "\n")

    cat("\n---- Latent process ----\n")
    cat("\nInitial probabilities:\n")
    names(x$delta) = paste("St", 1:nrow(x$betarTV), sep="")
    print(round(x$delta, digits))

    cat("\nTransition probabilities:\n")
    rownames(x$Gamma) = paste("fromSt", 1:nrow(x$betarTV), sep="")
    colnames(x$Gamma) = paste("toSt", 1:nrow(x$betarTV), sep="")
    print(round(x$Gamma, digits))
  }

  cat("\nLog-likelihood at convergence:", round(x$lk, digits))
  cat("\nNumber of observations:", x$nobs, "- Number of subjects:", x$nsbjs, "\n")

  if(x$miss == "non-monotone" & x$mod != "TC") message("Data affected by non-monotone missingness: parameter estimates may be biased.")

  invisible(x)
}
