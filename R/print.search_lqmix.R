#' Print a \code{search_lqmix} object
#'
#' Print an object of \code{\link{class}} \code{search_lqmix}
#'
#'
#' @param x a \code{search_lqmix} object
#' @param digits a non-null value for digits specifying the minimum number of significant digits to be printed
#' @param ... not used
#'
#' @return Return a \code{search_lqmix} object
#'
#' @export


print.search_lqmix = function(x, digits = max(3, getOption("digits") -3), ...){

    oo = x$optimal
    oo$call = NULL

    # if(!is.null(oo$se.scale)) print(summary(oo))
    # else
    print(oo)

  #if(oo$miss == "non-monotone" & oo$mod != "TC") message("Data affected by non-monotone missingness: parameter estimates may be biased.")
  invisible(x)
}
