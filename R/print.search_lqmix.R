#' Print a \code{search_lqmix} object
#'
#' Print an object of \code{\link{class}} \code{\link{search_lqmix}}.
#'
#'
#' @param x a \code{search_lqmix} object
#' @param digits a non-null value for digits specifying the minimum number of significant digits to be printed
#' @param ... not used
#'
#' @return Return a \code{\link{search_lqmix}} object.
#'
#' @export


print.search_lqmix = function(x, digits = max(3, getOption("digits") -3), ...){

    oo = x$optimal
    oo$call = NULL
    print(oo)

  invisible(x)
}
