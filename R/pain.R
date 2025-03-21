#' Pain Data
#'
#' The pain data frame consists of a total of 357 rows and 4 columns providing information on pain levels of 83 women in labor, followed for up 6 measurement occasions
#'
#' @docType data
#'
#' @importFrom Rdpack reprompt
#'
#' @usage data(pain)
#'
#' @details
#' The data set consists of repeated measurements of self-reported pain on n = 83 women.
#' 43 women were randomly assigned to a pain medication group and 40 to a placebo group.
#' The response was measured every 30 minutes on a 100-mm line: 0 means no pain and 100 means extreme pain.
#' The number of measurements for each woman varies from 1 to 6.
#' Data are severely skewed, and the skewness changes magnitude, and even sign, over time.
#'
#' @references
#' Davis, Charles S. "Semi-parametric and non-parametric methods for the analysis of repeated measurements with applications to clinical trials." Statistics in medicine 10.12 (1991): 1959-1980.
#'
#' @format
#' A data frame with 357 observations on the following 5 variables:
#'
#' \describe{
#'   \item{\code{id}}{woman id}
#'   \item{\code{meas}}{a numeric vector of self-reported pain scores on a 100mm line}
#'   \item{\code{trt}}{a dummy variable with values 1 for subjects who received a pain medication and 0 for subjects who received a placebo}
#'   \item{\code{time}}{a numeric vector of times (minutes since randomization) at which pain was measured}
#' }
#'
"pain"
