#' CD4 Data
#'
#' The cd4 data frame is made by a total of 2376 rows and 8 columns providing information on
#' CD4 cell counts of 369 subjects followed for a maximum of 12 measurement occasions.
#'
#' @docType data
#'
#' @importFrom Rdpack reprompt
#'
#' @usage data(cd4)
#'
#' @details
#' Multi-center AIDS Cohort Study providing a total of 2376 CD4+ cell counts of 369 HIV-infected men covering a period of approximately eight and half years.
#' The number of measurements for each subject varies from 1 to 12. The CD4+ cell data are highly unbalanced.
#'
#' @references{
#'   \insertRef{ref:cd4}{lqmix}
#' }
#'
#' @format
#' A data frame with 2376 observations on the following 8 variables:
#'
#' \describe{
#'   \item{\code{sbj.id}}{subject id}
#'   \item{\code{time.id}}{time id}
#'   \item{\code{count}}{CD4 count}
#'   \item{\code{lcount}}{log(CD4 count + 1)}
#'   \item{\code{time}}{years since seroconversion}
#'   \item{\code{age}}{age (yrs) centered around 30}
#'   \item{\code{packs}}{packs of cigarettes per day}
#'   \item{\code{partners}}{number of sexual partners}
#'   \item{\code{drugs}}{recreational drug use indicator}
#'   \item{\code{cesd}}{depression score}
#' }
#'
"cd4"
