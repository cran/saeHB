#' saeHB : Small Area Estimation using Hierarchical Bayesian Method
#'
#' Provides several functions for area level of small area estimation using hierarchical Bayesian (HB) method with several univariate distributions for variable of interest. The dataset that used in every function is generated accordingly in the Example. The 'rjags' package is employed to obtain parameter estimates. Model-based estimators involves the HB estimators which include the mean and the variation of mean. For the reference, see Rao and Molina (2015) <doi:10.1002/9781118735855>.
#'
#' @section Author(s):
#' Azka Ubaidillah \email{azka@@stis.ac.id}, Ika Yuni Wulansari \email{ikayuni@@stis.ac.id}, Zaza Yuda Perwira \email{221710086@@stis.ac.id}
#'
#' \strong{Maintainer}: Zaza Yuda Perwira \email{221710086@@stis.ac.id}
#'
#'
#' @section Functions:
#' \describe{
#'   \item{\code{Beta}}{Produces HB estimators, standard error, random effect variance, coefficient and plot under beta distribution}
#'   \item{\code{Binomial}}{Produces HB estimators, standard error, random effect variance, coefficient and plot under binomial distribution}
#'   \item{\code{Exponential}}{Produces HB estimators, standard error, random effect variance, coefficient and plot under exponential distribution}
#'   \item{\code{ExponentialDouble}}{Produces HB estimators, standard error, random effect variance, coefficient and plot double exponential distribution}
#'   \item{\code{Gamma}}{Produces HB estimators, standard error, random effect variance, coefficient and plot under Gamma distribution}
#'   \item{\code{Logistic}}{Produces HB estimators, standard error, random effect variance, coefficient and plot under logistic distribution}
#'   \item{\code{Lognormal}}{Produces HB estimators, standard error, random effect variance, coefficient and plot under lognormal distribution}
#'   \item{\code{NegativeBinomial}}{Produces HB estimators, standard error, random effect variance, coefficient and plot under negative binomial distribution}
#'   \item{\code{Normal}}{Produces HB estimators, standard error, random effect variance, coefficient and plot under normal distribution}
#'   \item{\code{Poisson}}{Produces HB estimators, standard error, random effect variance, coefficient and plot under poisson distribution}
#'   \item{\code{PoissonGamma}}{Produces HB estimators, standard error, random effect variance, coefficient and plot under poisson gamma distribution}
#'   \item{\code{Student_t}}{Produces HB estimators, standard error, random effect variance, coefficient and plot under student-t distribution}
#'   \item{\code{Student_tnc}}{Produces HB estimators, standard error, random effect variance, coefficient and plot under student-t (not central) distribution}
#'   \item{\code{Weibull}}{Produces HB estimators, standard error, random effect variance, coefficient and plot under weibull distribution}
#'}
#'
#' @section Reference:
#' \itemize{
#'   \item{Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd Edition. New York: John Wiley and Sons, Inc. <doi:10.1002/9781118735855>.}
#' }
#'
#'
#' @docType package
#' @name saeHB
#'
#' @import rjags
#' @import coda
#' @import stringr
#' @import stats
#' @import grDevices
#' @import graphics

NULL
