#' saeHB : Small Area Estimation using Hierarchical Bayesian Method
#'
#' Provides several functions for area level of small area estimation using hierarchical Bayesian (HB) method with Univariate Normal distribution and Univariate Beta distribution for variables of interest. Some dataset produced by a data generation are also provided. The 'rjags' package is employed to obtain parameter estimates. Model-based estimators involves the HB estimators which include the mean and the variation of mean. For the reference, see Rao and Molina (2015) <doi:10.1002/9781118735855>.
#'
#' @section Author(s):
#' Azka Ubaidillah \email{azka@@stis.ac.id}, Ika Yuni Wulansari \email{ikayuni@@stis.ac.id} and Zaza Yuda Perwira \email{221710086@@stis.ac.id}
#'
#' \strong{Maintainer}: Zaza Yuda Perwira \email{221710086@@stis.ac.id}
#'
#'
#' @section Functions:
#' \describe{
#'   \item{\code{Beta}}{Produces HB estimators, standard error, random effect variance, coefficient and plot under beta distribution}
#'   \item{\code{Normal}}{Produces HB estimators, standard error, random effect variance, coefficient and plot under normal distribution}
#'}
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