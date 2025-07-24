#' @title kasp_dat
#' @description A sample KASP genotype data for marker QC visualization.
#' @docType data
#' @usage data(kasp_dat)
#' @format A data frame with 1344 observations and 11 variables:
#' \describe{
#'   \item{\code{DaughterPlate}}{\code{character} KASP Daughter Plate ID.}
#'   \item{\code{MasterPlate}}{\code{character} KASP Master Plate ID.}
#'   \item{\code{MasterWell}}{\code{character} KASP Master Well ID.}
#'   \item{\code{Call}}{\code{character} KASP observed genotype calls.}
#'   \item{\code{X}}{\code{double} FAM fluorescence values.}
#'   \item{\code{Y}}{\code{double} HEX fluorescence values.}
#'   \item{\code{SNPID}}{\code{character} KASP SNP ID.}
#'   \item{\code{SubjectID}}{\code{character} KASP Subject ID.}
#'   \item{\code{CustomerID}}{\code{character} Customer unique ID for samples.}
#'   \item{\code{DaughterWell}}{\code{character} KASP Daughter Well ID.}
#'   \item{\code{Group}}{\code{character} Predicted genotype for positive controls.}
#'
#' }
"kasp_dat"
