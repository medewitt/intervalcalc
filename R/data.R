#' Serial Interval Data from Nishiura et al
#'
#' @source Supporting materials for Nishiura H, Linton NM, Akhmetzhanov AR 2020
#'     "Serial interval of novel coronavirus (COVID-19) infections"
#'     International journal of infectious diseases 93: 284–6
#'     (<doi:10.1016/j.ijid.2020.02.060>)
#' @format a data.frame with 28 rows and 10 columns
#' \describe{
#'     \item{V1}{Row number}
#'     \item{SIClassification}{Certainty of the Classification}
#'     \item{DiagnosisCounty}{County of Diagnosis}
#'     \item{InfectorOnset}{Onset for index case}
#'     \item{InfecteeOnset}{Onset for infected individual}
#'     \item{ER}{Symptom interval for infector}
#'     \item{EL}{Symptom interval for infector}
#'     \item{SL}{Symptom interval for infectee}
#'     \item{SR}{Symptom interval for infectee}
#'     \item{Source}{Source for the data}
#'
#' }

"nishiura"
