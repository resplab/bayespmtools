#' Isaric Dataset
#' 
#' Data from the International Severe Acute Respiratory and Emerging Infection Consortium regarding
#' Regions in the UK.
#' 
#' @name isaric
#' @format A data frame with 8 rows and 10 columns
#' \describe{
#'    \item{Region}{Region where the sample was drawn}
#'    \item{Sample_Size}{Raw number of total subjects available in the region's dataset}
#'    \item{n}{Number of subjects used in analysis after exclusions}
#'    \item{n_events}{Number of positive subjects}
#'    \item{cstat}{C-statistic}
#'    \item{cstat_l}{Lower bound for the confidence interval of the C-statistic}
#'    \item{cal_mean}{Calibraiton Mean}
#'    \item{cal_mean_l}{Lower bound for the confidence interval of the calibration mean}
#'    \item{cal_slope}{Calibration slope}
#'    \item{cal_slope_l}{Lower bound of the confidence interval of the calibration slope}
#' }
#' 
#' @source Simulated Data
#' @keywords datasets
"isaric"