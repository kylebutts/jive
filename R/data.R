#' Judge data from Stevenson (2018)
#'
#' This dataset is from Stevenson (2018). The dataset consists of individual court cases that are quasi-randomly assigned to 8 judges (magistrates) in Philadelphia. The data contains information on pretrial detention and several outcomes, including whether or not a defendant subsequently pleads guilty. 
#'
#' @format ## `stevenson`
#' A data frame with 331,971 rows and 50 columns:
#' \describe{
#'   \item{judge_pre}{Which judge the case was randomly-assigned to. Numbered 1 through 8.}
#'   \item{jail3}{Indicator variable. Equals one if the defendent is assigned to pretrial detention.}
#'   \item{guilt}{Indicator variable. Equals one if the defendant is guilty of at lesat one charge.}
#'   \item{bailDate}{Date of bail hearing.}
#'   \item{priorCases}{Number of prior cases the defendent has.}
#'   \item{prior_felChar}{Number of prior felony cases the defendent has.}
#'   \item{prior_guilt}{Number of prior cases where the defendent was found guilty.}
#'   \item{fel}{Indicator variable. Equals one if it is a felony charge.}
#'   \item{mis}{Indicator variable. Equals one if it is a misdomenor charge.}
#'   \item{robbery}{Indicator variable. Equals one if the case is a robery case.}
#'   \item{aggAss}{Indicator variable. Equals one if the case is a aggrevated assault case.}
#'   \item{possess}{Indicator variable. Equals one if the case is a drug posession case.}
#'   \item{drugSell}{Indicator variable. Equals one if the case is a drug selling case.}
#'   \item{DUI1st}{Indicator variable. Equals one if the case is a DUI first offense case.}
#'   \item{white}{Indicator varible. Equals one if the defendent is white.}
#'   \item{black}{Indicator varible. Equals one if the defendent is black.}
#'   \item{age}{Age of defendent.}
#'   \item{male}{Indicator variable. Equals one if the defendent is male.}
#'   \item{priorWI5}{Indicator variable. Equals one if the defendent has a prior arrest within five-years of the bail-hearing.}
#'   \item{morning}{Trials time of of day. Either morning, evening, grave, or weekend}
#'   \item{t}{The time-period}
#' }
#' @source [Stevenson 2018](https://academic.oup.com/jleo/article/34/4/511/5100740)
"stevenson"
