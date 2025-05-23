% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ijive.R
\name{ijive}
\alias{ijive}
\title{Estimate improved jackknife IV regression from Ackerberg and Devereux (2009)}
\usage{
ijive(
  formula,
  data,
  cluster = NULL,
  ssc = FALSE,
  lo_cluster = !is.null(cluster),
  return_leniency = FALSE
)
}
\arguments{
\item{formula}{Formula. following the syntax of the \code{fixest} package. In short, \code{y ~ exo | exo_FEs | endo ~ instrument | instrument_FEs}.}

\item{data}{Data.frame}

\item{cluster}{Character or formula. The cluster variable. For non-clustered robust standard errors, set to NULL.}

\item{ssc}{Logical. Should a small sample adjustment be made? Default is \code{TRUE}.}

\item{lo_cluster}{Logical. Should leave-out cluster be used for estimation? Default is \code{TRUE} if \code{cluster} is not NULL.}

\item{return_leniency}{Logical. Should the leave-out fitted value from the first-stage be returned? Default is \code{FALSE}.}
}
\value{
An object of class \code{jive_est} that contains the results. It contains the following information:
\item{beta}{Coefficient on the endogenous variable}
\item{se}{Standard error for the coefficient on the endogenous variable}
\item{F}{First-stage F statistic}
\item{Omega}{Estimate of the covariance matrix of reduced-form errors}
\item{Xi}{Estimate of the covariance matrix of the reduced form coefficients. See Kolesar (2012)}
\item{Sargan}{Sargan (1958) test for overidentifying restrictions. This is an empty list if only a single instrument is used.}
\item{CD}{Cragg and Donald (1993, 1997) test for overidentifying restrictions. This is an empty list if only a single instrument is used.}
\item{clustered}{Logical variable indicating whether standard errors are clustered.}
\item{n_clustered}{The number of clusters.}
\item{n}{The number of observations used.}
\item{n_instruments}{The number of non-collinear instruments included.}
\item{n_covariates}{The number of non-collinear covariates included.}
\item{That}{Leave-out fitted value from the first-stage. This is the fitted values of M_W * T on M_W * Z. Only returned if \code{return_leniency} is \code{TRUE}.}
}
\description{
Estimate the improved JIVE estimator from Ackerberg and Devereux (2009)
or the clustered version from Frandsen, Leslie, and McIntyre (2023).
Details can be found in
https://direct.mit.edu/rest/article-abstract/91/2/351/57771/Improved-JIVE-Estimators-for-Overidentified-Linear
and
https://sites.google.com/view/emilycleslie/research?authuser=0.
}
