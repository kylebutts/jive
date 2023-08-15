#' Estimate unbiased jackknife IV regression from Kolesár (2013)
#' 
#' @description 
#' This routine is for a single endogenous regressor following Kolesár (2013). 
#' The function produces a set of estimates using: OLS, two-stage least squares 
#' estimator (TSLS), the limited information maximum likelihood (LIML), 
#' the modified bias-correct 2SLS estimator (MB2SLS), the jackknife IV
#' estimator (JIVE), the unbiased jackknife IV estimator (UJIVE), and the 
#' reverse two-stage least squares estimator (RTSLS). Details can be found in 
#' https://www.princeton.edu/~mkolesar/papers/late_estimation.pdf.
#' 
#' For just UJIVE, see the dedicated function `ujive`.
#' 
#' @inheritParams jive
#' 
#' @return An object of class `jive_est` that contains the results. It contains the following information:
#' \item{beta}{Coefficient on the endogenous variable}
#' \item{se}{Standard error for the coefficient on the endogenous variable}
#' \item{F}{First-stage F statistic}
#' \item{Omega}{Estimate of the covariance matrix of reduced-form errors}
#' \item{Xi}{Estimate of the covariance matrix of the reduced form coefficients. See Kolesar (2012)}
#' \item{Sargan}{Sargan (1958) test for overidentifying restrictions. This is an empty list if only a single instrument is used.}
#' \item{CD}{Cragg and Donald (1993, 1997) test for overidentifying restrictions. This is an empty list if only a single instrument is used.}
#' \item{clustered}{Logical variable indicating whether standard errors are clustered.}
#' \item{n_clustered}{The number of clusters.}
#' \item{n}{The number of observations used.}
#' \item{n_instruments}{The number of non-collinear instruments included.}
#' \item{n_covariates}{The number of non-collinear covariates included.}
#' 
#' @export
ujive <- function(
  formula, data, cluster = NULL, ssc = FALSE
) { 
  
  # `formula` comes first, but flip if needed
  if (inherits(formula, "data.frame")) {
    tmp = formula
    formula = data
    data = tmp
  }

  check_args(data, formula, cluster, ssc)
  fml_parts = get_fml_parts(formula)

  if (!is.null(cluster)) {
    if (inherits(cluster, "character")) {
      cl = data[[cluster]]
    } else {
      cl = eval(cluster[[2]], data)
    }
    cl = as.numeric(as.factor(cl))
  }


  # Compute estimate -----------------------------------------------------------

  # TODO: Use only set of complete variables
  # To make sure full set of NAs are removed
  # if (is.null(est_ZW$obs_selection$obsRemoved)) {
  #   obs_to_keep = 1:nrow(data)
  # } else {
  #   obs_to_keep = est_ZW$obs_selection$obsRemoved
  # }
  est_ZW = fixest::feols(
    fixest::xpd(c(.[fml_parts$y_fml], .[fml_parts$T_fml]) ~ .[fml_parts$W_lin] + .[fml_parts$Z_lin] | .[fml_parts$W_FE] + .[fml_parts$Z_FE]), 
    data = data
  )

  Y = stats::model.matrix(est_ZW[[1]], "lhs", as.matrix = TRUE)
  T = stats::model.matrix(est_ZW[[2]], "lhs", as.matrix = TRUE)
  n = est_ZW[[2]]$nobs
  In = Matrix::Diagonal(n)
  D_ZW = block_diag_hatvalues(est_ZW[[2]])
  H_ZW_T = stats::predict(est_ZW[[2]])

  # First-stage fitted values
  That = Matrix::solve(In - D_ZW, H_ZW_T - (D_ZW %*% T))

  est_W = fixest::feols(
    fixest::xpd(c(.[fml_parts$y_fml], .[fml_parts$T_fml]) ~ .[fml_parts$W_lin] | .[fml_parts$W_FE]), 
    data = data
  )
  D_W = block_diag_hatvalues(est_W[[2]])
  H_W_T = stats::predict(est_W[[2]])

  Phat = That - Matrix::solve(In - D_W,  H_W_T - (D_W %*% T))

  # Point estimate
  est = Matrix::crossprod(Phat, Y) / Matrix::crossprod(Phat, T)

  # Standard error -------------------------------------------------------------
  epsilon = stats::resid(est_W[[1]]) - stats::resid(est_W[[2]]) * as.numeric(est)
  
  if (is.null(cluster)) {
    se = sqrt(sum(Phat^2 * epsilon^2)) / sum(Phat * T)
  } else {
    se_cl = lapply(
      split(1:length(cl), cl), 
      function(cl_idx) {
        sum((epsilon[cl_idx] * Phat[cl_idx])^2)
      }
    )
    se = sqrt(Reduce("+", se_cl)) / sum(Phat * T)
  }

  # Small-sample correction
  L = est_W[[2]]$nparams 
  K = est_ZW[[2]]$nparams - L
  G = max(cl)
  if (ssc) {
    if (is.null(cluster)) {
      se = sqrt(n / (n - L - 1)) * se
    } else {
      se = sqrt(((n - 1) / (n - L - 1)) * (G / (G - 1))) * se
    }
  }

  # F-stat ---------------------------------------------------------------------

  # [y_⊥ T_⊥]' [y_⊥ T_⊥] = [y T]' M_W [y T]
  YY = Matrix::crossprod(stats::resid(est_W))
  YMY = Matrix::crossprod(stats::resid(est_ZW))
  YPY = YY - YMY

  Sp = YMY / (n - K - L)
  S = YPY / n
  eigenvalues = eigen(Matrix::solve(Sp, S))$values
  mmin = min(eigenvalues)

  # first-stage F-statistic
  F = YPY[2, 2] / (K * Sp[2, 2])

  # Estimate of the covariance matrix of reduced-form errors
  Omega = Sp
  
  # Estimate of the covariance matrix of the reduced form coefficients
  Xi = YPY / n - (K / n) * Sp

  # Sargan and C
  Sargan = list()
  CD = list()
  if (K > 1) {
    Sargan$statistic = n * mmin / (1 - K / n - L / n + mmin)
    Sargan$pvalue = 1 - stats::pchisq(Sargan$statistic, K - 1)

    CD$statistic = n * mmin
    CD$pvalue = 1 - stats::pnorm(
      sqrt((n - K - L) / (n - L)) *
      stats::qnorm(stats::pchisq(CD$statistic, K - 1))
    )
  }

  beta = as.numeric(est)
  names(beta) <- all.vars(fml_parts$T_fml)[1]
  names(se) = all.vars(fml_parts$T_fml)[1]

  out = list(
    beta = beta, se = se,
    F = F, Omega = Omega, Xi = Xi, Sargan = Sargan, CD = CD,
    clustered = !is.null(cluster), n_cluster = G,
    n = n, n_instruments = K, n_covariates = L
  )
  class(out) <- c("UJIVE", "jive_est")
  return(out)
}
