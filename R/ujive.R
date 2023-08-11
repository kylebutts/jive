#' Many IV regressions from Kolesár (2013)
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
#' @return An object of class `manyiv_est` that contains the results.
#' 
#' @export
ujive <- function(
  data, y, exogenous, endogenous, instruments,
  cluster = NULL, ssc = FALSE
) { 
  
  check_args(data, y, exogenous, endogenous, instruments, cluster, ssc)

  # Process formula ------------------------------------------------------------
  fixest::setFixest_fml(
    ..y = y, 
    ..W = exogenous,
    ..T = endogenous, 
    ..Z = instruments
  )
  fml_full = fixest::xpd(..y ~ ..W | ..T ~ ..Z)
  fml_parts = get_fml_parts(fml_full)

  # Compute estimate -----------------------------------------------------------

  # TODO: Use only set of complete variables
  # To make sure full set of NAs are removed
  # if (is.null(est_ZW$obs_selection$obsRemoved)) {
  #   obs_to_keep = 1:nrow(data)
  # } else {
  #   obs_to_keep = est_ZW$obs_selection$obsRemoved
  # }
  est_ZW = fixest::feols(
    fixest::xpd(c(..y, ..T) ~ .[fml_parts$W_lin] + .[fml_parts$Z_lin] | .[fml_parts$W_FE] + .[fml_parts$Z_FE]), 
    data = data
  )

  T = stats::model.matrix(est_ZW[[2]], "lhs", as.matrix = TRUE)
  Y = stats::model.matrix(est_ZW[[1]], "lhs", as.matrix = TRUE)
  n = est_ZW[[2]]$nobs
  In = Matrix::Diagonal(n)
  D_ZW = Matrix::Diagonal(n, stats::hatvalues(est_ZW[[2]]))
  H_ZW_T = stats::predict(est_ZW[[2]])

  # First-stage fitted values
  That = Matrix::solve(In - D_ZW, H_ZW_T - (D_ZW %*% T))

  est_W = fixest::feols(
    fixest::xpd(c(..y, ..T) ~ .[fml_parts$W_lin] | .[fml_parts$W_FE]), 
    data = data
  )
  D_W = Matrix::Diagonal(n, stats::hatvalues(est_W[[2]]))
  H_W_T = stats::predict(est_W[[2]])

  Phat = That - Matrix::solve(In - D_W,  H_W_T - (D_W %*% T))

  # Point estimate
  est = Matrix::crossprod(Phat, Y) / Matrix::crossprod(Phat, T)

  # Standard error -------------------------------------------------------------
  epsilon = stats::resid(est_W[[1]]) - stats::resid(est_W[[2]]) * as.numeric(est)
  se = sqrt(sum(Phat^2 * epsilon^2)) / sum(Phat * T)

  # Small-sample correction
  L = est_W[[2]]$nparams 
  K = est_ZW[[2]]$nparams - L
  if (ssc) {
    # This matches ghpk-metrics/stata-manyiv but I'm not sure it's correct
    se = sqrt(n / (n - L - 1)) * se
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

  out = list(
    beta = as.numeric(est), se = se,
    F = F, Omega = Omega, Xi = Xi, Sargan = Sargan, CD = CD
  )
  class(out) <- c("UJIVE", "jive_est")
  return(out)
}
