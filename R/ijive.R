#' Estimate improved jackknife IV regression from Ackerberg and Devereux (2009)
#'
#' @description 
#' Estimate the improved JIVE estimator from Ackerberg and Devereux (2009)
#' or the clustered version from Frandsen, Leslie, and McIntyre (2023).
#' Details can be found in 
#' https://direct.mit.edu/rest/article-abstract/91/2/351/57771/Improved-JIVE-Estimators-for-Overidentified-Linear 
#' and 
#' https://sites.google.com/view/emilycleslie/research?authuser=0.
#' 
#' @inheritParams jive
#' @param return_leniency Logical. Should the leave-out fitted value from the first-stage be returned? Default is `FALSE`.
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
#' \item{That}{Leave-out fitted value from the first-stage. This is the fitted values of M_W * T on M_W * Z. Only returned if `return_leniency` is `TRUE`.}
#' 
#' @export
ijive <- function(
  formula, data, cluster = NULL, ssc = FALSE, lo_cluster = !is.null(cluster), return_leniency = FALSE
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

  # TODO: No covariates !!
  # TODO: Use only set of complete observations

  est_ZW = fixest::feols(
    fixest::xpd(c(.[fml_parts$y_fml], .[fml_parts$T_fml]) ~ .[fml_parts$W_lin] + .[fml_parts$Z_lin] | .[fml_parts$W_FE] + .[fml_parts$Z_FE]), 
    data = data
  )
  est_W = fixest::feols(
    fixest::xpd(c(.[fml_parts$y_fml], .[fml_parts$T_fml]) ~ .[fml_parts$W_lin] | .[fml_parts$W_FE]), 
    data = data
  )

  # Y = stats::model.matrix(est_ZW[[1]], "lhs", as.matrix = TRUE)
  # T = stats::model.matrix(est_ZW[[2]], "lhs", as.matrix = TRUE)
  Y_tilde = stats::resid(est_W[[1]])
  T_tilde = stats::resid(est_W[[2]])
  n = est_ZW[[2]]$nobs
  In = Matrix::Diagonal(n)
  
  # https://en.wikipedia.org/wiki/Projection_matrix#Blockwise_formula
  # Trick: 
  # H_{M_W Z} M_W T = (H_{[W Z]} - H_{W}) M_W T
  #                 = H_{[W Z]} M_W T
  #                 = H_{[W Z]} T - H_{[W Z]} H_W T
  #                 = H_{[W Z]} T - H_W T
  H_Ztilde_Ttilde = stats::predict(est_ZW[[2]]) - stats::predict(est_W[[2]])

  # From https://en.wikipedia.org/wiki/Projection_matrix#Blockwise_formula
  # D_{M_W Z} = D_{[W Z]} - D_{W}
  if (lo_cluster == TRUE) {
    # D_Ztilde * Ttilde = D_{[W Z]} Ttilde - D_{W} Ttilde
    D_Ztilde_Ttilde = 
      mult_DX_T(est_ZW[[2]], T_tilde, cl) -
      mult_DX_T(est_W[[2]], T_tilde, cl)

    # Calculate (I - D_{ZW})^{-1} block by block and multiply by H_ZW_T - D_ZW_T
    Phat = solve_ImDX1pDX2_T(
      est_ZW[[2]], est_W[[2]], H_Ztilde_Ttilde - D_Ztilde_Ttilde, cl
    )
    # That = Phat

  } else {
    D_ZW = block_diag_hatvalues(est_ZW[[2]])
    D_W  = block_diag_hatvalues(est_W[[2]])
    D_Ztilde = D_ZW - D_W
    D_Ztilde_Ttilde = D_Ztilde %*% T_tilde
    
    # First-stage fitted values
    Phat = Matrix::solve(In - D_ZW, H_Ztilde_Ttilde - D_Ztilde_Ttilde)
    # That = Phat
  }
  
  # Compute estimate
  est = Matrix::crossprod(Phat, Y_tilde) / Matrix::crossprod(Phat, T_tilde)

  # Standard error -------------------------------------------------------------
  # Y_tilde - T_tilde * \beta
  epsilon = stats::resid(est_W[[1]]) - stats::resid(est_W[[2]]) * as.numeric(est)

  if (is.null(cluster)) {
    se = sqrt(sum(Phat^2 * epsilon^2)) / sum(Phat * T_tilde)
  } else {
    se_cl = lapply(
      split(1:length(cl), cl), 
      function(cl_idx) {
        sum((epsilon[cl_idx] * Phat[cl_idx])^2)
      }
    )
    se = sqrt(Reduce("+", se_cl)) / sum(Phat * T_tilde)
  }

  # Small-sample correction
  L = est_W[[2]]$nparams 
  K = est_ZW[[2]]$nparams - L
  if (!is.null(cluster)) G = max(cl)
  if (ssc) {
    if (is.null(cluster)) {
      se = sqrt(n / (n - L - 1)) * se
    } else {
      se = sqrt(((n - 1) / (n - L - 1)) * (G / (G - 1))) * se
    }
  }
  
  # F-stat ---------------------------------------------------------------------

  # [y_⊥ T_⊥]' [y_⊥ T_⊥] = [y T]' M_W [y T]
  YY = Matrix::crossprod(stats::resid(est_W[1:2]))
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
    clustered = !is.null(cluster),
    n = n, n_instruments = K, n_covariates = L
  )
  if (!is.null(cluster)) out$n_cluster = G
  if (return_leniency) out$That = as.numeric(Phat)
  class(out) <- c("IJIVE", "jive_est")
  return(out)
}
