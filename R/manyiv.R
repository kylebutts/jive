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
manyiv <- function(
  data, y, exogenous, endogenous, instruments,
  cluster = NULL, ssc = FALSE
) { 
  # I could do y + exogenous as one & endogenous + instruments as another,
  # but I think that might mislead folks into thinking there are two regressions
  # occuring
  # Check arguments ------------------------------------------------------------
  dreamerr::check_arg(data, "data.frame")
  dreamerr::check_arg(y, cluster, "character var(data) | os formula var(data) right(1, 1)", .data = data)
  dreamerr::check_arg(
    endogenous,
    "os formula var(data) right(1, 1)",
    .data = data
  )
  dreamerr::check_arg(
    exogenous, instruments, 
    "os formula var(data) right(1, 2)",
    .data = data
  )
  dreamerr::check_arg(ssc, "logical scalar")

  # Get model matrices ---------------------------------------------------------
  fixest::setFixest_fml(
    ..y = y, 
    ..W = exogenous,
    ..T = endogenous, 
    ..Z = instruments
  )

  # Get each part of the formula into it's own part
  fml_full = fixest::xpd(..y ~ ..W | ..T ~ ..Z)
  fml_parts = get_fml_parts(fml_full)

  exo = fixest::feols(
    fixest::xpd(..y ~ .[fml_parts$W_lin] | .[fml_parts$W_FE]), 
    data = data
  )
  endo = fixest::feols(
    fixest::xpd(..T ~ .[fml_parts$Z_lin] | .[fml_parts$W_FE] + .[fml_parts$Z_FE]), 
    data = data
  )

  # y, W, and FE
  exo_smm = sparse_model_matrix(exo, data = data, type = c("lhs", "rhs", "fixef"), combine = FALSE)
  # T, Z, and IV_FE
  endo_smm = sparse_model_matrix(endo, data = data, type = c("lhs", "rhs", "fixef"), combine = FALSE)

  # Correctly drop  NAs
  # I've checked that this is robust to exo_na and endo_na being NULL
  # `union`, `setdiff` pass through nulls and `length(NULL) == 0`.
  exo_na = exo$obs_selection$obsRemoved
  endo_na = endo$obs_selection$obsRemoved
  to_remove = union(exo_na, endo_na)
  if (!is.null(exo_na)) {
    remaining_rows_exo = (1:nrow(data))[exo_na]
    yet_to_remove_exo = setdiff(endo_na, exo_na)
    if (length(yet_to_remove_exo) > 0) {
      to_remove_exo = which(remaining_rows_exo %in% (-1 * yet_to_remove_exo))
      exo_smm = lapply(exo_smm, function(x) x[-to_remove_exo, , drop = FALSE])
    }
  }
  if (!is.null(endo_na)) {
    remaining_rows_endo = (1:nrow(data))[endo_na]
    yet_to_remove_endo = setdiff(exo_na, exo_na)
    if (length(yet_to_remove_endo) > 0) {
      to_remove_endo = which(remaining_rows_endo %in% (-1 * yet_to_remove_endo))
      endo_smm = lapply(endo_smm, function(x) x[-to_remove_endo, , drop = FALSE])
    }
  }
  
  # Define matrices
  if (!is.null(cluster)) {
    if (inherits(cluster, "formula")) cluster = as.character(cluster)[2]
    cl = data[[cluster]]
    if (!is.null(to_remove)) cl[to_remove]
  }
  
  y = as.matrix(exo_smm$lhs)
  W = as.matrix(exo_smm$rhs)
  W_FE = exo_smm$fixef
  T = as.matrix(endo_smm$lhs)
  Z = as.matrix(endo_smm$rhs)
  Z_FE = endo_smm$fixef
  # rm(list(exo_smm, endo_smm))

  # Demean by fixed-effects ----------------------------------------------------

  # TODO: drop singletons (ugh)

  # Get FE variable for use in `demean`
  data_FE = data_IVFE = NULL

  # If no fixed-effects at all, remove second intercept column
  if (is.null(W_FE) & is.null(Z_FE)) {
    intercept_col = ("(Intercept)" == colnames(Z))
    if (any(intercept_col)) {
      Z = Z[, -intercept_col, drop = FALSE]
    }
  }


  if (!is.null(W_FE)) {
    data_FE = stats::model.matrix(exo, type = "fixef")
    yp = drop_0_cols(fixest::demean(y, data_FE))
    Tp = drop_0_cols(fixest::demean(T, data_FE))
    Wp = drop_0_cols(fixest::demean(W, data_FE))
  } else {
    yp = y
    Tp = T
    Wp = W
  }
  if (!is.null(Z_FE)) {
    data_IVFE = stats::model.matrix(endo, type = "fixef")
    yq = drop_0_cols(fixest::demean(y, data_IVFE))
    Tq = drop_0_cols(fixest::demean(T, data_IVFE))
    Wq = drop_0_cols(fixest::demean(W, data_IVFE))
    Zq = drop_0_cols(fixest::demean(Z, data_IVFE))
  } else {
    yq = yp
    Tq = Tp
    Wq = Wp
    
    if (is.null(W_FE)) {
      Zq = Z
    }
  }


  # TODO: 
  # Drop collinear from instruments
  # Get degree of freedom adjustments
  n = nrow(W)
  K = ncol(Z)
  L = ncol(W) 
  if (!is.null(W_FE)) L = L + ncol(W_FE)
  get_dof <- function(exo, endo, W, W_FE, Z, Z_FE) {

  }

  # Zq = drop_0_cols(Zq)
  # Z = drop_0_cols(Z)
  # Z = drop_0_cols(Z) 
  # print(K)
  # print(L)
  # print(n)

  # Point Estimates ------------------------------------------------------------
  # Note: n, K, L already defined

  # [y_⊥ T_⊥] = M_W [y T]
  MW_yT <- annihilator(Wp, cbind(yp, Tp))
  # [y_⊥ T_⊥] = M_W M_D [y T]
  MWD_yT <- annihilator(Wq, cbind(yq, Tq))
  # Z_⊥ = M_W M_D Z
  MWD_Z <- annihilator(Wq, Zq) 

  # [y_⊥ T_⊥]' [y_⊥ T_⊥] = [y T]' M_W [y T]
  YY <- Matrix::crossprod(MW_yT)

  # [solve(Z_⊥, y_⊥) solve(Z_⊥, T_⊥)] = Reduced form and First stage
  RFS = Matrix::solve(Matrix::crossprod(MWD_Z), Matrix::crossprod(MWD_Z, MWD_yT))

  # H_{Z_⊥} [y_⊥ T_⊥]
  # TODO: check ncol(Zq) > 0
  HZ_yT = (MWD_Z %*% RFS) + MW_yT - MWD_yT 

  # [y_⊥ T_⊥]' H_{Z_⊥} [y_⊥ T_⊥]
  YPY <- Matrix::crossprod(MW_yT, HZ_yT)
  # [y_⊥ T_⊥]' M_{Z_⊥} [y_⊥ T_⊥]
  YMY <- YY - YPY 

  ## k-class: OLS, TSLS, LIML, MBTLS -------------------------------------------

  # Note: These are all coded as
  # 
  #     (T_⊥' (I - k W_{Z_⊥}) y_⊥) / (T_⊥' (I - k W_{Z_⊥}) T_⊥)
  # 
  # So different values of k give different estimands.
  # 
  # - k = 0 -> YY[1, 2]  / YY[2, 2]  = (y_⊥' T_⊥) / (T_⊥' T_⊥)
  # - k = 1 -> YPY[1, 2] / YPY[2, 2] = (y_⊥' H_{Z_⊥} T_⊥) / (T_⊥' H_{Z_⊥} T_⊥)
  # - The other two give liml and mbtsls

  k <- c(
    0, 1, 
    min(eigen(mrsolve(YMY, YY))$values), 
    (1 - L / n) / (1 - (K - 1) / n - L / n)
  )
  beta_k <- 
    (YY[1, 2] - k * YMY[1, 2]) / 
    (YY[2, 2] - k * YMY[2, 2])

  # Data Check: Match Stata
  # round(beta_k, 4)

  ## JIVE, UJIVE ---------------------------------------------------------------

  # DW = diag(P_W) as a vector
  DW <- diag_hat(cbind(Wp, W_FE))
  DW <- Matrix::Diagonal(length(DW), DW)

  # (I - D_W)^{-1}
  iIDW <- 1 / (1 - Matrix::diag(DW))
  IDW = (Matrix::Diagonal(nrow(DW)) - DW)

  # DZW = diag(P_ZW) as a vector
  ZW = cbind(Wq, Zq)
  DZW <- diag_hat(cbind(ZW, Z_FE))
  DZW <- Matrix::Diagonal(length(DZW), DZW)

  # (I - D_ZW)^{-1}
  iIDZW <- 1 / (1 - Matrix::diag(DZW))
  IDZW = (Matrix::Diagonal(nrow(DZW)) - DZW)

  # JIVE and UJIVE 
  hatTjive <- (T - iIDW * annihilator(Wp, Tp))
  hatTujive <- T - iIDZW * annihilator(ZW, Tq)
  hatPujive <- hatTujive - hatTjive
  hatPjive <- as.matrix(annihilator(Wp, hatTujive))
  if (!is.null(data_FE)) {
    hatPjive <- fixest::demean(hatPjive, data_FE)
  }



  ## Store point estimates -----------------------------------------------------

  beta_labels <- c(
    'OLS', 
    'TSLS', 
    'LIML', 
    'MBTSLS', 
    'JIVE', 
    'UJIVE', 
    'RTSLS'
  )
  beta <- c(
    # OLS, TSLS, LIML, MBTSLS
    beta_k,
    # JIVE 
    as.numeric(
      (Matrix::t(hatPjive) %*% y) / (Matrix::t(hatPjive) %*% T)
    ),
    # UJIVE
    as.numeric(
      (Matrix::t(hatPujive) %*% y) / (Matrix::t(hatPujive) %*% T)
    ), 
    # RTSLS
    YPY[1, 1] / YPY[1, 2]
  )

  names(beta) <- beta_labels

  # Data Check: Match Stata
  # round(beta, 4)

  # Standard Errors ------------------------------------------------------------

  se = matrix(nrow = 5, ncol = 7)

  # e
  get_epsilon = function(beta, MW_yT) {
    cols = lapply(beta, function(b) {
      MW_yT[, 1, drop = FALSE] - MW_yT[, 2, drop = FALSE] * b
    })
    epsilon = do.call("cbind", cols)
    colnames(epsilon) = colnames(beta)
    return(epsilon)
  }

  ## Homoskedastic -------------------------------------------------------------

  epsilon = get_epsilon(beta[1:6], MW_yT)
  se[1, 1] = sqrt(
    (Matrix::crossprod(epsilon[, 1]) / n) / 
    as.numeric(Matrix::crossprod(MW_yT[, 2, drop = FALSE]))
  )
  se[1, 2] = sqrt(
    (Matrix::crossprod(epsilon[, 2]) / n) / 
    YPY[2, 2]
  )
  se[1, 3] = sqrt(
    (Matrix::crossprod(epsilon[, 3]) / n) / 
    YPY[2, 2]
  )
  se[1, 4] = sqrt(
    (Matrix::crossprod(epsilon[, 4]) / n) / 
    YPY[2, 2]
  )
  se[1, 5] = sqrt(
    (Matrix::crossprod(epsilon[, 5]) / n) * 
    as.numeric(Matrix::crossprod(hatPjive)) / 
    as.numeric(Matrix::crossprod(hatPjive, T))^2
  )
  se[1, 6] = sqrt(
    (Matrix::crossprod(epsilon[, 6]) / n) * 
    as.numeric(Matrix::crossprod(hatPujive)) / 
    as.numeric(Matrix::crossprod(hatPujive, T))^2
  )

  ## Heteroskedastic -----------------------------------------------------------

  hatP = cbind(MW_yT[, 2, drop = FALSE], HZ_yT[, 2], HZ_yT[, 2], HZ_yT[, 2], hatPjive, hatPujive)

  se[2, 1:6] = as.numeric(
    sqrt(Matrix::colSums((epsilon * hatP)^2)) / Matrix::crossprod(T, hatP)
  )
  # UJIVE: hatvalues(P_ujive) .* epsilon = M_W * y - M_w * T * \beta

  ## Clustered, if requested ---------------------------------------------------

  # using quf
  if (!is.null(cluster)) {
    cl = fixest:::quickUnclassFactor(cl)
    G = max(cl)
    sec = matrix(0, nrow = 1, ncol = 6)
    for (i in 1:G) {
      cl_i = which(cl == i)
      hatP_i = hatP[cl_i, ]
      epsilon_i = epsilon[cl_i, ]

      sec = sec + Matrix::colSums((epsilon_i * hatP_i)^2)
    }

    se[3, 1:6] = as.numeric(sqrt(sec) / Matrix::crossprod(T, hatP))
  }

  ## Many Instruments ----------------------------------------------------------

  # Notation
  Sp = YMY / (n - K - L)
  S = YPY / n
  eigenvalues = eigen(Matrix::solve(Sp, S))$values
  mmin = min(eigenvalues)

  # Hessian of random-effects
  lamre = max(eigenvalues) - K / n
  a = as.matrix(c(beta[3], 1), ncol = 1)
  b = as.matrix(c(1, -beta[3]), ncol = 1)

  helper_tsolve <- function(X, Y) {
    Matrix::t(Matrix::solve(
      Y %*% Matrix::t(Y), 
      Y %*% Matrix::t(X)
    ))
  }
  X = Matrix::tcrossprod(a)
  Y = as.numeric(helper_tsolve(t(a), Sp) %*% a)
  Z = lamre * 1 / Y^2 * as.matrix(Matrix::colSums(X * Y), ncol = 1)
  Omre = 
    (n - K - L) * Sp / (n - L) + 
    n / (n - L) * apply(S, 2, function(col) col - Z)

  Omre

  Qs = as.numeric(
    Matrix::crossprod(b, S) %*% b / 
    Matrix::crossprod(b, Omre) %*% b
  )
  c = lamre * Qs / ((1 - L / n) * (K / n + lamre))

  se[4, 3] = as.numeric(sqrt(
    -1 * as.numeric(Matrix::crossprod(b, Omre) %*% b) / 
    (n * lamre) * (lamre + K / n) /
    (
      Qs * Omre[2, 2] - S[2, 2] + 
      c / (1 - c) * Qs / (mrsolve(Omre, Matrix::t(a)) %*% a)
    )
  ))

  # mbtsls, using maximum URE likelihood plug-in estimator
  b = as.matrix(c(1, -1 * beta[4]), ncol = 1)
  Lam11 = max(0, Matrix::t(b) %*% (S - K/n * Sp) %*% b)

  if (mmin > K / n) {
    Lam22 = S[2, 2] - K/n * Sp[2, 2]
    Omure = Sp
  } else {
    Lam22 = as.numeric(lamre / helper_tsolve(t(a), Omre) %*% a)
    Omure = Omre
  }

  Gamma = matrix(c(1, 0, -1 * beta[4], 1), ncol = 2, byrow = TRUE)
  Sig = Matrix::crossprod(Gamma, Omure) %*% Gamma
  h = ((1 - L / n) * (K - 1) / n) / (1 - L / n - (K - 1) / n)

  Vvalid = Sig[1,1] / Lam22 + h * (Sig[1,1] * Sig[2,2] + Sig[1,2]^2) / Lam22^2
  Vinvalid = Vvalid + (Lam11 * Omure[2,2] + Lam11 * Lam22 * n/K) / Lam22^2

  se[4, 4] = sqrt(Vvalid / n)
  se[5, 4] = sqrt(Vinvalid / n)

  ## Small-sample adjustment, if requested -------------------------------------

  qc = if (ssc) sqrt(n / (n - L - 1)) else 1
  se[c(1, 2, 4, 5), ] = se[c(1, 2, 4, 5), ] * qc

  qc = if (ssc) sqrt(((n - 1) / (n - L - 1)) * (G / (G - 1))) else 1
  se[3, ] = se[3, ] * qc

  # Stats ----------------------------------------------------------------------

  # first-stage F
  F = YPY[2, 2] / (K * Sp[2, 2])

  # Xi
  Xi = YPY / n - (K / n) * Sp

  overid = vector("numeric", 2)
  pvalue = vector("numeric", 2)
  if (K > 1) {
    overid[1] = n * mmin / (1 - K / n - L / n + mmin)
    pvalue[1] = 1 - stats::pchisq(overid[1], K - 1)

    overid[2] = n * mmin
    pvalue[2] = 1 - stats::pnorm(
      sqrt((n - K - L) / (n - L)) *
      stats::qnorm(stats::pchisq(overid[2], K - 1))
    )
  }

  stats = list(
    F = F, Omega = Sp, Xi = Xi, 
    Sargan = c(overid[1], pvalue[1]),
    CD = c(overid[2], pvalue[2])
  )

  # Output --------------------------------------------------------------------

  out = list(
    beta = beta,
    se = se,
    F = stats$F,
    omega = stats$Omega,
    xi = stats$Xi,
    sargan = stats$Sargan,
    cd = stats$CD,
    nobs = n,
    n_instruments = K,
    n_covariates = L,
    clustered = !is.null(cluster)
  )
  class(out) <- "manyiv_est"
  
  return(out)
}

print.manyiv_est <- function(out) {
  if (out$clustered) {
    se = out$se[3, ]
    se_label = "Clustered SE"
  } else {
    se = out$se[2, ]
    se_label = "Robust SE"
  }
  beta = out$beta
  tstat = beta / se
  pval = 2 * (1 - stats::pnorm(abs(tstat)))

  est = data.frame(
    Estimator = names(beta),
    Estimate = sprintf("%.4f", beta)
  )
  se = c(sprintf("(%.3f)", se[1:6]), "(.)")
  est[[se_label]] = se

  # Print out nice table to console
  print(est)
  stringmagic::cat_magic(
    "\n{n ? n} observations, {n ? K} instruments, {n ? L} covariates, first-stage F = {%0.3f ? F}\n", 
    n = out$n, L = out$L, K = out$K, F = out$F
  )
}


# Solves xA' = B
mrsolve <- function(A, B) {
  # B A'^{-1} = (A'^{-1} B')'
  Matrix::t(Matrix::solve(Matrix::t(A), Matrix::t(B)))
}

# (I - A (A'A)^-1 A') B
annihilator <- function(A, B) {
  B - A %*% Matrix::solve(Matrix::crossprod(A), Matrix::crossprod(A, B))
}

# diag(X (X'X)^-1 X')
diag_hat <- function(X) {
  Matrix::diag(X %*% Matrix::solve(Matrix::crossprod(X), Matrix::t(X)))
}

# Drop columns of matrix that are all 0s
drop_0_cols <- function(mat) {
  keep_cols = apply(mat, 2, function(x) !all(x == 0))
  mat[, keep_cols, drop = FALSE]
}
