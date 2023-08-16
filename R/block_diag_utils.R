#' Computes D_X %*% T for diagonal block projection matrix
#' @param model A `fixest` model object.
#' @param T Numeric vector
#' @param cl A vector of cluster membership with the same length as the number
#' of observations used in `model`.
mult_DX_T <- function(model, T, cl = NULL) {
  if (!is.null(cl)) cl_split <- split(1:model$nobs, cl)
  n <- model$nobs

  # Only works for feglm/feols objects
  if (isTRUE(model$lean)) {
    stop("The method 'hatvalues.fixest' cannot be applied to 'lean' fixest objects. Please re-estimate with 'lean = FALSE'.")
  }

  if (!is.null(model$iv)) {
    stop("The method 'hatvalues.fixest' cannot be applied to IV models.")
  }

  method <- model$method_type

  if (!(method %in% c("feols", "feglm"))) {
    stop("'hatvalues' is currently not implemented for function ", method, ".")
  }

  # Outline: Blockwise Formula for Projection Matrix
  # https://en.wikipedia.org/wiki/Projection_matrix#Blockwise_formula
  # Part 1: fixed effects
  # Part 2: slope vars (demeaned by fixed effects)
  # Part 3: rhs vars (demeaned by fixed effects and slope vars)
  # Sum the three parts to get the projection matrix
  P_list <- list()
  mats <- list()
  mats$fixef <- fixest::sparse_model_matrix(
    model,
    type = c("fixef"),
    collin.rm = TRUE, na.rm = TRUE
  )

  # Check for slope.vars and move to seperate matrix
  if (!is.null(model$slope_flag)) {
    slope_var_cols <- which(
      grepl("\\[\\[", colnames(mats$fixef))
    )
    mats$slope_vars <- mats$fixef[, slope_var_cols]
    mats$fixef <- mats$fixef[, -slope_var_cols]
  }

  # Get weights
  if (method == "feols") {
    weights <- model$weights
  } else {
    weights <- model$irls_weights
  }

  if (!is.null(model$fixef_vars)) {
    f <- stats::model.matrix(model, type = "fixef")
    f <- subset(f, select = model$fixef_vars)

    # Matrix of fixed effects times weights
    if (!is.null(weights)) mats$fixef <- mats$fixef * sqrt(weights)

    # Demean X by FEs and slope_vars
    if (is.null(model$onlyFixef) & method == "feols") {
      mats$rhs <- fixest::demean(model)[, -1, drop = FALSE]
      if (!is.null(weights)) mats$rhs <- mats$rhs * sqrt(weights)
    }
    if (is.null(model$onlyFixef) & method == "feglm") {
      mats$rhs <- fixest::demean(
        stats::model.matrix(model, "rhs"),
        f,
        weights = weights
      )
      if (!is.null(weights)) mats$rhs <- mats$rhs * sqrt(weights)
    }

    # Demean slope_vars by FEs
    if (!is.null(mats$slope_vars)) {
      mats$slope_vars <- fixest::demean(
        as.matrix(mats$slope_vars),
        f,
        weights = weights
      )
      if (!is.null(weights)) mats$slope_vars <- mats$slope_vars * sqrt(weights)
    }
  } else if (is.null(model$onlyFixef)) {
    mats$rhs <- fixest::sparse_model_matrix(
      model,
      type = c("rhs"),
      collin.rm = TRUE, na.rm = TRUE
    )

    if (!is.null(weights)) mats$rhs <- mats$rhs * sqrt(weights)
  }

  # Caclulate P_ii's
  if (!is.null(mats$fixef)) {
    # https://stackoverflow.com/questions/39533785/how-to-compute-diagx-solvea-tx-efficiently-without-taking-matrix-i
    U <- Matrix::t(Matrix::chol(Matrix::crossprod(mats$fixef)))

    if (is.null(cl)) {
      # Solve U' X
      Z <- Matrix::solve(U, Matrix::t(mats$fixef))
      P_list$fixef <- Matrix::colSums(Z^2) * T
    } else {
      # Cluster-by-cluster, solve U' X_c
      P_fixef_list <- lapply(cl_split, function(cl_idx) {
        T_cl <- T[cl_idx]
        Z_cl <- Matrix::solve(
          U, Matrix::t(mats$fixef)[, cl_idx, drop = FALSE]
        )
        P_T_cl <- Matrix::crossprod(Z_cl) %*% T_cl

        # Extract those elements into the larger P matrix col and row indices
        Matrix::sparseVector(
          x = P_T_cl[, 1],
          i = cl_idx,
          length = n
        )
      })
      P_list$fixef <- Reduce("+", P_fixef_list)
    }
  }

  # Note: This might fail if looking at too large of a matrix
  if (!is.null(mats$slope_vars)) {
    if (is.null(cl)) {
      Z <- Matrix::solve(
        Matrix::crossprod(mats$slope_vars),
        Matrix::t(mats$slope_vars)
      )
      P_list$slope_vars <- Matrix::rowSums(mats$slope_vars * Matrix::t(Z)) * T
    } else {
      tSlSl <- Matrix::crossprod(mats$slope_vars)
      P_slope_vars_list <- lapply(cl_split, function(cl_idx) {
        T_cl <- T[cl_idx]
        tSl_cl <- Matrix::t(mats$slope_vars)[, cl_idx, drop = FALSE]
        P_T_cl <- Matrix::t(tSl_cl) %*% Matrix::solve(tSlSl, tSl_cl %*% T_cl)

        # Extract those elements into the larger P matrix col and row indices
        Matrix::sparseVector(
          x = P_T_cl[, 1],
          i = cl_idx,
          length = n
        )
      })
      P_list$slope_vars <- Reduce("+", P_slope_vars_list)

      # head(diag(mats$slope_vars %*% Matrix::solve(Matrix::crossprod(mats$slope_vars), Matrix::t(mats$slope_vars))))
      # head(diag(P_list$slope_vars))
    }
  }

  if (!is.null(mats$rhs)) {
    tXXinv <- model$cov.iid
    if (method == "feols") tXXinv <- tXXinv / model$sigma2

    if (is.null(cl)) {
      P_list$rhs <- Matrix::diag(mats$rhs %*% tXXinv %*% Matrix::t(mats$rhs)) * T
    } else {
      P_rhs_list <- lapply(cl_split, function(cl_idx) {
        T_cl <- T[cl_idx]
        tX_cl <- Matrix::t(mats$rhs)[, cl_idx, drop = FALSE]
        P_T_cl <- Matrix::t(tX_cl) %*% tXXinv %*% tX_cl %*% T_cl

        Matrix::sparseVector(
          x = P_T_cl[, 1],
          i = cl_idx,
          length = n
        )
      })
      P_list$rhs <- Reduce("+", P_rhs_list)

      # head(diag(mats$rhs %*% (tXXinv %*% Matrix::t(mats$rhs))))
      # head(diag(P_list$rhs))
    }
  }

  P <- Reduce("+", P_list)
  if (!inherits(P, "dsparseVector")) {
    P <- Matrix::sparseVector(x = P, i = 1:n, length = n)
  }
  return(P)
}

#' Computes (I - D_X)^{-1} T for diagonal block projection matrix
#' @param model A `fixest` model object.
#' @param T Numeric vector
#' @param cl A vector of cluster membership with the same length as the number
#' of observations used in `model`.
solve_ImDX_T <- function(model, T, cl = NULL) {
  
  if (!is.null(cl)) cl_split <- split(1:model$nobs, cl)
  n <- model$nobs
  method <- "feols"

  # Outline: Blockwise Formula for Projection Matrix
  # https://en.wikipedia.org/wiki/Projection_matrix#Blockwise_formula
  # Part 1: fixed effects
  # Part 2: slope vars (demeaned by fixed effects)
  # Part 3: rhs vars (demeaned by fixed effects and slope vars)
  # Sum the three parts to get the projection matrix
  P_list <- list()
  mats <- list()
  mats$fixef <- fixest::sparse_model_matrix(
    model,
    type = c("fixef"),
    collin.rm = TRUE, na.rm = TRUE
  )

  # Check for slope.vars and move to seperate matrix
  if (!is.null(model$slope_flag)) {
    slope_var_cols <- which(
      grepl("\\[\\[", colnames(mats$fixef))
    )
    mats$slope_vars <- mats$fixef[, slope_var_cols]
    mats$fixef <- mats$fixef[, -slope_var_cols]
  }

  # Get weights
  weights <- model$weights

  if (!is.null(model$fixef_vars)) {
    f <- stats::model.matrix(model, type = "fixef")
    f <- subset(f, select = model$fixef_vars)

    # Matrix of fixed effects times weights
    if (!is.null(weights)) mats$fixef <- mats$fixef * sqrt(weights)

    # Demean X by FEs and slope_vars
    if (is.null(model$onlyFixef) & method == "feols") {
      mats$rhs <- fixest::demean(model)[, -1, drop = FALSE]
      if (!is.null(weights)) mats$rhs <- mats$rhs * sqrt(weights)
    }
    if (is.null(model$onlyFixef) & method == "feglm") {
      mats$rhs <- fixest::demean(
        stats::model.matrix(model, "rhs"),
        f,
        weights = weights
      )
      if (!is.null(weights)) mats$rhs <- mats$rhs * sqrt(weights)
    }

    # Demean slope_vars by FEs
    if (!is.null(mats$slope_vars)) {
      mats$slope_vars <- fixest::demean(
        as.matrix(mats$slope_vars),
        f,
        weights = weights
      )
      if (!is.null(weights)) mats$slope_vars <- mats$slope_vars * sqrt(weights)
    }
  } else if (is.null(model$onlyFixef)) {
    mats$rhs <- fixest::sparse_model_matrix(
      model,
      type = c("rhs"),
      collin.rm = TRUE, na.rm = TRUE
    )

    if (!is.null(weights)) mats$rhs <- mats$rhs * sqrt(weights)
  }

  # Caclulate P_ii's
  if (is.null(cl)) {
    if (!is.null(mats$fixef)) {
      U <- Matrix::t(Matrix::chol(Matrix::crossprod(mats$fixef)))
      Z <- Matrix::solve(U, Matrix::t(mats$fixef))
      P_list$fixef <- Matrix::colSums(Z^2)
    }
    if (!is.null(mats$slope_vars)) {
      Z <- Matrix::solve(
        Matrix::crossprod(mats$slope_vars),
        Matrix::t(mats$slope_vars)
      )
      P_list$slope_vars <- Matrix::rowSums(mats$slope_vars * Matrix::t(Z))
    }
    if (!is.null(mats$rhs)) {
      tXXinv <- model$cov.iid / model$sigma2
      P_list$rhs <- Matrix::diag(mats$rhs %*% tXXinv %*% Matrix::t(mats$rhs))
    }
    P = Reduce("+", P_list)
    solve = T / (1 - P)
    return(solve)

  } else {
    if (!is.null(mats$fixef)) {
      U <- Matrix::t(Matrix::chol(Matrix::crossprod(mats$fixef)))
    }
    if (!is.null(mats$slope_vars)) {
      tSlSl <- Matrix::crossprod(mats$slope_vars)
    }
    if (!is.null(mats$rhs)) {
      tXXinv <- model$cov.iid / model$sigma2
    }

    solve_list <- lapply(cl_split, function(cl_idx) {
      T_cl <- T[cl_idx]

      P_list_cl <- list()
      if (!is.null(mats$fixef)) {
        Z_cl <- Matrix::solve(
          U, Matrix::t(mats$fixef)[, cl_idx, drop = FALSE]
        )
        P_list_cl$fixef <- Matrix::crossprod(Z_cl)
      }
      if (!is.null(mats$slope_vars)) {
        tSl_cl <- Matrix::t(mats$slope_vars)[, cl_idx, drop = FALSE]
        P_list_cl$slope_vars <- Matrix::t(tSl_cl) %*% Matrix::solve(tSlSl, tSl_cl)
      }
      if (!is.null(mats$rhs)) {
        tX_cl <- Matrix::t(mats$rhs)[, cl_idx, drop = FALSE]
        P_list_cl$rhs <- Matrix::t(tX_cl) %*% tXXinv %*% tX_cl
      }
      P_cl <- Reduce("+", P_list_cl)
      ImD_cl <- Matrix::Diagonal(nrow(P_cl)) - P_cl
      solve_ImD_cl_T_cl <- Matrix::solve(ImD_cl, T_cl)

      # Extract those elements into the larger P matrix col and row indices
      Matrix::sparseVector(
        x = as.numeric(solve_ImD_cl_T_cl),
        i = cl_idx,
        length = n
      )
    })

    solve <- Reduce("+", solve_list)
    if (!inherits(solve, "dsparseVector")) {
      solve <- Matrix::sparseVector(x = solve, i = 1:n, length = n)
    }
    return(solve)
  }
  
}

#' Computes (I - (D_X1 - D_X2))^{-1} T for diagonal block projection matrix
#' @param model1 A `fixest` model object. This should be the complete model (W Z)
#' @param model2 A `fixest` model object. This should be the simpler model (W)
#' @param T Numeric vector
#' @param cl A vector of cluster membership with the same length as the number
#' of observations used in `model`.
solve_ImDX1pDX2_T <- function(model1, model2, T, cl = NULL) {
  
  if (!is.null(cl)) cl_split <- split(1:model1$nobs, cl)
  n <- model1$nobs
  method <- "feols"

  # Outline: Blockwise Formula for Projection Matrix
  # https://en.wikipedia.org/wiki/Projection_matrix#Blockwise_formula
  # Part 1: fixed effects
  # Part 2: slope vars (demeaned by fixed effects)
  # Part 3: rhs vars (demeaned by fixed effects and slope vars)
  # Sum the three parts to get the projection matrix
  P_list <- list()
  mats1 <- list()
  mats1$fixef <- fixest::sparse_model_matrix(
    model1,
    type = c("fixef"),
    collin.rm = TRUE, na.rm = TRUE
  )
  if (!is.null(model1$slope_flag)) {
    slope_var_cols <- which(
      grepl("\\[\\[", colnames(mats1$fixef))
    )
    mats1$slope_vars <- mats1$fixef[, slope_var_cols]
    mats1$fixef <- mats1$fixef[, -slope_var_cols]
  }
  weights1 <- model1$weights
  if (!is.null(model1$fixef_vars)) {
    f1 <- stats::model.matrix(model1, type = "fixef")
    f1 <- subset(f1, select = model1$fixef_vars)

    # Matrix of fixed effects times weights
    if (!is.null(weights1)) mats1$fixef <- mats1$fixef * sqrt(weights1)

    # Demean X by FEs and slope_vars
    if (is.null(model1$onlyFixef)) {
      mats1$rhs <- fixest::demean(model1)[, -1, drop = FALSE]
      if (!is.null(weights1)) mats1$rhs <- mats1$rhs * sqrt(weights1)
    }

    # Demean slope_vars by FEs
    if (!is.null(mats1$slope_vars)) {
      mats1$slope_vars <- fixest::demean(
        as.matrix(mats1$slope_vars),
        f1,
        weights = weights1
      )
      if (!is.null(weights1)) mats1$slope_vars <- mats1$slope_vars * sqrt(weights1)
    }
  } else if (is.null(model1$onlyFixef)) {
    mats1$rhs <- fixest::sparse_model_matrix(
      model1,
      type = c("rhs"),
      collin.rm = TRUE, na.rm = TRUE
    )

    if (!is.null(weights1)) mats1$rhs <- mats1$rhs * sqrt(weights1)
  }

  mats2 <- list()
  mats2$fixef <- fixest::sparse_model_matrix(
    model2,
    type = c("fixef"),
    collin.rm = TRUE, na.rm = TRUE
  )
  weights2 <- model2$weights
  if (!is.null(model2$slope_flag)) {
    slope_var_cols <- which(
      grepl("\\[\\[", colnames(mats2$fixef))
    )
    mats2$slope_vars <- mats2$fixef[, slope_var_cols]
    mats2$fixef <- mats2$fixef[, -slope_var_cols]
  }
  if (!is.null(model2$fixef_vars)) {
    f2 <- stats::model.matrix(model2, type = "fixef")
    f2 <- subset(f2, select = model2$fixef_vars)

    # Matrix of fixed effects times weights
    if (!is.null(weights2)) mats2$fixef <- mats2$fixef * sqrt(weights2)

    # Demean X by FEs and slope_vars
    if (is.null(model2$onlyFixef)) {
      mats2$rhs <- fixest::demean(model2)[, -1, drop = FALSE]
      if (!is.null(weights2)) mats2$rhs <- mats2$rhs * sqrt(weights2)
    }

    # Demean slope_vars by FEs
    if (!is.null(mats2$slope_vars)) {
      mats2$slope_vars <- fixest::demean(
        as.matrix(mats2$slope_vars),
        f2,
        weights = weights2
      )
      if (!is.null(weights2)) mats2$slope_vars <- mats2$slope_vars * sqrt(weights2)
    }
  } else if (is.null(model2$onlyFixef)) {
    mats2$rhs <- fixest::sparse_model_matrix(
      model2,
      type = c("rhs"),
      collin.rm = TRUE, na.rm = TRUE
    )

    if (!is.null(weights2)) mats2$rhs <- mats2$rhs * sqrt(weights2)
  }

  # Caclulate P_ii's
  if (is.null(cl)) {
    if (!is.null(mats1$fixef)) {
      U <- Matrix::t(Matrix::chol(Matrix::crossprod(mats1$fixef)))
      Z <- Matrix::solve(U, Matrix::t(mats1$fixef))
      P_list$fixef1 <- Matrix::colSums(Z^2)
    }
    if (!is.null(mats1$slope_vars)) {
      Z <- Matrix::solve(
        Matrix::crossprod(mats1$slope_vars),
        Matrix::t(mats1$slope_vars)
      )
      P_list$slope_vars1 <- Matrix::rowSums(mats1$slope_vars * Matrix::t(Z))
    }
    if (!is.null(mats1$rhs)) {
      tXXinv <- model1$cov.iid / model1$sigma2

      P_list$rhs1 <- Matrix::diag(mats1$rhs %*% tXXinv %*% Matrix::t(mats1$rhs))
    }
    if (!is.null(mats2$fixef)) {
      U <- Matrix::t(Matrix::chol(Matrix::crossprod(mats2$fixef)))
      Z <- Matrix::solve(U, Matrix::t(mats2$fixef))
      P_list$fixef2 <- - Matrix::colSums(Z^2)
    }
    if (!is.null(mats2$slope_vars)) {
      Z <- Matrix::solve(
        Matrix::crossprod(mats2$slope_vars),
        Matrix::t(mats2$slope_vars)
      )
      P_list$slope_vars2 <- - Matrix::rowSums(mats2$slope_vars * Matrix::t(Z))
    }
    if (!is.null(mats2$rhs)) {
      tXXinv <- model2$cov.iid / model2$sigma2

      P_list$rhs2 <- - Matrix::diag(mats2$rhs %*% tXXinv %*% Matrix::t(mats2$rhs))
    }
    P = Reduce("+", P_list)
    solve = T / (1 - P)
    return(solve)

  } else {
    if (!is.null(mats1$fixef)) {
      U1 <- Matrix::t(Matrix::chol(Matrix::crossprod(mats1$fixef)))
    }
    if (!is.null(mats1$slope_vars)) {
      tSlSl1 <- Matrix::crossprod(mats1$slope_vars)
    }
    if (!is.null(mats1$rhs)) {
      tXXinv1 <- model1$cov.iid / model1$sigma2
    }
    if (!is.null(mats2$fixef)) {
      U2 <- Matrix::t(Matrix::chol(Matrix::crossprod(mats2$fixef)))
    }
    if (!is.null(mats2$slope_vars)) {
      tSlSl2 <- Matrix::crossprod(mats2$slope_vars)
    }
    if (!is.null(mats2$rhs)) {
      tXXinv2 <- model2$cov.iid / model2$sigma2
    }

    solve_list <- lapply(cl_split, function(cl_idx) {
      T_cl <- T[cl_idx]

      P_list_cl <- list()
      if (!is.null(mats1$fixef)) {
        Z_cl <- Matrix::solve(
          U1, Matrix::t(mats1$fixef)[, cl_idx, drop = FALSE]
        )
        P_list_cl$fixef1 <- Matrix::crossprod(Z_cl)
      }
      if (!is.null(mats1$slope_vars)) {
        tSl_cl <- Matrix::t(mats1$slope_vars)[, cl_idx, drop = FALSE]
        P_list_cl$slope_vars1 <- Matrix::t(tSl_cl) %*% Matrix::solve(tSlSl1, tSl_cl)
      }
      if (!is.null(mats1$rhs)) {
        tX_cl <- Matrix::t(mats1$rhs)[, cl_idx, drop = FALSE]
        P_list_cl$rhs1 <- Matrix::t(tX_cl) %*% tXXinv1 %*% tX_cl
      }
      if (!is.null(mats2$fixef)) {
        Z_cl <- Matrix::solve(
          U2, Matrix::t(mats2$fixef)[, cl_idx, drop = FALSE]
        )
        P_list_cl$fixef2 <- -Matrix::crossprod(Z_cl)
      }
      if (!is.null(mats2$slope_vars)) {
        tSl_cl <- Matrix::t(mats2$slope_vars)[, cl_idx, drop = FALSE]
        P_list_cl$slope_vars2 <- -Matrix::t(tSl_cl) %*% Matrix::solve(tSlSl2, tSl_cl)
      }
      if (!is.null(mats2$rhs)) {
        tX_cl <- Matrix::t(mats2$rhs)[, cl_idx, drop = FALSE]
        P_list_cl$rhs2 <- -Matrix::t(tX_cl) %*% tXXinv2 %*% tX_cl
      }
      P_cl <- Reduce("+", P_list_cl)
      ImD_cl <- Matrix::Diagonal(nrow(P_cl)) - P_cl
      solve_ImD_cl_T_cl <- Matrix::solve(ImD_cl, T_cl)

      # Extract those elements into the larger P matrix col and row indices
      Matrix::sparseVector(
        x = as.numeric(solve_ImD_cl_T_cl),
        i = cl_idx,
        length = n
      )
    })

    solve <- Reduce("+", solve_list)
    if (!inherits(solve, "dsparseVector")) {
      solve <- Matrix::sparseVector(x = solve, i = 1:n, length = n)
    }
    return(solve)
  }
  
}
