#' Block Diagonal Hat Matrix
#' 
#' Calculates the block diagonal hat matrix for a fixest model where the blocks 
#' correspond to clusters (cl_i == cl_j) otherwise the (i,j) element is 
#' zeroed out. Note the data doesn't have to be sorted for this.
#' 
#' @param model A `fixest` model object.
#' @param cl A vector of cluster membership with the same length as the number 
#' of observations used in `model`.
#' 
#' @return A `nobs` by `nobs` sparse matrix of the block diagonal hat matrix.
#' 
#' @export
block_diag_hatvalues <- function(model, cl = NULL) {
  
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
    U <- Matrix::chol(Matrix::crossprod(mats$fixef))

    if (is.null(cl)) {
      # Solve U' X
      Z <- Matrix::solve(Matrix::t(U), Matrix::t(mats$fixef))
      P_list$fixef <- Matrix::Diagonal(n, Matrix::colSums(Z^2))
    } else {
      # Cluster-by-cluster, solve U' X_c
      P_fixef_list <- lapply(cl_split, function(cl_idx) {
        Z_cl <- Matrix::solve(
          Matrix::t(U), Matrix::t(mats$fixef)[, cl_idx, drop = FALSE]
        )
        P_cl <- Matrix::crossprod(Z_cl)

        # Extract those elements into the larger P matrix col and row indices
        trip <- methods::as(P_cl, "TsparseMatrix")
        Matrix::sparseMatrix(
          i = cl_idx[trip@i + 1],
          j = cl_idx[trip@j + 1],
          x = trip@x,
          symmetric = TRUE,
          dims = c(n, n)
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
      # P_ii = diag(X %*% Z) = rowSums(X * Z)
      P_list$slope_vars <- Matrix::Diagonal(
        n, Matrix::rowSums(mats$slope_vars * Matrix::t(Z))
      )
    } else {
      tSlSl <- Matrix::crossprod(mats$slope_vars)
      P_slope_vars_list <- lapply(cl_split, function(cl_idx) {
        tSl_cl <- Matrix::t(mats$slope_vars)[, cl_idx, drop = FALSE]
        P_cl <- Matrix::t(tSl_cl) %*% Matrix::solve(tSlSl, tSl_cl)

        # Extract those elements into the larger P matrix col and row indices
        trip <- methods::as(P_cl, "TsparseMatrix")
        Matrix::sparseMatrix(
          i = cl_idx[trip@i + 1],
          j = cl_idx[trip@j + 1],
          x = trip@x,
          symmetric = TRUE,
          dims = c(n, n)
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
      P_list$rhs <- Matrix::Diagonal(
        n, Matrix::rowSums(mats$rhs * Matrix::t(tXXinv %*% Matrix::t(mats$rhs)))
      )
    } else {
      P_rhs_list <- lapply(cl_split, function(cl_idx) {
        tX_cl <- Matrix::t(mats$rhs)[, cl_idx, drop = FALSE]
        P_cl <- Matrix::t(tX_cl) %*% tXXinv %*% tX_cl

        # Extract those elements into the larger P matrix col and row indices
        trip <- methods::as(P_cl, "TsparseMatrix")
        Matrix::sparseMatrix(
          i = cl_idx[trip@i + 1],
          j = cl_idx[trip@j + 1],
          x = trip@x,
          # This is an optimization for when there are no FEs
          symmetric = !is.null(mats$fixef),
          dims = c(n, n)
        )
      })
      P_list$rhs <- Reduce("+", P_rhs_list)

      # head(diag(mats$rhs %*% (tXXinv %*% Matrix::t(mats$rhs))))
      # head(diag(P_list$rhs))
    }
  }

  P <- Reduce("+", P_list)
  # For efficient operations, mark this (correctly) as a sparse matrix
  if (!is.null(cl)) P <- methods::as(P, "symmetricMatrix")
  
  return(P)
}
