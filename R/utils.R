#' Split formula into terms
#' 
#' @param formula Full formula following `fixest` syntax:
#' `y ~ W | W_FE | T ~ Z | Z_FE`.
#' @param parts_as_formula Logical. If `TRUE`, then each part will be a 
#' right-hand side formula. Default is `FALSE`
#' 
#' @return List of expressions/formula for each part of the formula. It will be of type `symbol`/`language` unless `parts_as_formula = TRUE`. Can be used with 
#' `fixest::xpd` and the dot bracket syntax to create formula. Any missing
#' elements will be given a value of `NULL`. The list contains the following:
#' \item{y_fml}{The LHS}
#' \item{W_lin}{The linear part of the exogenous variables}
#' \item{W_FE}{The fixed effects part of the exogenous variables}
#' \item{T_fml}{The endogenous variable}
#' \item{Z_lin}{The linear part of the instruments}
#' \item{Z_FE}{The fixed effects part of the instruments}
#' 
#' @export
get_fml_parts <- function(formula, parts_as_formula = FALSE) {
  has_lhs <- !is_rhs_only(formula)
  fml_split_tilde <- fml_breaker(formula, "~")

  res <- list(
    y_fml = NULL, W_lin = NULL, W_FE = NULL, T_fml = NULL, Z_lin = NULL, Z_FE = NULL
  )

  # LHS
  if (has_lhs) {
    res$y_fml <- fml_split_tilde[[length(fml_split_tilde)]]
    # Drop y and treat everything as RHS only formula
    fml_split_tilde <- fml_split_tilde[-length(fml_split_tilde)]
  }

  # OLS
  if (length(fml_split_tilde) == 1) {
    fml_split <- fml_breaker(fml_split_tilde[[1]], "|")
    if (length(fml_split) == 2) {
      res$W_lin <- fml_split[[2]]
      res$W_FE <- fml_split[[1]]
    } else {
      res$W_lin <- fml_split[[1]]
    }
  }

  # IV
  if (length(fml_split_tilde) == 2) {
    fml_Z_split <- fml_breaker(fml_split_tilde[[1]], "|")
    fml_W_T_split <- fml_breaker(fml_split_tilde[[2]], "|")

    if (length(fml_Z_split) == 2) {
      res$Z_FE <- fml_Z_split[[1]]
      res$Z_lin <- fml_Z_split[[2]]
    } else {
      res$Z_lin <- fml_Z_split[[1]]
    }

    # This works b/c we know there is an IV
    if (length(fml_W_T_split) == 3) {
      res$T_fml <- fml_W_T_split[[1]]
      res$W_FE <- fml_W_T_split[[2]]
      res$W_lin <- fml_W_T_split[[3]]
    } else {
      res$T_fml <- fml_W_T_split[[1]]
      res$W_lin <- fml_W_T_split[[2]]
    }
  }

  if (parts_as_formula) res = lapply(res, rhs_fml_maker)

  return(res)
}

#' Break apart formula (from right to left) based on a symbole (`~` or `|`)
#'
#' @param fml Formula following `fixest` syntax.
#' @param op String. Either `~` or `|`
#' @return list of `symbol` or `language` from right to left that are split at each occurence of `op`.
fml_breaker <- function(fml, op) {
  res <- list()
  k <- 1
  while (is_operator(fml, op) & length(fml) == 3) {
    res[[k]] <- fml[[3]]
    k <- k + 1
    fml <- fml[[2]]
  }

  if (length(fml) == 2) {
    res[[k]] <- fml[[2]]
  } else {
    res[[k]] <- fml
  }
  res
}

# Makes `symbol`/`language` into rhs formula
rhs_fml_maker <- function(rhs) {
  res <- ~.
  res[[2]] <- rhs
  return(res)
}

# Checks if formula is only rhs
is_rhs_only <- function(fml) {
  # e.g. fml = ~ x
  if (length(fml) == 2) {
    return(TRUE)
  }
  # e.g. fml = ~ x | t ~ z
  if (length(fml[[2]]) == 2) {
    return(TRUE)
  }
  return(FALSE)
}

# From fixest
is_operator <- function(x, op) {
  if (length(x) <= 1) {
    return(FALSE)
  } else {
    return(x[[1]] == op)
  }
}

# helper to check arguments from jive
check_args <- function(
    data, formula,
    cluster = NULL, ssc = FALSE) {
  # Check arguments ------------------------------------------------------------
  dreamerr::check_arg(data, "data.frame")
  dreamerr::check_arg(formula, "ts formula var(data)", .data = data)
  if (!is.null(cluster)) {
    dreamerr::check_arg(cluster, "character var(data) | os formula var(data) right(1, 1)", .data = data)
  }
  dreamerr::check_arg(ssc, "logical scalar")
}

