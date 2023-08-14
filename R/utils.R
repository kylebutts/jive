#' Break apart formula (from right to left) based on a symbole (`~` or `|`)
#' 
#' @param fml Formula following `fixest` syntax. 
#' @param op String. Either `~` or `|`
#' @return list of `symbol` or `language` from right to left that are split at each occurence of `op`.
#' @export
fml_breaker <- function(fml, op) {
  res = list()
  k = 1
  while (is_operator(fml, op) & length(fml) == 3) {
    res[[k]] = fml[[3]]
    k = k + 1
    fml = fml[[2]]
  }

  if (length(fml) == 2) {
    res[[k]] = fml[[2]]
  } else {
    res[[k]] = fml
  }
  res
}

# Makes `symbol`/`language` into rhs formula
rhs_fml_maker <- function(rhs) {
  res = ~ .
  res[[2]] = rhs
  return(res)
}

#' Break apart `fml_full` into formula: 
#' y ~ W | W_FE | T ~ Z | Z_FE  
get_fml_parts <- function(fml_full) {
  # fml_full = ~ x
  has_lhs = !is_rhs_only(fml_full)
  fml_split_tilde = fml_breaker(fml_full, "~")

  res = list()
  # LHS
  if (has_lhs) {
    res$y_fml = fml_split_tilde[[length(fml_split_tilde)]]
    # Drop y and treat everything as RHS only formula
    fml_split_tilde = fml_split_tilde[-length(fml_split_tilde)]
  }

  # OLS
  if (length(fml_split_tilde) == 1) {
    fml_split = fml_breaker(fml_split_tilde[[1]], "|")
    if (length(fml_split) == 2) { 
      res$W_lin = fml_split[[2]]
      res$W_FE = fml_split[[1]]
    } else {
      res$W_lin = fml_split[[1]]
    }
  }

  # IV
  if (length(fml_split_tilde) == 2) {
    fml_Z_split = fml_breaker(fml_split_tilde[[1]], "|")
    fml_W_T_split = fml_breaker(fml_split_tilde[[2]], "|")
    
    if (length(fml_Z_split) == 2) { 
      res$Z_FE = fml_Z_split[[1]]
      res$Z_lin = fml_Z_split[[2]]
    } else {
      res$Z_lin = fml_Z_split[[1]]
    }

    # This works b/c we know there is an IV
    if (length(fml_W_T_split) == 3) {
      res$T_fml = fml_W_T_split[[1]]
      res$W_FE = fml_W_T_split[[2]]
      res$W_lin = fml_W_T_split[[3]]
    } else {
      res$T_fml = fml_W_T_split[[1]]
      res$W_lin = fml_W_T_split[[2]]
    }
  }

  return(res)
}

# merge rhs of formula
merge_rhs_formula <- function(fml1, fml2) {
  if (is.null(fml1)) fml1 = ~ 0
  if (is.null(fml2)) fml2 = ~ 0
	return(fixest::xpd(~ .[fml1] + .[fml2]))
}

# helper to check arguments from jive
check_args <- function(
  data, y, exogenous, endogenous, instruments,
  cluster = NULL, ssc = FALSE
) { 

  # I could do y + exogenous as one & endogenous + instruments as another,
  # but I think that might mislead folks into thinking there are two regressions
  # occuring
  # Check arguments ------------------------------------------------------------
  dreamerr::check_arg(data, "data.frame")
  dreamerr::check_arg(y, "character var(data) | os formula var(data) right(1, 1)
  ", .data = data)
  if (!is.null(cluster)) {
    dreamerr::check_arg(cluster, "character var(data) | os formula var(data) right(1, 1)", .data = data)
  }
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
}
