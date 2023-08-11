# Break apart `fml_full` into formula: 
# y ~ W | W_FE | T ~ Z | Z_FE  
get_fml_parts <- function(fml_full) {
  fml_split_tilde = fixest:::fml_breaker(fml_full, "~")
  instruments_split = fixest:::fml_breaker(fml_split_tilde[[1]], "|")
  exo_endo_split = fixest:::fml_breaker(fml_split_tilde[[2]], "|")
  y = fml_split_tilde[[3]]
  if (length(instruments_split) == 2) {
    Z_FE = fixest:::fml_maker(instruments_split[[1]])
    Z_lin = fixest:::fml_maker(instruments_split[[2]])
  } else {
    Z_FE = ~ 0
    Z_lin = fixest:::fml_maker(instruments_split[[1]])
  }
  if (length(exo_endo_split) == 3) {
    T_fml = fixest:::fml_maker(exo_endo_split[[1]])
    W_FE = fixest:::fml_maker(exo_endo_split[[2]])
    W_lin = fixest:::fml_maker(exo_endo_split[[3]])
  } else {
    T_fml = fixest:::fml_maker(exo_endo_split[[1]])
    W_FE = ~ 0
    W_lin = fixest:::fml_maker(exo_endo_split[[2]])
  }
  y_fml = fixest:::fml_maker(y)

  return(list(
    y_fml = y_fml,
    W_lin = W_lin,
    W_FE = W_FE,
    T_fml = T_fml,
    Z_lin = Z_lin,
    Z_FE = Z_FE
  ))
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
