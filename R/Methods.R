#' A print method for `jive_est` objects
#'
#' @param x A `jive_est` object from `jive`/`ujive`/`ijive`
#' @param ... Unused
#'
#' @export
print.jive_est <- function(x, ...) {
  se_label = if (x$clustered) "Clustered SE" else "Robust SE"
  beta <- x$beta
  se <- x$se
  tstat <- beta / se
  pval <- 2 * (1 - stats::pnorm(abs(tstat)))

  # Print coefficients
  cat("Coefficients: \n")  

  cmat <- cbind(as.matrix(beta), as.matrix(se), as.matrix(tstat), as.matrix(pval))
  colnames(cmat) <- c("Estimate", se_label, "Z value", "Pr(>z)")
  stats::printCoefmat(cmat, signif.legend = TRUE)

  stringmagic::cat_magic(
    "{n ? n} observations, {n ? K} instruments, {n ? L} covariates",
    n = x$n, L = x$n_covariates, K = x$n_instruments
  )
  stringmagic::cat_magic(
    "\nFirst-stage F: stat = {%0.3f ? F}",
    F = x$F
  )
  if (length(x$Sargan) > 0) {
    stringmagic::cat_magic(
      "\n       Sargan: stat = {%0.3f ? stat}, p = {%0.3f ? pval}", stat = x$Sargan$statistic, pval = x$Sargan$pvalue
    )
  }
  if (length(x$CD) > 0) {
    stringmagic::cat_magic(
      "\n           CD: stat = {%0.3f ? stat}, p = {%0.3f ? pval}", stat = x$CD$statistic, pval = x$CD$pvalue
    )
  }
}

#' A print method for `manyiv_est` objects
#'
#' @param x A `manyiv_est` object from `manyiv`
#' @param ... Unused
#'
#' @export
print.manyiv_est  <- function(x, ...) {
  if (x$clustered) {
    se = x$se[3, ]
    se_label = "Clustered SE"
  } else {
    se = x$se[2, ]
    se_label = "Robust SE"
  }
  beta = x$beta
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
    n = x$n,
    L = x$L,
    K = x$K,
    F = x$F
  )
  return(NULL)
}

#' A print method for `jive_est` objects
#' @param object A `jive_est` object from `jive`/`ujive`/`ijive`
#' @param ... Unused
#' @export
summary.jive_est <- function(object, ...) {
  print(object)
}

#' A print method for `manyiv_est` objects
#' @param object A `manyiv_est` object from `manyiv`
#' @param ... Unused
#' @export
summary.manyiv_est <- function(object, ...) {
  print(object)
}

# #' Extract coefficient from a `jive_est` object
# #' @inheritParams summary.jive_est
# #' @export
# coef.jive_est <- function(object, ...) {
#   return(object$beta)
# }

# #' Extract standard error from a `jive_est` object
# #' @inheritParams summary.jive_est
# #' @export
# se.jive_est <- function(object, ...) {
#   return(object$se)
# }

