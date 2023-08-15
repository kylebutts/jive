test_that("Covariates", {
  model <- feols(mpg ~ hp + wt, mtcars)
  X <- sparse_model_matrix(model, c("rhs", "fixef"))
  P <- (X %*% solve(crossprod(X), t(X)))
  P_block <- block_diag_hatvalues(model, cl = mtcars$cyl)
  
  expect_equal(
    as.matrix(P_block[mtcars$cyl == 6, mtcars$cyl == 6]),
    as.matrix(P[mtcars$cyl == 6, mtcars$cyl == 6])
  )
})

test_that("Covariates & FEs", {
  model <- feols(mpg ~ hp | cyl, mtcars)
  X <- sparse_model_matrix(model, c("rhs", "fixef"))
  P <- (X %*% solve(crossprod(X), t(X)))
  P_block <- block_diag_hatvalues(model, cl = mtcars$cyl)
  
  expect_equal(
    as.matrix(P_block[mtcars$cyl == 6, mtcars$cyl == 6]),
    as.matrix(P[mtcars$cyl == 6, mtcars$cyl == 6])
  )
})

test_that("Covariates, FEs, & slope vars", {
  model <- feols(mpg ~ wt | cyl[hp], mtcars)
  X <- sparse_model_matrix(model, c("rhs", "fixef"))
  P <- (X %*% solve(crossprod(X), t(X)))
  P_block <- block_diag_hatvalues(model, cl = mtcars$cyl)
  
  expect_equal(
    as.matrix(P_block[mtcars$cyl == 6, mtcars$cyl == 6]),
    as.matrix(P[mtcars$cyl == 6, mtcars$cyl == 6])
  )
})

test_that("cl = NULL gives D matrix", {
  model <- feols(mpg ~ wt | cyl[hp], mtcars)
  X <- sparse_model_matrix(model, c("rhs", "fixef"))
  P <- (X %*% solve(crossprod(X), t(X)))
  P_block <- block_diag_hatvalues(model)
  
  expect_equal(
    Matrix::Diagonal(nrow(P_block), diag(P_block)), 
    P_block
  )
  expect_equal(
    diag(P_block),
    diag(P)
  )
})

