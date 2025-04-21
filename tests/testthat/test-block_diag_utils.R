library(fixest)
library(Matrix)

test_that("mult_DX_T works with Covariates", {
  model <- feols(mpg ~ hp + wt, mtcars)
  T <- mtcars$drat

  X <- sparse_model_matrix(model, c("rhs", "fixef"))
  P <- (X %*% Matrix::solve(Matrix::crossprod(X), Matrix::t(X)))

  prod_manual <- Matrix::diag(P) * T
  prod <- mult_DX_T(model, T)
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )

  prod <- mult_DX_T(model, T, cl = 1:nrow(mtcars))
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )

  cl <- mtcars$am
  D_cl <- P
  D_cl[mtcars$am == 1, mtcars$am == 0] <- 0
  D_cl[mtcars$am == 0, mtcars$am == 1] <- 0
  prod_manual <- D_cl %*% T
  prod <- mult_DX_T(model, T, cl = mtcars$am)
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )
})

test_that("mult_DX_T works with Covariates & FEs", {
  model <- feols(mpg ~ wt | cyl, mtcars)
  T <- mtcars$drat

  X <- sparse_model_matrix(model, c("rhs", "fixef"))
  P <- (X %*% Matrix::solve(Matrix::crossprod(X), Matrix::t(X)))
  prod_manual <- Matrix::diag(P) * T
  prod <- mult_DX_T(model, T)
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )

  prod <- mult_DX_T(model, T, cl = 1:nrow(mtcars))
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )

  cl <- mtcars$am
  D_cl <- P
  D_cl[mtcars$am == 1, mtcars$am == 0] <- 0
  D_cl[mtcars$am == 0, mtcars$am == 1] <- 0
  prod_manual <- D_cl %*% T
  prod <- mult_DX_T(model, T, cl = mtcars$am)
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )
})

test_that("mult_DX_T works with Covariates & FEs", {
  model <- feols(mpg ~ wt | cyl[hp], mtcars)
  T <- mtcars$drat

  X <- sparse_model_matrix(model, c("rhs", "fixef"))
  P <- (X %*% Matrix::solve(Matrix::crossprod(X), Matrix::t(X)))
  prod_manual <- Matrix::diag(P) * T
  prod <- mult_DX_T(model, T)
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )

  prod <- mult_DX_T(model, T, cl = 1:nrow(mtcars))
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )

  cl <- mtcars$am
  D_cl <- P
  D_cl[mtcars$am == 1, mtcars$am == 0] <- 0
  D_cl[mtcars$am == 0, mtcars$am == 1] <- 0
  prod_manual <- D_cl %*% T
  prod <- mult_DX_T(model, T, cl = mtcars$am)
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )
})

test_that("solve_ImDX_T works with Covariates", {
  model <- feols(mpg ~ hp + wt, mtcars)
  T <- mtcars$drat

  X <- sparse_model_matrix(model, c("rhs", "fixef"))
  P <- (X %*% Matrix::solve(Matrix::crossprod(X), Matrix::t(X)))
  D <- Matrix::Diagonal(nrow(P), Matrix::diag(P))
  ImD <- Matrix::Diagonal(nrow(P)) - D

  prod_manual <- Matrix::solve(ImD, T)
  prod <- solve_ImDX_T(model, T)
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )

  prod <- solve_ImDX_T(model, T, cl = 1:nrow(mtcars))
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )

  cl = mtcars$am
  D_cl = P
  D_cl[mtcars$am == 1, mtcars$am == 0] = 0
  D_cl[mtcars$am == 0, mtcars$am == 1] = 0
  ImD_cl = Matrix::Diagonal(nrow(D_cl)) - D_cl
  prod_manual <- Matrix::solve(ImD_cl, T)
  prod <- solve_ImDX_T(model, T, cl = mtcars$am)
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )
})

test_that("solve_ImDX_T works with Covariates & FEs", {
  model <- feols(mpg ~ hp + wt | cyl, mtcars)
  T <- mtcars$drat

  X <- sparse_model_matrix(model, c("rhs", "fixef"))
  P <- (X %*% Matrix::solve(Matrix::crossprod(X), Matrix::t(X)))
  ImD <- Matrix::Diagonal(nrow(P)) - Matrix::Diagonal(nrow(P), Matrix::diag(P))

  prod_manual <- T / (1 - Matrix::diag(P))
  prod <- solve_ImDX_T(model, T)
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )

  prod <- solve_ImDX_T(model, T, cl = 1:nrow(mtcars))
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )

  cl = mtcars$am
  D_cl = P
  D_cl[mtcars$am == 1, mtcars$am == 0] = 0
  D_cl[mtcars$am == 0, mtcars$am == 1] = 0
  ImD_cl = Matrix::Diagonal(nrow(D_cl)) - D_cl
  prod_manual <- Matrix::solve(ImD_cl, T)
  prod <- solve_ImDX_T(model, T, cl = mtcars$am)
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )
})

test_that("solve_ImDX_T works with Covariates & FEs & Slope FEs", {
  model <- feols(mpg ~ wt | cyl[hp], mtcars)
  T <- mtcars$drat

  X <- sparse_model_matrix(model, c("rhs", "fixef"))
  P <- (X %*% Matrix::solve(Matrix::crossprod(X), Matrix::t(X)))
  ImD <- Matrix::Diagonal(nrow(P)) - Matrix::Diagonal(nrow(P), Matrix::diag(P))

  prod_manual <- T / (1 - Matrix::diag(P))
  prod <- solve_ImDX_T(model, T)
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )

  prod <- solve_ImDX_T(model, T, cl = 1:nrow(mtcars))
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )

  cl = mtcars$am
  D_cl = P
  D_cl[mtcars$am == 1, mtcars$am == 0] = 0
  D_cl[mtcars$am == 0, mtcars$am == 1] = 0
  ImD_cl = Matrix::Diagonal(nrow(D_cl)) - D_cl
  prod_manual <- Matrix::solve(ImD_cl, T)
  prod <- solve_ImDX_T(model, T, cl = mtcars$am)
  expect_equal(
    as.numeric(prod_manual),
    as.numeric(prod)
  )
})

test_that("solve_ImDX1mDX2_T works", {
  model1 <- feols(mpg ~ hp + wt | cyl, mtcars)
  model2 <- feols(mpg ~ 0 + hp + wt, mtcars)
  T <- mtcars$drat

  X1 <- sparse_model_matrix(model1, c("rhs", "fixef"))
  X2 <- sparse_model_matrix(model2, c("rhs", "fixef"))
  P1 <- (X1 %*% Matrix::solve(Matrix::crossprod(X1), Matrix::t(X1)))
  P2 <- (X2 %*% Matrix::solve(Matrix::crossprod(X2), Matrix::t(X2)))
  I <- Matrix::Diagonal(nrow(P2))
  M2 <- I - P2
  MX2_Z <- (M2 %*% X1)[, 3:5]
  P_MX2_Z <- (MX2_Z %*%
    Matrix::solve(Matrix::crossprod(MX2_Z), Matrix::t(MX2_Z)))
  # Block decomposition of projection matrix check
  expect_equal(
    (P1 - P2),
    P_MX2_Z
  )

  D1 <- Matrix::Diagonal(nrow(P1), Matrix::diag(P1))
  D2 <- Matrix::Diagonal(nrow(P2), Matrix::diag(P2))
  ImDX1pDX2_manual <- Matrix::solve(I - (D1 - D2), T)
  ImDX1pDX2 <- solve_ImDX1pDX2_T(model1, model2, T)
  expect_equal(
    as.numeric(ImDX1pDX2_manual),
    as.numeric(ImDX1pDX2)
  )

  ImDX1pDX2 <- solve_ImDX1pDX2_T(model1, model2, T, cl = 1:nrow(mtcars))
  expect_equal(
    as.numeric(ImDX1pDX2_manual),
    as.numeric(ImDX1pDX2)
  )

  cl = mtcars$am
  diff = (P1 - P2)
  diff[mtcars$am == 1, mtcars$am == 0] = 0
  diff[mtcars$am == 0, mtcars$am == 1] = 0
  ImD_cl = Matrix::Diagonal(nrow(diff)) - diff
  ImDX1pDX2_manual <- Matrix::solve(ImD_cl, T)
  ImDX1pDX2 <- solve_ImDX1pDX2_T(model1, model2, T, cl = mtcars$am)
  expect_equal(
    as.numeric(ImDX1pDX2_manual),
    as.numeric(ImDX1pDX2)
  )
})
