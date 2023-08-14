test_that("fml_breaker works with all types of formula", {
  expect_equal(
    fml_breaker(~ x, "~"), 
    list((~x)[[2]])
  )
  expect_equal(
    fml_breaker(~ x | fe, "~"),
    list((~ x | fe)[[2]])
  )
  expect_equal(
    fml_breaker(~ x | t ~ z, "~"),
    list((~ z)[[2]], (~ x | t)[[2]])
  )
  expect_equal(
    fml_breaker(~ x | fe | t ~ z, "~"),
    list((~ z)[[2]], (~ x | fe | t)[[2]])
  )
  expect_equal(
    fml_breaker(~ x | fe | t ~ z | z_fe, "~"),
    list((~ z | z_fe)[[2]], (~ x | fe | t)[[2]])
  )
  expect_equal(
    fml_breaker(y ~ x, "~"),
    list((~x)[[2]], (~y)[[2]])
  )
  expect_equal(
    fml_breaker(y ~ x | fe, "~"),
    list((~x | fe)[[2]], (~y)[[2]])
  )
  expect_equal(
    fml_breaker(y ~ x | t ~ z, "~"),
    list((~z)[[2]], (~x | t)[[2]], (~y)[[2]])
  )
  expect_equal(
    fml_breaker(y ~ x | t ~ z | z_fe, "~"),
    list((~z | z_fe)[[2]], (~x | t)[[2]], (~y)[[2]])
  )
  expect_equal(
    fml_breaker(y ~ x | fe | t ~ z, "~"),
    list((~z)[[2]], (~x | fe | t)[[2]], (~y)[[2]])
  )
  expect_equal(
    fml_breaker(y ~ x | fe | t ~ z | z_fe, "~"),
    list((~z | z_fe)[[2]], (~x | fe | t)[[2]], (~y)[[2]])
  )
})

