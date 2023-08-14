test_that("get_fml_parts works with all types of formula", {
  expect_equal(
    get_fml_parts(~ x), 
    list(
      W_lin = (~x)[[2]]
    )
  )
  expect_equal(
    get_fml_parts(~ x | fe),
    list(
      W_lin = (~x)[[2]],
      W_FE = (~ fe)[[2]]
    )
  )
  expect_equal(
    get_fml_parts(~ x | t ~ z),
    list(
      Z_lin = (~ z)[[2]], 
      T_fml = (~ t)[[2]],
      W_lin = (~x)[[2]]
    )
  )
  expect_equal(
    get_fml_parts(~ x | fe | t ~ z),
    list(
      Z_lin = (~ z)[[2]], 
      T_fml = (~ t)[[2]],
      W_FE = (~ fe)[[2]],
      W_lin = (~x)[[2]]
    )
  )
  expect_equal(
    get_fml_parts(~ x | fe | t ~ z | z_fe),
    list(
      Z_FE = (~ z_fe)[[2]], 
      Z_lin = (~ z)[[2]], 
      T_fml = (~ t)[[2]],
      W_FE = (~ fe)[[2]],
      W_lin = (~x)[[2]]
    )
  )
  expect_equal(
    get_fml_parts(y ~ x),
    list(
      y_fml = (~ y)[[2]],
      W_lin = (~x)[[2]]
    )
  )
  expect_equal(
    get_fml_parts(y ~ x | fe),
    list(
      y_fml = (~ y)[[2]],
      W_lin = (~x)[[2]],
      W_FE = (~ fe)[[2]]
    )
  )
  expect_equal(
    get_fml_parts(y ~ x | t ~ z),
    list(
      y_fml = (~ y)[[2]],
      Z_lin = (~ z)[[2]], 
      T_fml = (~ t)[[2]],
      W_lin = (~x)[[2]]
    )
  )
  expect_equal(
    get_fml_parts(y ~ x | t ~ z | z_fe),
    list(
      y_fml = (~ y)[[2]],
      Z_FE = (~ z_fe)[[2]], 
      Z_lin = (~ z)[[2]], 
      T_fml = (~ t)[[2]],
      W_lin = (~x)[[2]]
    )
  )
  expect_equal(
    get_fml_parts(y ~ x | fe | t ~ z),
    list(
      y_fml = (~ y)[[2]],
      Z_lin = (~ z)[[2]], 
      T_fml = (~ t)[[2]],
      W_FE = (~ fe)[[2]],
      W_lin = (~x)[[2]]
    )
  )
  expect_equal(
    get_fml_parts(y ~ x | fe | t ~ z | z_fe),
    list(
      y_fml = (~ y)[[2]],
      Z_FE = (~ z_fe)[[2]], 
      Z_lin = (~ z)[[2]], 
      T_fml = (~ t)[[2]],
      W_FE = (~ fe)[[2]],
      W_lin = (~x)[[2]]
    )
  )
})
