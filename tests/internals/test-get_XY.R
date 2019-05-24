test_that("get_XY works for all combinations of X,Y,data, and formula", {
  expect_identical(get_XY(list( formula = f_Ya, data=mae_data)), get_XY(list(X=X, Y=Ya,  data=mae_data)))
  expect_identical(get_XY(list( formula = f_Ya, data=mae_data)), get_XY(list(X=Xm_Ya, Y=Yam)))
})

