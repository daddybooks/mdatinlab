test_that("check this function is OK or not", {
  testthat::expect_true(
    all(check_ifout(HPLC_dat))
  )
})
