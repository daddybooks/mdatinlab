data("HPLC_dat")
test_that("this checking function is OK",{
  testthat::expect_equal(

    check_ifdup(HPLC_dat),"OK!"
  )
})

