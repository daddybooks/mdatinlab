data("HPLC_dat")
test_that("check function", {
  expect_equal(
    check_nsample(x = dat,dis = 6,last = "G129(11.5)"),"last sample name is indentity, you can go on"
  )
})
