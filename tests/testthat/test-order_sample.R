data("HPLC_dat")
test_that("multiplication works", {
  expect_equal(order_sample(HPLC_dat,samplenames = HPLC_dat$samplenames) %>% nrow(),99)
})
