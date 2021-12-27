data("dat")
test_that("multiplication works", {
  expect_equal(convert2frame(dat,dis = 6,npeak = 3) %>% nrow,99)
})
