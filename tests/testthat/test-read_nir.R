path <- "inst/test-data/"
test_that("nir reading function", {
 expect_equal( read_nir(path) %>% nrow(), 10)
})
