data(tspectra)
data(pheno_suger)
test_that("check this funcion is ok or not ", {
  testthat::expect_output(preprocess_compare(x = tspectra[,1:100],y = pheno_suger[,10],cv = 10,ncomp = c(2,3,4,5)))
})
