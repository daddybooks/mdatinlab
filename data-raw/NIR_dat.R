## code to prepare `NIR_dat` dataset goes here
path <- "inst/test-data/"
NIR_dat<- read_nir(path)
usethis::use_data(NIR_dat, overwrite = TRUE)

