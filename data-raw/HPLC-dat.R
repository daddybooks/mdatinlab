## code to prepare `HPLC_dat` dataset goes here
load(file = "./data/HPLC_dat.rda")
usethis::use_data(HPLC_dat, overwrite = TRUE)

HPLC_dat_frame<- convert2frame(HPLC_dat,dis = 6,npeak = 3)
usethis::use_data(HPLC_dat_frame)
