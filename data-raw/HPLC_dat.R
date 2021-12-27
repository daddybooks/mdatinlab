## code to prepare `HPLC_dat` dataset goes here
dat<- read.csv("C:/Users/Administrator/Desktop/bioimformatics/2019南瓜GWAS分析数据/tmp12suger.csv",header = F)
HPLC_dat<- convert2frame(dat,dis = 6,npeak = 3)
usethis::use_data(HPLC_dat, overwrite = TRUE)
