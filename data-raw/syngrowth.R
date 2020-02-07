## code to prepare `syngrowth` dataset goes here

syngrowth <- read.csv('data-raw/syngrowth.csv')

usethis::use_data(syngrowth, compress="gzip", overwrite=T)
