## code to prepare `nishiura` dataset goes here
##

nishiura <- data.table::fread("https://raw.githubusercontent.com/aakhmetz/COVID19SerialInterval/master/data/supplemetary%20table.csv")

usethis::use_data(nishiura, overwrite = TRUE)
