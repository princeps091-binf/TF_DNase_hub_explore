library(renv)

renv::init()

renv::install("tidyverse")
renv::install("furrr")
renv::install("vroom")
renv::install("RcppRoll")
renv::install("Rcpp")

renv::install("bioc::rtracklayer")
renv::install("bioc::GenomicRanges")

renv::snapshot()
