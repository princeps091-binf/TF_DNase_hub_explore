library(renv)

renv::init()

renv::install("tidyverse")
renv::install("furrr")

renv::install("bioc::rtracklayer")
renv::install("bioc::GenomicRanges")

renv::snapshot()
