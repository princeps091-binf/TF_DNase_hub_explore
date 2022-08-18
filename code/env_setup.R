library(renv)

renv::init()

renv::install("tidyverse")
renv::install("bioc::rtracklayer")
renv::install("bioc::GenomicRanges")

renv::snapshot()
