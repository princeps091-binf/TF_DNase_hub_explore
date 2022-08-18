library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
options(scipen = 999999999)
res_set <- c('1Mb','500kb','100kb','50kb','10kb','5kb')
res_num <- c(1e6,5e5,1e5,5e4,1e4,5e3)
names(res_num)<-res_set

#-----------------------------------------

#Utils. Fn
data_tbl_load_fn<-function(file){
  out_tbl<-get(base::load(file))
  tmp_obj<-names(mget(base::load(file)))
  rm(list=tmp_obj)
  rm(tmp_obj)
  return(out_tbl)
}
#-----------------------------------------

candidate_hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/GM12878_union_top_trans_res_dagger_tbl.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"
DNase_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/DNase/ENCSR000EMT_DNase_GM12878.bed"
# Hotspot have a more expected size distribution
Dnase_hotspot_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/DNase/ENCSR000EMT_DNase_GM12878_hotspots.bed"
PAX5_bigWig_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/PAX5/ENCFF267GFQ_PAX5_GM12878.bigWig"
#-----------------------------------------
DNase_GRange<-import.bed(Dnase_hotspot_file)

bwf_manual <-BigWigFile(PAX5_bigWig_file)


#Read coverage over a region from a bigWig file
sel <- rtracklayer::BigWigSelection(DNase_GRange)
coverage_ranges <- rtracklayer::import.bw(PAX5_bigWig_file, selection = sel)
coverage_rle = GenomicRanges::coverage(coverage_ranges, weight = GenomicRanges::score(coverage_ranges))

plot(as.numeric(coverage_rle$chr1[start(DNase_GRange)[300]:end(DNase_GRange)[300]]))

