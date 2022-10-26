library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(furrr)
library(vroom)
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
produce_rle_from_bigwig_fn<-function(grange_set,bigWig_file){
  sel <- rtracklayer::BigWigSelection(grange_set)
  coverage_ranges <- rtracklayer::import.bw(bigWig_file, selection = sel)
  coverage_rle <- GenomicRanges::coverage(coverage_ranges, weight = GenomicRanges::score(coverage_ranges))
  return(coverage_rle)
}

#-----------------------------------------
tf_bigWig_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/PAX5/ENCFF267GFQ_PAX5_FC_GM12878.bigWig"
tf_pval_bigWig_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/PAX5/ENCFF955FMV_PAX5_pval_GM12878.bigWig"
tf_bed_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/PAX5/ENCFF199FYD_PAX5_GM12878.bed.gz"
dat_file<-"~/Documents/multires_bhicect/data/GM12878/"

tf_bed_GRange<-import.bedGraph(tf_bed_file)
bwf_manual <-BigWigFile(tf_bigWig_file)

tf_rle<-produce_rle_from_bigwig_fn(tf_bed_GRange,tf_bigWig_file)
pval_rle<-produce_rle_from_bigwig_fn(tf_bed_GRange,tf_pval_bigWig_file)
chr_set<-unique(as.character(seqnames(tf_bed_GRange)))
peak_tbl<-do.call(bind_rows,lapply(chr_set,function(chromo){
  message(chromo)
  chr_GRange<-tf_bed_GRange[seqnames(tf_bed_GRange)==chromo]
  tmp_fc<-viewApply(Views(as.list(tf_rle@listData)[[chromo]], start=start(chr_GRange), end=end(chr_GRange)),as.numeric)
  tmp_pval<-viewApply(Views(as.list(pval_rle@listData)[[chromo]], start=start(chr_GRange), end=end(chr_GRange)),as.numeric)
  
  return(tibble(chr=chromo,start=start(chr_GRange),end=end(chr_GRange),peak.fc=tmp_fc,peak.pval=tmp_pval))
}))

peak_tbl %>% 
  dplyr::select(peak.pval) %>% 
  unnest(cols=c(peak.pval)) %>% 
  ggplot(.,aes(peak.pval))+
  geom_density()+
  scale_x_log10()

peak_tbl %>% 
  mutate(pos=pmap(list(start,end),function(start,end)1:(end-start+1)),
         ID=str_c(chr,start,end,sep="_"),
         w=end-start) %>% 
  filter(w<=310) %>% 
  dplyr::select(ID,pos,peak.fc) %>% 
  unnest(cols=c(pos,peak.fc)) %>% 
  mutate(ID=reorder(ID,peak.fc)) %>% 
  ggplot(.,aes(x=pos,y=ID,fill=log10(peak.fc+1)))+
  geom_tile()+
  scale_fill_viridis_c()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())
