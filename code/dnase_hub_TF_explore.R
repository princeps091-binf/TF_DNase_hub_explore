library(rtracklayer)
library(GenomicRanges)
library(tidyverse)
library(furrr)
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
Build_GRange_fn<-function(chromo,res,bins,res_num){
  inter_cl_Grange<-   GRanges(seqnames=chromo,
                              ranges = IRanges(start=bins,
                                               end=bins + res_num[res]-1
                              ))
  inter_cl_Grange<-GenomicRanges::reduce(inter_cl_Grange)
  return(inter_cl_Grange)
  
}
Build_coord_fn<-function(top_compound_hub_5kb_tbl,spec_res_file){
  coord_tbl<-do.call(bind_rows,map(unique(top_compound_hub_5kb_tbl$chr),function(chromo){
    message(chromo)
    base::load(paste0(spec_res_file,chromo,"_spec_res.Rda"))
    tmp_tbl<-top_compound_hub_5kb_tbl %>% 
      filter(chr==chromo) %>% 
      mutate(bins=chr_spec_res$cl_member[node]) %>% 
      mutate(bins=map(bins,as.numeric)) 
    
  }))
  return(coord_tbl)
}
produce_rle_from_bigwig_fn<-function(grange_set,bigWig_file){
  sel <- rtracklayer::BigWigSelection(tmp_hub)
  coverage_ranges <- rtracklayer::import.bw(bigWig_file, selection = sel)
  coverage_rle <- GenomicRanges::coverage(coverage_ranges, weight = GenomicRanges::score(coverage_ranges))
  return(coverage_rle)
}
#-----------------------------------------

candidate_hub_file<-"~/Documents/multires_bhicect/Bootstrapp_fn/data/DAGGER_tbl/trans_res/GM12878_union_top_trans_res_dagger_tbl.Rda"
spec_res_file<-"~/Documents/multires_bhicect/data/GM12878/spec_res/"
DNase_bigwig<-"~/Documents/multires_bhicect/data/epi_data/GM12878/DNase/ENCFF901GZH_DNase_GM12878.bigWig"
H3K27ac_bigwig<-"~/Documents/multires_bhicect/data/epi_data/GM12878/H3K27ac/ENCFF180LKW_GM12878_FC.bigWig"
# Hotspot have a more expected size distribution
Dnase_hotspot_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/DNase/ENCSR000EMT_DNase_GM12878_hotspots.bed"
PAX5_bigWig_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/PAX5/ENCFF267GFQ_PAX5_GM12878.bigWig"
PAX5_bed_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/PAX5/ENCFF199FYD_PAX5_GM12878.bed.gz"

#-----------------------------------------
DNase_GRange<-import.bed(Dnase_hotspot_file)
PAX5_GRange<-import.bedGraph(PAX5_bed_file)
bwf_manual <-BigWigFile(PAX5_bigWig_file)
bwf_dnase <-BigWigFile(DNase_bigwig)
bwf_H3K27ac <-BigWigFile(H3K27ac_bigwig)

#Read coverage over a region from a bigWig file
sel <- rtracklayer::BigWigSelection(DNase_GRange)
coverage_ranges <- rtracklayer::import.bw(PAX5_bigWig_file, selection = sel)
coverage_rle <- GenomicRanges::coverage(coverage_ranges, weight = GenomicRanges::score(coverage_ranges))


hub_tbl<-data_tbl_load_fn(candidate_hub_file) %>% 
  mutate(res=str_split_fixed(node,"_",2)[,1])

plan(multisession,workers=3)
hub_tbl<-Build_coord_fn(hub_tbl,spec_res_file) %>% 
  mutate(GRange=future_pmap(list(chr,res,bins),function(chromo,res,bins){
    Build_GRange_fn(chromo,res,bins,res_num)
  }))
plan(sequential)
hub_tbl %>% 
  filter(res=="5kb" & chr == "chr19")
tmp_hub<-(hub_tbl %>% 
  filter(res=="5kb" & chr == "chr19" & node == "5kb_20_190_7750000_7845000") %>% 
  dplyr::select(GRange) %>% 
  unlist)[[1]]
chr_DNase_GRange<-DNase_GRange[seqnames(DNase_GRange)=="chr19"]
chr_TF_GRange<-PAX5_GRange[seqnames(PAX5_GRange)=="chr19"]

tf_rle<-produce_rle_from_bigwig_fn(tmp_hub,PAX5_bigWig_file)
dnase_rle<-produce_rle_from_bigwig_fn(tmp_hub,DNase_bigwig)
H3K27ac_rle<-produce_rle_from_bigwig_fn(tmp_hub,H3K27ac_bigwig)

dnase_inter_tbl<-as_tibble(chr_DNase_GRange[unique(queryHits(findOverlaps(chr_DNase_GRange,tmp_hub[1])))],set="dnase")
TF_inter_tbl<-as_tibble(chr_TF_GRange[unique(queryHits(findOverlaps(chr_TF_GRange,tmp_hub[1])))],set="TF")

tibble(pos=start(tmp_hub)[1]:end(tmp_hub)[1],value=as.numeric(dnase_rle$chr1[start(tmp_hub)[1]:end(tmp_hub)[1]]),set="dnase") %>% 
  bind_rows(.,tibble(pos=start(tmp_hub)[1]:end(tmp_hub)[1],value=as.numeric(tf_rle$chr1[start(tmp_hub)[1]:end(tmp_hub)[1]]),set="TF")) %>% 
  ggplot(.,aes(pos,sqrt(value)))+
  geom_line(size=0.05)+
  geom_segment(data = dnase_inter_tbl,aes(x=start,y=9,xend=end,yend=9),size=5)+
  geom_segment(data = TF_inter_tbl,aes(x=start,y=9.5,xend=end,yend=9.5),size=5)+
  geom_hline(yintercept = 1)+
  facet_grid(set~.)
tibble(pos=start(tmp_hub)[1]:end(tmp_hub)[1],value=as.numeric(dnase_rle$chr1[start(tmp_hub)[1]:end(tmp_hub)[1]]),set="dnase") %>% 
  ggplot(.,aes(x=1,log10(value+1)))+
  geom_boxplot()
tibble(pos=start(tmp_hub)[1]:end(tmp_hub)[1],value=as.numeric(tf_rle$chr1[start(tmp_hub)[1]:end(tmp_hub)[1]]),set="TF") %>% 
  ggplot(.,aes(pos,sqrt(value-1)))+
  geom_line(size=0.05)+
  geom_segment(data = dnase_inter_tbl,aes(x=start,y=9,xend=end,yend=9),size=5)+
  geom_segment(data = TF_inter_tbl,aes(x=start,y=9.5,xend=end,yend=9.5),size=5)+
  geom_hline(yintercept = 0)


tibble(pos=start(tmp_hub)[1]:end(tmp_hub)[1],value=as.numeric(H3K27ac_rle$chr19[start(tmp_hub)[1]:end(tmp_hub)[1]]),set="H3K27ac") %>% 
  bind_rows(.,tibble(pos=start(tmp_hub)[1]:end(tmp_hub)[1],value=as.numeric(tf_rle$chr19[start(tmp_hub)[1]:end(tmp_hub)[1]]),set="TF")) %>% 
  ggplot(.,aes(pos,sqrt(value)))+
  geom_line(size=0.05)+
  geom_segment(data = TF_inter_tbl,aes(x=start,y=6.5,xend=end,yend=6.5),size=10)+
  geom_hline(yintercept = 1)+
  facet_grid(set~.)
