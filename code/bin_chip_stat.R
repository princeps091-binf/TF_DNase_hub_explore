library(rtracklayer)
library(GenomicRanges)
library(RcppRoll)
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

hic_dat_in<-function(dat_file,cl_res,chromo){
  chr_dat<-vroom(paste0(dat_file,cl_res,"/",chromo,".txt"),delim = "\t",col_names = F,trim_ws = T,escape_double = F)
  return(chr_dat%>%mutate(X3=as.numeric(X3))%>%filter(!(is.na(X3)))%>%filter(X1!=X2))
}

Build_GRange_fn<-function(chromo,res,bins,res_num){
  inter_cl_Grange<-   GRanges(seqnames=chromo,
                              ranges = IRanges(start=bins,
                                               end=bins + res_num[res]-1
                              ))
  inter_cl_Grange<-GenomicRanges::reduce(inter_cl_Grange)
  return(inter_cl_Grange)
  
}
produce_rle_from_bigwig_fn<-function(grange_set,bigWig_file){
  sel <- rtracklayer::BigWigSelection(grange_set)
  coverage_ranges <- rtracklayer::import.bw(bigWig_file, selection = sel)
  coverage_rle <- GenomicRanges::coverage(coverage_ranges, weight = GenomicRanges::score(coverage_ranges))
  return(coverage_rle)
}
#-----------------------------------------------------------------
dat_file<-"~/Documents/multires_bhicect/data/GM12878/"
tf_bigWig_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/PAX5/ENCFF267GFQ_PAX5_FC_GM12878.bigWig"
tf_bed_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/PAX5/ENCFF199FYD_PAX5_GM12878.bed.gz"

tf_bed_GRange<-import.bedGraph(tf_bed_file)
bwf_manual <-BigWigFile(tf_bigWig_file)
dat_res<-"50kb"
chromo<-"chr19"
chr_dat<-hic_dat_in(dat_file,dat_res,chromo)
bin_tbl<-tibble(chr=chromo,res=dat_res,bins=unique(c(chr_dat$X1,chr_dat$X2)))
rm(chr_dat)
plan(multisession,workers=5)
bin_tbl<-bin_tbl %>% 
  mutate(GRange=future_pmap(list(chr,res,bins),function(chromo,res,bins){
    Build_GRange_fn(chromo,res,bins,res_num)
  }))
plan(sequential)
plan(multisession,workers=5)
bin_tbl<-bin_tbl %>%
  #  slice(1:5) %>% 
  mutate(fc_stat=future_pmap_dbl(list(chr,GRange),function(chromo,GRange){
    tf_rle<-produce_rle_from_bigwig_fn(GRange,tf_bigWig_file)[[chromo]]
    value_tbl<-tibble(pos=start(GRange):end(GRange),value=as.numeric(tf_rle[start(GRange):end(GRange)]))
    peak_tbl<-tibble(chr=chromo,start=start(IRanges::intersect(GRange,tf_bed_GRange)),end=end(IRanges::intersect(GRange,tf_bed_GRange)))
    if(nrow(peak_tbl)>0){
      peak_s<-sum(unlist(lapply(1:nrow(peak_tbl),function(i){
        sum(value_tbl$value[which(findInterval(value_tbl$pos,c(peak_tbl$start[i],peak_tbl$end[i]))==1)])
      })))
      
    } else{
      peak_s<-0
    }
    return(peak_s)
    
    
  }))
plan(sequential)

plan(multisession,workers=5)
bin_tbl<-bin_tbl %>%
  #  slice(1:5) %>% 
  mutate(fc_val=future_pmap_dbl(list(chr,GRange),function(chromo,GRange){
    tf_rle<-produce_rle_from_bigwig_fn(GRange,tf_bigWig_file)[[chromo]]
    value_tbl<-tibble(pos=start(GRange):end(GRange),value=as.numeric(tf_rle[start(GRange):end(GRange)]))
    return(value_tbl %>%filter(value > 1) %>% summarise(s=sum(value)) %>% unlist)
    
  }))
plan(sequential)
value_tbl %>% ggplot(.,aes(pos,value))+geom_line()+
  geom_segment(data=peak_tbl,aes(x=start,xend=end,y=15,yend=15),size=5)+
  geom_hline(yintercept = 1)
return(value_tbl %>% summarise(sign.fc=sum(value>1),n=n()))
