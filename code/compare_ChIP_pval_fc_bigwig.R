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

tf_fc_bigWig_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/PAX5/ENCFF267GFQ_PAX5_FC_GM12878.bigWig"
tf_pv_bigWig_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/PAX5/ENCFF955FMV_PAX5_pval_GM12878.bigWig"
dat_file<-"~/Documents/multires_bhicect/data/GM12878/"

bwfc_manual <-BigWigFile(tf_fc_bigWig_file)
bwpv_manual <-BigWigFile(tf_pv_bigWig_file)

dat_res<-"50kb"
chromo<-"chr20"

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
  #    slice(1:5) %>% 
  #  filter(bins==17450000)
  mutate(fc_d=future_pmap(list(chr,GRange),function(chromo,GRange){
    tf_rle<-produce_rle_from_bigwig_fn(GRange,tf_fc_bigWig_file)[[chromo]]
    bin_bw_dat_tbl<-tibble(pos=start(GRange):end(GRange),value=as.numeric(tf_rle[start(GRange):end(GRange)]))
    sign_pos<-bin_bw_dat_tbl$pos[which(bin_bw_dat_tbl$value>1)]
    border_pos<-which(diff(sign_pos)>1)
    tmp_d<-tibble(chr=chromo,start=c(sign_pos[1],sign_pos[border_pos+1]),end=c(sign_pos[border_pos],sign_pos[length(sign_pos)])) 
    return(tmp_d)
    
  }))
plan(sequential)

sign_dna_tbl<-do.call(bind_rows,bin_tbl$fc_d) 

odd_Grange<-   GRanges(seqnames=sign_dna_tbl$chr,
                       ranges = IRanges(start=sign_dna_tbl$start,
                                        end=sign_dna_tbl$end
                       ))
odd_Grange<-GenomicRanges::reduce(odd_Grange)

tf_chr_fc_rle<-produce_rle_from_bigwig_fn(odd_Grange,tf_fc_bigWig_file)
tf_chr_pv_rle<-produce_rle_from_bigwig_fn(odd_Grange,tf_pv_bigWig_file)

fc_value<-Views(tf_chr_fc_rle$chr20, start=start(odd_Grange), end=end(odd_Grange))
pv_value<-Views(tf_chr_pv_rle$chr20, start=start(odd_Grange), end=end(odd_Grange))

fc_val_vec<-unlist(lapply(fc_value[1:50],as.numeric))
pv_val_vec<-unlist(lapply(pv_value[1:50],as.numeric))
plot(fc_val_vec,pv_val_vec)