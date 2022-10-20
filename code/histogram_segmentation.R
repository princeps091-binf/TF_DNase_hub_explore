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

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# 2006 Author Hideaki Shimazaki
# Shimazaki and Shinomoto, Neural Comput 19 1503-1527, 2007
# Department of Physics, Kyoto University
# shimazaki at ton.scphys.kyoto-u.ac.jp
sshist <- function(x){
  
  N <- 2: length(unique(x))
  C <- numeric(length(N))
  D<-sapply(N,function(x)diff_range_x/x)
  diff_range_x<-diff(range(x))
  # MAKE parallel
  for (i in 1:length(N)) {
    message(i/length(N))

    edges <- seq(min(x),max(x),length=N[i])
    hp <- hist(x, breaks = edges, plot=FALSE )
    ki <- hp$counts
    
    k <- mean(ki)
    v <- sum((ki-k)^2)/N[i]
    
    C[i] <- (2*k-v)/D[i]^2	#Cost Function
  }
  
  idx <- which.min(C)
  optD <- D[idx]
  
  edges <- seq(min(x),max(x),length=N[idx])
  h = hist(x, breaks = edges )
  rug(x)
  
  return(h)
}

#-----------------------------------------
tf_bigWig_file<-"~/Documents/multires_bhicect/data/epi_data/GM12878/PAX5/ENCFF267GFQ_PAX5_FC_GM12878.bigWig"
dat_file<-"~/Documents/multires_bhicect/data/GM12878/"

bwf_manual <-BigWigFile(tf_bigWig_file)

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
    tf_rle<-produce_rle_from_bigwig_fn(GRange,tf_bigWig_file)[[chromo]]
    bin_bw_dat_tbl<-tibble(pos=start(GRange):end(GRange),value=as.numeric(tf_rle[start(GRange):end(GRange)]))
    sign_pos<-bin_bw_dat_tbl$pos[which(bin_bw_dat_tbl$value>0)]
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
tf_chr_rle<-produce_rle_from_bigwig_fn(odd_Grange,tf_bigWig_file)
set_val_rle<-runValue(tf_chr_rle$chr20)
set_width_rle<-width(tf_chr_rle$chr20)
set_val_tbl<-tibble(
  rep=set_width_rle[which(set_val_rle!=0)],
  value=set_val_rle[which(set_val_rle!=0)]
)

bp_value_tbl<-tibble(value=rep(set_val_tbl$value,times=set_val_tbl$rep))

bw <- 2 * IQR(log10(bp_value_tbl$value)) / length(bp_value_tbl$value)^(1/3)

sshist(bp_value_tbl$value) 
