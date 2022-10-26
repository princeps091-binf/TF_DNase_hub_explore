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
  test<-viewApply(Views(as.list(tf_rle@listData)[[chromo]], start=start(chr_GRange), end=end(chr_GRange)),as.numeric)
  return(tibble(chr=chromo,start=start(chr_GRange),end=end(chr_GRange),peak.pval=test))
}))

peak_tbl %>% 
  dplyr::select(peak.fc) %>% 
  unnest(cols=c(peak.fc)) %>% 
  ggplot(.,aes(peak.fc))+
  geom_density()+
  scale_x_log10()
ggsave(paste0("~/Documents/multires_bhicect/TF_DNase_hub_explore/img/",chromo,"PAX5_peak_signal_intensity_distribtion.png"))

dat_res<-"50kb"
chromo<-"chr22"

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

sign_dna_tbl %>% 
  ggplot(.,aes(x=start,xend=end,y=1,yend=1))+geom_segment(size=20)
ggsave(paste0("~/Documents/multires_bhicect/TF_DNase_hub_explore/img/",chromo,"PAX5_over_ctrl.png"))
sign_dna_tbl%>%
  mutate(d=(end-start) + 1) %>% 
  ggplot(.,aes(d))+
  geom_density()+
  scale_x_log10()+
  theme_classic()
ggsave(paste0("~/Documents/multires_bhicect/TF_DNase_hub_explore/img/",chromo,"PAX5_signal_size_distribtion.png"))
# Notice odd concentration of stretches of that length
sign_dna_tbl%>%
  mutate(d=(end-start) + 1) %>% 
  group_by(d) %>% 
  summarise(n=n()) %>% 
  arrange(desc(n))

odd_set_tbl<-sign_dna_tbl%>%
  mutate(d=(end-start) + 1) %>% 
  filter(d==120)

odd_Grange<-   GRanges(seqnames=sign_dna_tbl$chr,
                            ranges = IRanges(start=sign_dna_tbl$start,
                                             end=sign_dna_tbl$end
                            ))
odd_Grange<-GenomicRanges::reduce(odd_Grange)
tf_chr_rle<-produce_rle_from_bigwig_fn(odd_Grange,tf_bigWig_file)
set_val_rle<-runValue(tf_chr_rle$chr22)
set_width_rle<-width(tf_chr_rle$chr22)
set_val_tbl<-tibble(
  rep=set_width_rle[which(set_val_rle!=0)],
value=set_val_rle[which(set_val_rle!=0)]
)

plot(density(log10(rep(set_val_tbl$value,times=set_val_tbl$rep))))


### for example #####

bp_value_tbl<-tibble(value=rep(set_val_tbl$value,times=set_val_tbl$rep))
hist(log10(bp_value_tbl$value), breaks="FD")

bw <- 2 * IQR(log10(bp_value_tbl$value)) / length(bp_value_tbl$value)^(1/3)
length(bp_value_tbl$value)^(1/3)/bw


bp_value_tbl %>% 
  ggplot(.,aes(log10(value)))+
  geom_histogram(binwidth = bw)
bp_value_tbl%>% 
  ggplot(.,aes(value))+
  geom_density()+
  scale_x_log10()+
  theme_classic()
ggsave(paste0("~/Documents/multires_bhicect/TF_DNase_hub_explore/img/",chromo,"PAX5_signal_intensity_distribtion.png"))

plan(multisession,workers=5)
sign_dna_tbl<-sign_dna_tbl %>% 
  mutate(peak.val=future_pmap_dbl(list(chr,start,end),function(chromo,start,end){
    tmp<-as.list(tf_chr_rle@listData)[[chromo]]
    #tibble(chr=chromo,pos=start:end,value=as.numeric(tmp[start:end]))
    sd(as.numeric(tmp[start:end]))
  }))
plan(sequential)


sign_dna_tbl %>% 
  filter(!(is.na(peak.val))) %>%
  mutate(d=(end-start)+1) %>% 
  filter(d==125) %>% 
  ggplot(.,aes(peak.val))+
  geom_density()
sign_dna_tbl %>% 
  filter(!(is.na(peak.val))) %>%
  mutate(d=(end-start)+1) %>% 
  filter(peak.val==0) %>% 
  ggplot(.,aes(d))+
  geom_density()+
  scale_x_log10()

sign_dna_tbl %>% 
  dplyr::select(peak.val) %>% 
  unnest(cols=c(peak.val)) %>% 
  ggplot(.,aes(peak.val))+
  geom_density()+
#  geom_vline(xintercept = c(1.55,1.40,2.32))+
  scale_x_log10()

chr19_sign_tbl<-do.call(bind_rows,sign_dna_tbl$peak.val)

chr19_sign_tbl %>% 
  ggplot(.,aes(pos,value))+geom_point(size=0.01)
