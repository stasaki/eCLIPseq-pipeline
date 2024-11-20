
options(stringsAsFactors = FALSE)
library(tidyverse)
library(data.table)
library(rtracklayer)
library(plyranges)
library("optparse")

option_list = list(
  make_option(c("-i", "--sites_count"), type="character", default=NULL, 
              help="sites count", metavar="character"),
  make_option(c("-s", "--input_stats"), type="character", default=NULL, 
              help="input stats", metavar="character"),
  make_option(c("-p", "--ip_stats"), type="character", default=NULL, 
              help="ip stats", metavar="character"),
  make_option(c("-r", "--out_rds"), type="character", default=NULL, 
              help="fisher output", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


ds = data.table::fread(opt$sites_count )
input_total =  read_lines(opt$input_stats )%>%.[2]%>%gsub(",.+","",.)%>%as.numeric()
ip_total =  read_lines(opt$ip_stats)%>%.[2]%>%gsub(",.+","",.)%>%as.numeric()

colnames(ds)[9]="INPUT"
colnames(ds)[8]="IP"

ds%>%
  mutate(dummy=paste0(V1,":",V2))%>%
  filter(!duplicated(dummy))%>%
  dplyr::select(-dummy) -> ds

ds$IP_Total = ip_total
ds$INPUT_Total = input_total

ds$fold = (ds$IP/ds$IP_Total)/((ds$INPUT+1)/ds$INPUT_Total)
ds$logP = -phyper(ds$IP-1,ds$IP_Total,ds$INPUT_Total,ds$IP+ds$INPUT,lower.tail = FALSE,log.p = T)
saveRDS(ds,file=opt$out_rds)
