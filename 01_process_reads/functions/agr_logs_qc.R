

library(tidyverse)
library(data.table)
library("optparse")
options(stringsAsFactors = FALSE)

option_list = list(
  make_option(c("-s", "--star_log"), type="character", default=NULL, 
              help="star_log", metavar="character"),
  make_option(c("-d", "--dedup_log"), type="character", default=NULL, 
              help="dedup_log", metavar="character"),
  make_option(c("-f", "--fillter_mapped_log"), type="character", default=NULL, 
              help="fillter_mapped_log", metavar="character"),
  make_option(c("-o", "--out_rds"), type="character", default=NULL, 
              help="output", metavar="character"),
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

ds = read_lines(opt$star_log)%>%
  .[grepl("\t",.)]%>%
  fread(text=.,sep = "\t",header = F)%>%
  mutate(V1=gsub("\\|","",V1))%>%
  mutate(V1=str_trim(V1))%>%
  filter(!grepl(":",V2))%>%
  mutate(V2=gsub("%","",V2)%>%as.numeric())%>%
  dplyr::rename(Stat=V1,Value=V2)%>%
  mutate(task_name = "star",
         file_name = basename(opt$star_log))

read_lines(opt$dedup_log)%>%
  .[grepl("Reads:|out:|deduplicated:|position:",.)]%>%
  gsub(".+INFO ","",.)%>%
  fread(text=.,sep = ":",header = F)%>%
  dplyr::rename(Stat=V1,Value=V2)%>%
  mutate(task_name = "dedup",
         file_name = basename(opt$dedup_log))%>%
  bind_rows(ds,.) -> ds

read_lines(opt$fillter_mapped_log)%>%.[-1]%>%
  fread(text=.,sep = "\t",header = F)%>%
  dplyr::rename(Stat=V1,Value=V2)%>%
  mutate(task_name = "fillter_mapped",
         file_name = basename(opt$fillter_mapped_log))%>%
  bind_rows(ds,.) -> ds


saveRDS(ds,file=opt$out_rds)

