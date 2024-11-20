options(stringsAsFactors = FALSE)
library(tidyverse)
library(data.table)
library(rtracklayer)
library(plyranges)
library("optparse")

option_list = list(
  make_option(c("-i", "--input_count"), type="character", default=NULL, 
              help="input count", metavar="character"),
  make_option(c("-s", "--input_summary"), type="character", default=NULL, 
              help="input summary", metavar="character"),
  make_option(c("-p", "--ip_count"), type="character", default=NULL, 
              help="ip count", metavar="character"),
  make_option(c("-m", "--ip_summary"), type="character", default=NULL, 
              help="ip summary", metavar="character"),
  make_option(c("-r", "--out_rds"), type="character", default=NULL, 
              help="fisher output", metavar="character"),
  make_option(c("-f", "--out_bed_fdr"), type="character", default=NULL, 
              help="fisher output with only FDR cutoff", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#premRNA = "False"
#bed_file ="len100.demultiplex.Aligned.out.sorted.transcriptome.clipper.bed"
#gtf_file = "/mnt/new_disk/ROSMAP_eCLIP_pilot/pipeline_reference/DLPFC_major_transcripts.gtf"
#gtf_file ="/mnt/new_disk/ROSMAP_eCLIP_pilot/pipeline_reference/DLPFC_major_genes.gtf"
#Rscript --vanilla convert_clipper_bed_to_saf.R -b len100.demultiplex.Aligned.out.sorted.transcriptome.clipper.bed -r /mnt/new_disk/ROSMAP_eCLIP_pilot/pipeline_reference/DLPFC_major_transcripts.gtf -p False Rscript --vanilla convert_clipper_bed_to_saf.R -b len100.demultiplex.Aligned.out.sorted.transcriptome.clipper.bed -r /mnt/new_disk/ROSMAP_eCLIP_pilot/pipeline_reference/DLPFC_major_transcripts.gtf -p False
input_count = opt$input_count
input_summary = opt$input_summary
ip_count = opt$ip_count
ip_summary = opt$ip_summary

input = fread(input_count)
colnames(input)[7]="INPUT"

ip = fread(ip_count)
colnames(ip)[7]="IP"

input_summary=fread(input_summary)
colnames(input_summary)[2]="Reads"
input_total=input_summary$Reads%>%sum
ip_summary=fread(ip_summary)
colnames(ip_summary)[2]="Reads"
ip_total=ip_summary$Reads%>%sum

ds = inner_join(ip,input,by=c("Geneid","Chr","Start","End","Strand","Length"))
ds$IP_Total = ip_total
ds$INPUT_Total = input_total

rm(ip)
rm(input)

ds$fold = (ds$IP/ds$IP_Total)/(ds$INPUT/ds$INPUT_Total)
ds$logP = -phyper(ds$IP-1,ds$IP_Total,ds$INPUT_Total,ds$IP+ds$INPUT,lower.tail = FALSE,log.p = T)
#ds$p.adj = qvalue::qvalue(exp(-ds$logP))$qvalues
#ds$p.adj = stats::p.adjust(p=exp(-ds$logP),method = "BH")

ds$p.adj =check = tryCatch({
  qvalue::qvalue(exp(-ds$logP))$qvalues
}, error = function(e){
  stats::p.adjust(p=exp(-ds$logP),method = "BH")
})

ds%>%
  filter(exp(-logP)<0.05,p.adj<0.05)%>%
  mutate(Dummy=0)%>%
  dplyr::select(Chr,Start,End,Geneid,Dummy,Strand)%>%
  mutate(Start=Start-1)%>%
  data.table::fwrite(.,file = opt$out_bed_fdr,scipen =999,
                     append = F,quote = F,sep = "\t",row.names = F,col.names = F)
saveRDS(ds,file=opt$out_rds)
