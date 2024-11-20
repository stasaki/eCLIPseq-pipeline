
options(stringsAsFactors = FALSE)
library(tidyverse)
library(data.table)
library(rtracklayer)
library(plyranges)
library("optparse")
option_list = list(
  make_option(c("-b", "--bed"), type="character", default=NULL, 
              help="clipper bed file", metavar="character"),
  make_option(c("-r", "--ref"), type="character", default=NULL, 
              help="reference gtf file", metavar="character"),
  make_option(c("-p", "--premRNA"), type="character", default=NULL, 
              help="reference gtf file", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, 
              help="SAF output", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default="False", 
              help="mapping to genome", metavar="character")
)

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

#premRNA = "False"
#bed_file ="len100.demultiplex.Aligned.out.sorted.transcriptome.clipper.bed"
#gtf_file = "/mnt/new_disk/ROSMAP_eCLIP_pilot/pipeline_reference/DLPFC_major_transcripts.gtf"
#gtf_file ="/mnt/new_disk/ROSMAP_eCLIP_pilot/pipeline_reference/DLPFC_major_genes.gtf"
#Rscript --vanilla convert_clipper_bed_to_saf.R -b len100.demultiplex.Aligned.out.sorted.transcriptome.clipper.bed -r /mnt/new_disk/ROSMAP_eCLIP_pilot/pipeline_reference/DLPFC_major_transcripts.gtf -p False Rscript --vanilla convert_clipper_bed_to_saf.R -b len100.demultiplex.Aligned.out.sorted.transcriptome.clipper.bed -r /mnt/new_disk/ROSMAP_eCLIP_pilot/pipeline_reference/DLPFC_major_transcripts.gtf -p False

bed_file = opt$bed
gtf_file = opt$ref
premRNA = opt$premRNA
genome = opt$genome

gtf=import(gtf_file)
if(premRNA=="False"){
  # transcript
  gtf=gtf%>%
    filter(type=="exon")
  if(genome=="False"){
    gtf$gene_id=gtf$transcript_id
  }
}else{
  # premRNA
  gtf=gtf%>%
    filter(type=="gene")
}


bed = fread(text = bed_file)
bed$GeneID = bed$V4%>%str_split(string = .,pattern = "_")%>%lapply(.,function(x){paste0(x[1],"_",x[2])})%>%unlist()

saf = bed%>%
  dplyr::select(GeneID,Chr=V1,Start=V2,End=V3,Strand=V6)%>%
  mutate(Start=Start+1) # change to 1-based both inclusive
#GeneID		Chr	Start	End	Strand
#497097		chr1	3204563	3207049	-

# remove peaks located at introns
saf$gene_id=saf$GeneID%>%gsub("_.+","",.)
res = findOverlaps(GRanges(saf),gtf)

data.frame(saf_id = saf$GeneID[queryHits(res)], gtf_id = gtf$gene_id[subjectHits(res)])%>%
  unique()%>%#filter(saf_id%in%remove_saf_id)
  group_by(saf_id)%>%
  summarise(N=dplyr::n())%>%
  filter(N>1)%>%.$saf_id -> remove_saf_id

indx = saf$gene_id[queryHits(res)] == gtf$gene_id[subjectHits(res)]
saf = saf[saf$GeneID%in%c(saf$GeneID[queryHits(res)][indx]),]
saf = saf[!saf$GeneID%in%remove_saf_id,]

if(IRanges::coverage(GRanges(saf[saf$GeneID%in%c(saf$GeneID[queryHits(res)][indx]),])%>%
                     filter(strand=="-"))%>%max%>%max >1){
  stop()
}

if(IRanges::coverage(GRanges(saf[saf$GeneID%in%c(saf$GeneID[queryHits(res)][indx]),])%>%
                     filter(strand=="+"))%>%max%>%max >1){
  stop()
}

data.table::fwrite(saf,file = opt$out,
                   scipen =999,
                   append = F,quote = F,sep = "\t",row.names = F,col.names = T)

