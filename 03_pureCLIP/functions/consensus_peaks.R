

library(tidyverse)
library(rtracklayer)
options(stringsAsFactors = FALSE)

# clipper_tr_pilot ####
outdir ="pureclip_genome_exB2"
dir.create(paste0("result/",outdir))
min_samples = 10

config = yaml::read_yaml(file = "configs/config.yaml")


list.files(path = "pureclip/",
           pattern = "demultiplex.Aligned.out.sorted.rmDup.sorted.featureCounts.stripped_5p_mismatches.sorted.m6A-IP.crosslink_sites_count.fisher.rds",
           recursive = T,full.names = T)%>%.[(dirname(.)%>%dirname()%>%basename())%in%names(config$iids)]%>%
  lapply(., function(x){
    ds <- readRDS(x)
    ds%>%
      mutate(p.adj=p.adjust(exp(-logP),method = "BH"))%>%
      filter(exp(-logP)<0.05,p.adj<0.05)%>%
      dplyr::select(Chr=V1,Start=V2,End=V3,Strand=V6)%>%
      mutate(Geneid=Chr,
             sample=dirname(x)%>%dirname()%>%basename())%>%
      return()
  })%>%bind_rows() -> ds

# add geneid
"gsutil cp gs://pipeline_references/release_27/gencode.v27.primary_assembly.annotation.gtf /mnt/new_disk/tmp/"%>%
  system()
gtf = import("/mnt/new_disk/tmp/gencode.v27.primary_assembly.annotation.gtf")
gtf = gtf%>%
  as.data.frame()%>%
  filter(type=="exon")%>%
  ungroup()%>%
  dplyr::select(seqnames,start,end, strand,gene_id)%>%
  GRanges()

indx=findOverlaps(GRanges(ds),gtf,ignore.strand=FALSE)
indx%>%as.data.frame()%>%
  mutate(subjectHits=gtf$gene_id[subjectHits])%>%
  group_by(queryHits)%>%
  summarise(subjectHits=list(unique(subjectHits)),
            subjectHits_n=length(subjectHits%>%unlist()))%>%
  filter(subjectHits_n==1)%>%
  unnest() -> indx

ds$Geneid=NA
ds$Geneid[indx$queryHits]=indx$subjectHits
# ds = ds%>%
#   filter(!is.na(Geneid))
# ds%>%
#   filter(is.na(Geneid))%>%.$Chr%>%table
rm(gtf)
rm(indx)
WGCNA::collectGarbage()

# make BED file
ds%>%
  group_by(Chr,Start,End,Strand,Geneid)%>%
  summarise(N=dplyr::n(),source=paste0(sample,collapse = ";"))%>%
  filter(N>=min_samples)%>%
  group_by(Geneid)%>%
  mutate(Geneid=paste0(Geneid,"_",1:dplyr::n()))%>%
  dplyr::select(Chr, Start, End, Geneid,N, Strand,source)%>%
  data.table::fwrite(.,file = paste0("result/",outdir,"/consensus_peaks.bed"),
                     scipen =999,
                     append = F,quote = F,sep = "\t",row.names = F,col.names = F)

paste0("/usr/bin/gsutil cp result/",outdir,"/consensus_peaks.bed gs://shinya_test/",
       config$output_dir,"consensus/",outdir,"/consensus_peaks.bed")%>%system()


for (N_th in c(10,20,30,40,50,60)){
  ds%>%
    group_by(Chr,Start,End,Strand,Geneid)%>%
    summarise(N=dplyr::n(),source=paste0(sample,collapse = ";"))%>%
    filter(N>=N_th)%>%
    group_by(Geneid)%>%
    mutate(Geneid=paste0(Geneid,"_",1:dplyr::n()))%>%
    dplyr::select(Chr, Start, End, Geneid,N, Strand,source)%>%
    data.table::fwrite(.,file = paste0("result/",outdir,"/consensus_peaks_",N_th,".bed"),
                       scipen =999,
                       append = F,quote = F,sep = "\t",row.names = F,col.names = F)
  
  paste0("gsutil cp result/",outdir,"/consensus_peaks_",N_th,".bed gs://shinya_test/",
         config$output_dir,"consensus/",outdir,"/consensus_peaks_",N_th,".bed")%>%system()
}

# # looked only around TSS ####
# library(Guitar)
# bed_file="result/pureclip_genome_exB2/consensus_peaks_60.bed"
# bed = data.table::fread(bed_file)%>%
#   mutate(V2=V2+1)%>%
#   dplyr::select(seqnames=V1,start=V2,end=V3,strand=V6)
# 
# n_sample = min(c(10000,nrow(bed)))
# bed = bed%>%
#   sample_n(n_sample)%>%
#   #filter(seqnames %in% c(gtf%>%as.data.frame()%>%filter(strand=="+")%>%.$transcript_id%>%as.character()%>%unique()))%>%
#   GRanges()
# 
# system("gsutil cp gs://pipeline_references/release_27/gencode.v27.primary_assembly.annotation.gtf /mnt/new_disk/tmp/")
# major_transcripts.gtf = "/mnt/new_disk/tmp/gencode.v27.primary_assembly.annotation.gtf"
# p <- GuitarPlot(txGTF =major_transcripts.gtf,
#                 stGRangeLists = GRangesList(bed),
#                 miscOutFilePrefix = NA,
#                 txGuitarTxdbSaveFile="gencode.v27.primary_assembly.annotation")
# 
# library(plyranges)
# load("GuitarTxdb-gencode.v27.primary_assembly.annotation-20220825")
# 
# library(parallel)
# mclapply(names(guitarTxdb$tx$tx), function(gene){
#   
#   gene_range = guitarTxdb$tx$tx[[gene]]
#   
#   if(strand(gene_range)%>%as.character()%>%.[1] =="+"){
#     window_start = 901
#     window_end = 1100
#   }else{
#     window_start = width(gene_range)%>%sum - 1100
#     window_end = width(gene_range)%>%sum - 901
#   }
#   
#   gene_range%>%width%>%cumsum()  -> cs_vec
#   cs_vec_end = cs_vec - window_end
#   bin_one_before_end = which(cs_vec_end<0)%>%max()
#   if(is.infinite(bin_one_before_end)){
#     # at the first bin
#     window_to_take_end = window_end
#     bin_one_before_end=0
#   }else{
#     window_to_take_end = abs(cs_vec_end[bin_one_before_end])
#   }
#   
#   cs_vec_start = cs_vec - window_start
#   bin_one_before_start = which(cs_vec_start<0)%>%max()
#   if(is.infinite(bin_one_before_start)){
#     # at the first bin
#     window_to_take_start = window_start
#     bin_one_before_start=0
#   }else{
#     window_to_take_start = abs(cs_vec_start[bin_one_before_start])
#   }
#   
#   if(bin_one_before_start == bin_one_before_end ){
#     gene_range[bin_one_before_start+1]%>%
#       mutate(end=start+window_to_take_end,
#              start=start+window_to_take_start)%>%
#       mutate(Geneid = gene)%>%
#       as.data.frame()%>%
#       return()
#   }else if ((bin_one_before_start+1)==bin_one_before_end){
#     c(gene_range[bin_one_before_start+1]%>%
#         mutate(start=start+window_to_take_start),
#       gene_range[bin_one_before_end+1]%>%
#         mutate(end=start+window_to_take_end))%>%
#       mutate(Geneid = gene)%>%
#       as.data.frame()%>%
#       return()
#   }else{
#     c(gene_range[bin_one_before_start+1]%>%
#         mutate(start=start+window_to_take_start),
#       gene_range[(bin_one_before_start+2):bin_one_before_end],
#       gene_range[bin_one_before_end+1]%>%
#         mutate(end=start+window_to_take_end))%>%
#       mutate(Geneid = gene)%>%
#       as.data.frame()%>%
#       return()
#   }
# },mc.cores = 4)%>%bind_rows() -> tss_range
# saveRDS(tss_range,file="result/tss_range.rds")
# 
# indx = findOverlaps(GRanges(ds),GRanges(tss_range))
# indx = unique(queryHits(indx))
# for (N_th in c(10,20,30,40,50,60)){
#   ds[indx,]%>%
#     group_by(Chr,Start,End,Strand,Geneid)%>%
#     summarise(N=dplyr::n(),source=paste0(sample,collapse = ";"))%>%
#     filter(N>=N_th)%>%
#     group_by(Geneid)%>%
#     mutate(Geneid=paste0(Geneid,"_",1:dplyr::n()))%>%
#     dplyr::select(Chr, Start, End, Geneid,N, Strand,source)%>%
#     data.table::fwrite(.,file = paste0("result/",outdir,"/consensus_peaks_TSS_",N_th,".bed"),
#                        scipen =999,
#                        append = F,quote = F,sep = "\t",row.names = F,col.names = F)
#   
#   paste0("gsutil cp result/",outdir,"/consensus_peaks_TSS_",N_th,".bed gs://shinya_test/",
#          config$output_dir,"consensus/",outdir,"/consensus_peaks_TSS_",N_th,".bed")%>%system()
# }
# 
# # looked only around the end of CDS ####
# 
# load("GuitarTxdb-gencode.v27.primary_assembly.annotation-20220825")
# mclapply(names(guitarTxdb$mrna$tx), function(gene){
#   
#   gene_range = guitarTxdb$mrna$tx[[gene]]
#   
#   if(strand(gene_range)%>%as.character()%>%.[1] =="+"){
#     window_start = guitarTxdb$mrna$componentWidth[gene,1:3]%>%sum() -100 + 1
#     window_end = window_start+200 -1 
#   }else{
#     window_end = width(gene_range)%>%sum - (guitarTxdb$mrna$componentWidth[gene,1:3]%>%sum() -100 + 1)
#     window_start = window_end - 200 + 1
#   }
#   
#   gene_range%>%width%>%cumsum()  -> cs_vec
#   cs_vec_end = cs_vec - window_end
#   bin_one_before_end = which(cs_vec_end<0)%>%max()
#   if(is.infinite(bin_one_before_end)){
#     # at the first bin
#     window_to_take_end = window_end
#     bin_one_before_end=0
#   }else{
#     window_to_take_end = abs(cs_vec_end[bin_one_before_end])
#   }
#   
#   cs_vec_start = cs_vec - window_start
#   bin_one_before_start = which(cs_vec_start<0)%>%max()
#   if(is.infinite(bin_one_before_start)){
#     # at the first bin
#     window_to_take_start = window_start
#     bin_one_before_start=0
#   }else{
#     window_to_take_start = abs(cs_vec_start[bin_one_before_start])
#   }
#   
#   if(bin_one_before_start == bin_one_before_end ){
#     gene_range[bin_one_before_start+1]%>%
#       mutate(end=start+window_to_take_end,
#              start=start+window_to_take_start)%>%
#       mutate(Geneid = gene)%>%
#       as.data.frame()%>%
#       return()
#   }else if ((bin_one_before_start+1)==bin_one_before_end){
#     c(gene_range[bin_one_before_start+1]%>%
#         mutate(start=start+window_to_take_start),
#       gene_range[bin_one_before_end+1]%>%
#         mutate(end=start+window_to_take_end))%>%
#       mutate(Geneid = gene)%>%
#       as.data.frame()%>%
#       return()
#   }else{
#     c(gene_range[bin_one_before_start+1]%>%
#         mutate(start=start+window_to_take_start),
#       gene_range[(bin_one_before_start+2):bin_one_before_end],
#       gene_range[bin_one_before_end+1]%>%
#         mutate(end=start+window_to_take_end))%>%
#       mutate(Geneid = gene)%>%
#       as.data.frame()%>%
#       return()
#   }
# },mc.cores = 4)%>%bind_rows() -> cds_range
# saveRDS(cds_range,file="result/cds_range.rds")
# 
# indx = findOverlaps(GRanges(ds),GRanges(cds_range))
# indx = unique(queryHits(indx))
# for (N_th in c(10,20,30,40,50,60)){
#   ds[indx,]%>%
#     group_by(Chr,Start,End,Strand,Geneid)%>%
#     summarise(N=dplyr::n(),source=paste0(sample,collapse = ";"))%>%
#     filter(N>=N_th)%>%
#     group_by(Geneid)%>%
#     mutate(Geneid=paste0(Geneid,"_",1:dplyr::n()))%>%
#     dplyr::select(Chr, Start, End, Geneid,N, Strand,source)%>%
#     data.table::fwrite(.,file = paste0("result/",outdir,"/consensus_peaks_CDS_",N_th,".bed"),
#                        scipen =999,
#                        append = F,quote = F,sep = "\t",row.names = F,col.names = F)
#   
#   paste0("gsutil cp result/",outdir,"/consensus_peaks_CDS_",N_th,".bed gs://shinya_test/",
#          config$output_dir,"consensus/",outdir,"/consensus_peaks_CDS_",N_th,".bed")%>%system()
# }
# 
