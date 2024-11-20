
library(tidyverse)
library(rtracklayer)
library(GenomicRanges)
library(data.table)
library(plyranges)
library(yaml)
options(stringsAsFactors = FALSE)
dir.create("clipper/consensus")

config = yaml::read_yaml(file = "configs/config.yaml")

# collect peak info #
lapply(names(config$iids), function(x){
  paste0("clipper/demultiplex.Aligned.out.sorted.clipper.fisher.rds")%>%
    readRDS(.)%>%
    filter(exp(-logP)<0.05,p.adj<0.05)%>%
    dplyr::select(Chr,Start,End,Geneid,Strand,fold)%>%
    mutate(sample=x)%>%
    return()
})%>%
  bind_rows()->  radc_peaks

lapply(c(20,40), function(min_samples){
  
  lapply(c("+","-"), function(st){
    radc_peaks = radc_peaks%>%
      filter(Strand==st)%>%
      group_by(sample)%>%
      nest()%>%
      rowwise()%>%
      mutate(data=list(GRanges(data)))
    names(radc_peaks$data)=radc_peaks$sample
    peak_grangeslist = GRangesList(radc_peaks$data)
    peak_coverage <- GenomicRanges::coverage(peak_grangeslist)
    covered_ranges <- IRanges::slice(peak_coverage, lower=min_samples, rangesOnly=T)
    covered_granges <- GRanges(covered_ranges)
    
    # add gene_id
    radc_peaks = radc_peaks%>%
      rowwise()%>%
      mutate(data=list(as.data.frame(data)))%>%
      unnest()%>%
      GRanges()
    indx = findOverlaps(query = covered_granges,subject = radc_peaks)
    data.frame(rowid = queryHits(indx),
               Geneid=radc_peaks$Geneid[subjectHits(indx)],
               fold=radc_peaks$fold[subjectHits(indx)])%>%
      group_by(rowid)%>%
      summarise(Geneid=list(Geneid%>%gsub("_.+","",.)%>%unique()),
                fold=mean(fold))%>%
      rowwise()%>%
      mutate(Geneid_n=length(Geneid))%>%
      filter(Geneid_n==1)%>%
      unnest()%>%
      dplyr::select(-Geneid_n)%>%
      unique() -> add_gene_id
    covered_granges$gene_id = NA
    covered_granges$gene_id[add_gene_id$rowid]=add_gene_id$Geneid
    covered_granges$fold[add_gene_id$rowid]=add_gene_id$fold
    
    # make SAF file
    covered_granges%>%
      as.data.frame()%>%
      group_by(gene_id)%>%
      mutate(GeneID=paste0(gene_id,"_",1:n()),
             Strand=st)%>%
      dplyr::select(GeneID,Chr=seqnames,	Start=start,	End=end,	Strand,	gene_id,fold) -> saf
    
    return(saf)
  })%>%
    bind_rows()-> saf
  
  return(saf)
}) -> consensus_peaks

# extend peak range
consensus_peaks_min20=consensus_peaks[[1]]
consensus_peaks_min40=consensus_peaks[[2]]

res_fo = findOverlaps(query = consensus_peaks_min40%>%
                        GRanges(),subject = consensus_peaks_min20%>%
                        GRanges(),type = "within")
consensus_peaks_min40_ext = consensus_peaks_min20[subjectHits(res_fo)%>%unique(),]
consensus_peaks_min40_ext%>%
  dplyr::select(-fold)%>%
  data.table::fwrite(.,file = "clipper/consensus/consensus_peaks_min40ext.saf",
                     scipen =999,
                     append = F,quote = F,sep = "\t",row.names = F,col.names = T)
consensus_peaks_min40_ext%>%
  saveRDS(.,file = "clipper/consensus/consensus_peaks_min40ext.rds")
