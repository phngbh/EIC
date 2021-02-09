#!/usr/bin/env Rscript

library(ggplot2)
library(dplyr)
library(biomaRt)
library(stringr)

#Load 100 SNPs
data_snp = read.csv("100snps.txt", sep = "\t", header = F)

#Get info from biomart
snp_mart = useMart(biomart="ENSEMBL_MART_SNP", 
                   host="grch37.ensembl.org", path="/biomart/martservice", 
                   dataset="hsapiens_snp")
data_snp_mart = getBM(attributes=c(
  "refsnp_id", "chr_name","allele", "allele_1", "minor_allele",
  "minor_allele_freq", "synonym_name", "variation_names"),
  filters="snp_filter", values=data_snp$V1,
  mart=snp_mart, uniqueRows=TRUE)
  
#Get coordinates for 100 SNPs (for genotype imputation using other tools)
data_snp_coord = getBM(attributes=c(
  "refsnp_id", "chr_name", "chrom_start", "chrom_end"),
  filters="snp_filter", values=rsID_ethc,
  mart=snp_mart, uniqueRows=TRUE)
data_snp_coord = data_snp_coord[-grep("H",data_snp_coord$chr_name),]
write.table(data_snp_coord, file = "snp.txt", sep = "\t", row.names = F, quote = F)