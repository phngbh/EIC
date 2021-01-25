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
write.table(data_snp_coord, file = "snp.txt", sep = "\t", row.names = F)

#Import imputed data
for (chr in seq(1:21)){
  file=paste0("snp_imputed_",chr,".txt")
  if (file.exists(file)){
    assign(paste0("snp_",chr,"_imputed"),read.delim(file, header = T, sep = "\t", row.names = NULL))
  }
}
cname=colnames(snp_1_imputed)
for (chr in seq(1:21)){
  name=paste0('snp_',chr,'_imputed')
  if (exists(name)) {
    obj=get(name)
    obj=obj[,cname]
    assign(name,obj)
  }
}
snp_5_imputed = snp_5_imputed[-4,]
snp_11_imputed = snp_11_imputed[-2,]
snp_imputed= do.call(rbind, lapply( paste0("snp_",c(1:6,8:20),"_imputed") , get) )
snp_imputed$ID= as.character(snp_imputed$ID)
index=match(snp_imputed$POS,data_snp_coord$chrom_start)
for (i in 1:length(index)){
  snp_imputed$ID[i]=data_snp_coord$refsnp_id[index[i]]
}

snp_imputed= snp_imputed[,-c(1,2,4:9)]
rownames(snp_imputed)=snp_imputed$ID
data_final= as.matrix(snp_imputed[,-1])

#Label cell lines with estimated ethnicity probability using frequencies from 1000Gen
size_1kg = read.csv("sample_size_1kg.csv", header = T)
data_1kg = read.delim("genotypes2.txt", header = F, sep = "",
                         col.names = c("rsID", "population", "genotype", "frequency"))
data_1kg$sup_pop = vector(mode = "character", length = nrow(data_1kg))
data_1kg$size= vector(mode = "numeric", length = nrow(data_1kg))
index = match(data_1kg$population, size_1kg$Population)

for (i in 1:nrow(data_1kg)){
  data_1kg$size[i]=size_1kg$Number.of.genotypes[index[i]]
  if (data_1kg$population[i]=="ALL"){
    data_1kg$sup_pop[i]="ALL"
  } else if (data_1kg$population[i] %in% c("ACB","AFR","ASW","ESN","GWD","LWK","MSL","YRI")){
    data_1kg$sup_pop[i]="AFR"
  } else if (data_1kg$population[i] %in% c("AMR","CLM","MXL","PEL","PUR")){
    data_1kg$sup_pop[i]="AMR"
  } else if (data_1kg$population[i] %in% c("EAS","CDX","CHB","CHX","JPT","KHV")){
    data_1kg$sup_pop[i]="EAS"
  } else if (data_1kg$population[i] %in% c("EUR","CEU","FIN","GBR","IBS","TSI")){
    data_1kg$sup_pop[i]="EUR"
  } else if (data_1kg$population[i] %in% c("SAS","BEB","GIH","ITU","PJL","STU")){
    data_1kg$sup_pop[i]="SAS"
  } else {
    data_1kg$sup_pop[i]="NA"
  }
}
data_1kg=data_1kg[-grep("NA",data_1kg$sup_pop),] 
data_1kg=data_1kg[!data_1kg$population %in% c("ALL","AFR","AMR","EAS","EUR","SAS"),] %>%
  droplevels()

prob = vector(mode = "list", length = ncol(data_final))
names(prob) = colnames(data_final)

for (n in 1:ncol(data_final)){
  prob[[n]] = matrix(nrow = nrow(data_final), ncol = nlevels(data_1kg$population),
                     dimnames = list(rownames(data_final),levels(data_1kg$population)))
  for (i in 1:nrow(data_final)) {
    for (j in 1:ncol(prob[[n]])) {
      sub_data = data_1kg[as.character(data_1kg$rsID) == rownames(data_final)[i] &
                             as.character(data_1kg$population) == colnames(prob[[n]])[j], ]
      temp_index = match(data_final[i,n], sub_data$genotype)
      if (is.na(temp_index)){
        tmp = strsplit(data_final[i,n],NULL)[[1]]
        tmp_rev = rev(tmp)
        data_final[i,n] = paste(tmp_rev,collapse = '')
        temp_index = match(data_final[i,n], sub_data$genotype)
      }
      temp_snp_freq = sub_data$frequency[temp_index]
      temp_pop_size = sub_data$size[temp_index]
      prob[[n]][i,j] = as.numeric(temp_snp_freq) * as.numeric(temp_pop_size)  
    } 
  }
  prob[[n]] = apply(prob[[n]], 2, prod, na.rm=T) 
  prob[[n]] = names(which.max(prob[[n]]))
}

prob_vec = unlist(prob)
cell_ancestry = data.frame(cell=names(prob_vec), ancestry=prob_vec)
cell_ancestry$sup_pop=match(cell_ancestry$ancestry,data_1kg$population)
cell_ancestry$sup_pop=data_1kg$sup_pop[cell_ancestry$sup_pop]
cell_ancestry$cell = str_replace_all(cell_ancestry$cell, c("\\."="-", "^X"="")) 
cell_info=read.csv("cell_line.csv",header = T)
cell_ancestry$cancer=match(cell_ancestry$cell,cell_info$Name)
cell_ancestry$cancer=cell_info$Cancer[cell_ancestry$cancer]




