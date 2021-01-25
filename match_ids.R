#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

map = read.csv(args[1],header = F, sep = "\t")
snplist = read.csv(args[2],header = F, sep = "\t")
index = match(snplist$V3,map$V4)
snplist$id = map$V2[index]
snplist = na.omit(snplist)
write.table(snplist, file = paste0("./scratch/snp_",args[3],"_withIDs.txt"), row.names = F, col.names = F, quote = F)