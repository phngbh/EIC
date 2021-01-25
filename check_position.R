#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

sample = read.table(args[1], header = F, sep = "\t")
ref = read.table(args[2], header = F)
index = match(sample$V4, ref$V1)
index_unmatched = is.na(index)  
unmatched_position = sample[index_unmatched,] 
write.table(unmatched_position, file = "qc/unmatched_position.txt", row.names = F, col.names = F, quote = F)