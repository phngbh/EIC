#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

sample = read.table(args[1], header = F, sep = "\t", stringsAsFactors = F)
ref = read.table(args[2], header = F, stringsAsFactors = F)
ref = ref[match(sample$V4, ref$V1),]

for (i in 1:nrow(sample)){
  sam_ref = sample$V6[i]
  sam_alt = sample$V5[i]
  ref_ref = ref$V2[i]
  ref_alt = ref$V3[i]
  
  if (sam_ref == ref_ref && sam_alt == ref_alt){
    next
  } else if ( sam_ref == ref_alt && sam_alt == ref_ref){
    next
  } else {
    cat(paste(sample[i,2]), file="./qc/unmatched_alleles.txt", append=TRUE, sep = "\n")
  }
}