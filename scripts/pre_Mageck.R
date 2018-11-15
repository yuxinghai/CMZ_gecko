#!/usr/bin/env Rscript
args <- commandArgs(T)
## input 
dir <- args[1]
cond <- unlist(strsplit(args[2], split = ":"))

## output
sampleInfo <- args[3]
tsv_count <- args[4]

sample_Info <- read.csv(sampleInfo,header = T,sep=",", stringsAsFactors=F)
saps <- unique(sample_Info[sample_Info$condition %in% cond,]$samples)
NameList <- file.path(dir, paste("count/",saps,".screen_A.csv",sep=""))
Lib <- lapply(NameList,function(x){
    read.table(x, header = F,sep = "\t")})
name <- unlist(lapply(NameList, function(x) {
    unlist(strsplit(basename(x),split=".screen_A.csv"))[1]}))
Lib_df <- Reduce(function(x, y) merge(x, y, by=c("V1","V2","V3")),Lib)
colnames(Lib_df) <- c("gene", "sgRNA", "seq", name )
Lib_df_final <- Lib_df[,c("sgRNA", "gene", name)]
write.table(Lib_df_final, tsv_count, sep = "\t", row.names = F, quote = F)

