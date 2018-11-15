args <- commandArgs(T)
dir <- args[1]
NameList <- unlist(strsplit(args[2], split = ":"))

library("dplyr")
library("ggplot2")
library("reshape2")
library("RColorBrewer")

## output 
pdf_cdf <- args[3]
tsv_skew <- args[4]
Lib <- lapply(NameList,function(x){
    read.table(x, header = F,sep = "\t")})
name <- unlist(lapply(NameList, function(x) {
    unlist(strsplit(basename(x),split=".screen_A.csv"))[1]}))
Lib_df <- Reduce(function(x, y) merge(x, y, by=c("V1","V2","V3")),Lib)
colnames(Lib_df) <- c("Gene", "SgRNA", "Seq", name)
Lib_df_m <- melt(Lib_df,id.vars = c("Gene", "SgRNA", "Seq"))
colnames(Lib_df_m) <- c("Gene", "SgRNA", "Seq","Sample","CountNum")
Normlize_df <- Lib_df_m %>% group_by(Sample) %>% mutate(
    NormCount=CountNum/(sum(CountNum)/sum(CountNum)[1])+1)

## run_qc data
sap_mat <- Lib_df[,c(name)]
avg_count <- colMeans(sap_mat)
skew_ratio <- apply(sap_mat,2,function(x){
    return(quantile(x,prob=c(0.9))/quantile(x,prob=c(0.1)))
})
df_meta <- as.data.frame(cbind(avg_count,skew_ratio))
write.table(df_meta,tsv_skew,sep = "\t",quote = F)

## CDF
munual_col <- brewer.pal(9,"Set1")
p <- ggplot(Normlize_df, aes(x=log10(NormCount),y=ecdf(log10(Normlize_df$NormCount))(log10(Normlize_df$NormCount)),colour = Sample)) + 
    stat_ecdf(geom = "step") +
    scale_x_continuous(breaks = seq(0,max(log10(Normlize_df$NormCount)),1)) +
    scale_y_continuous(breaks = seq(0,1,0.1)) +
    scale_color_manual(values= munual_col) +
    labs(x='log10(Normlized count num)', y='CDF') +
    theme_bw() +
    theme(
      axis.text.y = element_text(size=12),
      axis.text.x = element_text(size=12))
pdf(pdf_cdf, width = 6, height = 4)
print(p)
dev.off()
