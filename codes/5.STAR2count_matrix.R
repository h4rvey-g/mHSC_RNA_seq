# used the get and merge the raw counts from the STAR outputs
# call DEGs using DESeq2
library("data.table")

indir = "3.STARalign_clean"
outdir = "5.DEG"
dir.create(outdir, showWarnings = T, recursive = T)
#get the sample infomation
sample_meta = fread("samples.txt", header = F)
sample_names = sample_meta$V4

df = fread(paste0(indir, "/",sample_names[1],"/",sample_names[1],"ReadsPerGene.out.tab"))
countData = df[,c(1, 2)]
for (i in 2:length(sample_names)) {
  df = fread(paste0(indir, "/",sample_names[i],"/",sample_names[i],"ReadsPerGene.out.tab"))
  count = df[,2]
  countData = cbind(countData, count)
}
###Skip first 4 lines, count data starts on the 5th line
countData = countData[c(5:nrow(countData)), ]

# Noet: all STAR outputs have the same GeneID order defined by the STARindex file “geneInfo.tab”
# so we can cbind them directly and use the index file to get the geneNames

geneInfo = fread("STARindex/geneInfo.tab")

countData = cbind(geneInfo, countData[, -1])

colnames(countData) = c("GeneID", "GeneName", "GeneType", sample_names)

# remove the genes with counts < 5 in more than half of the samples
filter = apply(countData[,-c(1:3)] , 1, function(x) sum(x >= 5) >= length(sample_names)*0.5)
countData = countData[filter, ]

# output the raw count matrix
write.table(
  countData,
  file = paste0(outdir, "/star2counts.tsv") ,
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = T
)

# convert count matrix to TPM FPKM etc

# get the gene lenght

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GenomicFeatures")
library(GenomicFeatures)


## make TxDb from GTF file 
txdb <- makeTxDbFromGFF('gencode.vM31.primary_assembly.annotation.gtf')

## get gene information
all.genes <- genes(txdb)

## import your list of gene names
my.genes <- countData$GeneID

## get the length of each of those genes
my.genes.lengths <- width(all.genes[my.genes])
## put the names back on the lengths
names(my.genes.lengths) <- my.genes

library(DGEobj.utils)
BiocManager::install("edgeR")
library(edgeR)

count_matrix = as.matrix(countData[,-c(1:3)])
tpm = convertCounts(count_matrix, "TPM", my.genes.lengths)
fpkm = convertCounts(count_matrix, "FPKM", my.genes.lengths)
cpm = convertCounts(count_matrix, "CPM", my.genes.lengths)

# output the TPM matrix
write.table(
  cbind(countData[,c(1:3)],as.data.frame(tpm)),
  file = paste0(outdir, "/star2TPM.tsv") ,
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = T
)

# output the RPKM matrix
write.table(
  cbind(countData[,c(1:3)],as.data.frame(fpkm)),
  file = paste0(outdir, "/star2FPKM.tsv") ,
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = T
)

# output the CPM matrix
write.table(
  cbind(countData[,c(1:3)],as.data.frame(cpm)),
  file = paste0(outdir, "/star2CPM.tsv") ,
  quote = F,
  sep = "\t",
  row.names = F,
  col.names = T
)



# DEG analysis
BiocManager::install("DESeq2")
library(DESeq2)

group_list = c("D", "D", "B", "B", "A", "C", "A", "C") # all group labels corresponding to columns of countData

geneID = countData$GeneID
expMat = as.matrix(countData[,-c(1:3)])
rownames(expMat) = geneID

#colData <- data.frame(row.names = colnames(expMat),group_list=group_list)


# A vs B
compare = "B/A"
compare = "D/C"
compare = "C/A"
compare = "D/B"
{
  group = unlist(strsplit(compare,"/"))
  BA = which(group_list %in% group)
  expMat_BA = expMat[,BA]
  colData_BA = data.frame(row.names = colnames(expMat_BA),group_list=group_list[BA])
  
  dds <- DESeqDataSetFromMatrix(countData = expMat_BA,
                                colData = colData_BA,
                                design = ~ group_list)
  dds <- DESeq(dds,quiet = F) 
  
  res <- results(dds,contrast=c("group_list", group))  #指定提取为B/A结果
  
  
  resOrdered <- res[order(res$padj),]
  DEG_DEseq2 <- as.data.frame(resOrdered)
  DEG_DEseq2 <- na.omit(DEG_DEseq2)
  DEG_DEseq2$GeneID = rownames(DEG_DEseq2)
  
  #add gene symbols to the DEG results
  DEG_DEseq2 = merge(DEG_DEseq2, countData[,c(1:2)], by = "GeneID")
  
  # output the CPM matrix
  colname = c("GeneID",	"GeneName", "baseMean",	"log2FoldChange",	"lfcSE",	"stat",	"pvalue",	"padj")
  write.table(
    DEG_DEseq2[, colname],
    file = paste0(outdir, "/DEG_DESeq2_",paste(group, collapse = ""),".tsv") ,
    quote = F,
    sep = "\t",
    row.names = F,
    col.names = T
  )
  
  ###valcano_DEseq2
  library(ggplot2)
  library(ggrepel)
  logFC_t <- with(DEG_DEseq2,mean(abs(log2FoldChange)) + 2*sd(abs(log2FoldChange)))
  logFC_t <- round(logFC_t,2)
  logFC_t = 2 #ifelse(logFC_t >2, 2, logFC_t) 
  
  DEG_DEseq2$change = as.factor(ifelse(DEG_DEseq2$pvalue < 0.05 & abs(DEG_DEseq2$log2FoldChange) > logFC_t,
                                       ifelse(DEG_DEseq2$log2FoldChange > logFC_t, 'UP','DOWN'),'STABLE'))
  stats = table(DEG_DEseq2$change)
  
  pdf(file = paste0(outdir, "/DEG_DESeq2_",paste(group, collapse = ""),"_logFC_",logFC_t,"_plots.pdf"))
  # add top genes to show names in plot
  DEG_DEseq2$label <- ifelse(DEG_DEseq2$pvalue < 0.00001 & abs(DEG_DEseq2$log2FoldChange) >= 2*logFC_t, DEG_DEseq2$GeneName, "")
  p = ggplot(DEG_DEseq2,aes(x=log2FoldChange, y=-log10(pvalue),color=change))+
    geom_point(alpha=0.4,size=2) +
    theme_classic() +
    xlab(paste0("log2FC(",compare,"; auto_cutoff=",logFC_t,")")) +
    ylab("-log10(p-value)") +
    theme(plot.title = element_text(size = 15,hjust = 0.5)) +
    scale_colour_manual(values = c('#33a3dc','#72777b','#c85d44'),
                        labels=paste(paste(names(stats),stats, sep = " ( n="),")")) +
    geom_hline(yintercept = -log10(0.05),lty =2 ) +
    geom_vline(xintercept = c(-logFC_t,logFC_t),lty =2) #+
    # geom_label_repel(data = DEG_DEseq2,aes(label =label),
    #                  size =3,box.padding = unit(2,"lines"),
    #                  point.padding = unit(0.5,"lines"),
    #                  segment.color ="#ffd400",
    #                  show.legend = FALSE,max.overlaps = 20000)
  print(p)
  
  library(pheatmap)
  ###heatmap_DEseq2
  dat <- expMat_BA
  FC <- DEG_DEseq2$log2FoldChange
  names(FC) <- rownames(DEG_DEseq2)
  top = 100 # top 100 DEGs for both directions, in total 2*top genes
  DEG_Deseq2_top <- as.numeric(c(names(head(sort(FC),top)),names(tail(sort(FC),top)))) 
  dat <- t(scale(t(dat[DEG_Deseq2_top,])))
  data <- na.omit(dat)
  group <- data.frame(group=group_list[BA])
  rownames(group)=colnames(dat)
  pheatmap(dat, cluster_cols = T, show_colnames =T,show_rownames =F,
           cluster_rows = T,
           annotation_col =group,
           color =colorRampPalette(c("#00ae9d","white","#E85827"))(100))
  dev.off()
}
