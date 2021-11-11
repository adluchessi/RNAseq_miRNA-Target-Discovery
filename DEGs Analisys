#DEGs Analisys
library(DESeq2)
getwd()

cts = read.delim("nat_bw2-ReadCount.tab", sep="\t", header = T)
head(cts)
dim(cts)

colnames(cts)= c("arvc1","arvc2","arvc3","control1","control2","arvc4","arvc5","arvc6","control3","control4","arvc7","arvc8","control5","arvc9")

Name=colnames(cts)
Time=c(rep("Case", 3), rep("Control", 2), rep("Case",3), rep("Control", 2), rep("Case",2), "Control", "Case")
coldata=data.frame(Name,Time)
head(coldata,14)

log_cts<- log(cts+1, 10) # This is log base 10 + 1 for "0"
table(rowSums(log_cts>log10(3))>=11)

keep.exprs<- rowSums(log_cts>log10(3))>=11
cts_filt<-cts[keep.exprs,]
dim(cts_filt)
head(cts_filt,10)

dds <- DESeqDataSetFromMatrix(countData = round(cts_filt), colData = coldata, design = ~ Time)
dds$Time <- factor(dds$Time, levels = c("Control","Case"))
head(dds)

dds <- DESeq(dds)
res <- results(dds)
head(res)
dim(res)

resOrdered_padj <- res[which(res$padj < 0.0001),]
dim(resOrdered_padj) #436 DEGs

#Find TOP20
resOrdered_padj_fc_neg <- resOrdered_padj[which(resOrdered_padj$log2FoldChange <= -2.5),]
dim(resOrdered_padj_fc_neg)
resOrdered_padj_fc_neg <- resOrdered_padj_fc_neg[order(resOrdered_padj_fc_neg$log2FoldChange),]
head(resOrdered_padj_fc_neg, 12)

resOrdered_padj_fc_pos <- resOrdered_padj[which(resOrdered_padj$log2FoldChange >= 2.5),]
resOrdered_padj_fc_pos <- resOrdered_padj_fc_pos[order(resOrdered_padj_fc_pos$log2FoldChange),]
head(resOrdered_padj_fc_pos, 15)

top20=rbind(resOrdered_padj_fc_neg,resOrdered_padj_fc_pos)
head(top20, 25)
dim(top20)
top20 <- resOrdered_padj[rownames(resOrdered_padj) %in% select_var, ] 
dim(top20)

write.csv(res,"DEGs_arvc.csv", row.names = T)
