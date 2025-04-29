#long-short-survival project
#MacOS 15.4.1
#R 4.4.2
#Jiacheng Dai

###DEG analysis
library(nlme)
library(sva)
library(reshape2)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)

rm(list=ls())
load("01data/01ourdata/workingdata.RData")
datExpr = as.matrix(datExpr)
datMeta$Survival = factor(datMeta$Survival,levels=c("short","middle","long"))

meta1 = matrix(NA, nrow=nrow(datExpr), ncol=3) #middle
meta2 = matrix(NA, nrow=nrow(datExpr), ncol=3) #long
for(i in 1:nrow(datExpr)) {
  if(i%%100==0) print(i)
  expr = datExpr[i,]
  tryCatch({
    res = summary(lme(expr~Survival, data = datMeta, random=~1|names))
    meta1[i,] = res$tTable[2,c(1,2,5)]
    meta2[i,] = res$tTable[3,c(1,2,5)]
  }, error=function(e){})
}
meta1 = as.data.frame(meta1)
meta2 = as.data.frame(meta2)
colnames(meta1) = colnames(meta2) = c("beta", "SE", "p")
rownames(meta1) = rownames(meta2) = rownames(datExpr)
meta1$fdr = p.adjust(meta1$p, "fdr")
meta2$fdr = p.adjust(meta2$p, "fdr")
write.csv(meta1, "02result/01DEG/meta1.csv")
write.csv(meta2, "02result/01DEG/meta2.csv")

summary=data.frame(meta1$beta,meta1$p, meta2$beta, meta2$p)
rownames(summary) = rownames(meta1)
colnames(summary) = c("beta.middle","p.middle","beta.long","p.long")
write.csv(summary, file = "02result/01DEG/lme-summary.csv")



###immune cell GSEA
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
library(dplyr)
library(pheatmap)
library(tinyarray)

rm(list=ls())
load("01data/01ourdata/workingdata.RData")
#1
geneset = read.csv("01data/gene_reference.csv",header=T)
geneset = geneset[,1:2]
list = split(as.matrix(geneset)[,1], geneset[,2])
#2
geneset = read.csv("01data/NKcellTangCell2023.csv",header=T)
list = split(as.matrix(geneset)[,1], geneset[,2])
#3
geneset = read.csv("01data/TcellGuoNM2018.csv",header=T)
list = split(as.matrix(geneset)[,1], geneset[,2])

datMeta1 = datMeta[datMeta$Survival %in% c("short","long"),]
datExpr1 = datExpr[,rownames(datMeta1)]
data.gsva = gsva(as.matrix(datExpr), list, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)

data.gsva.1 = scale(data.gsva)
data.gsva.1 = apply(data.gsva, 1, scale)
rownames(data.gsva.1) = colnames(data.gsva)
data.gsva.1 = t(data.gsva.1)

normalization = function(x) {
    return((x-min(x))/(max(x)-min(x)))
}
normdata.gsva.1 = normalization(data.gsva.1)

meta1 = datMeta[datMeta$Survival=="short",]
meta3 = datMeta[datMeta$Survival=="long",]
meta = rbind(meta1, meta3)
annotation_col = data.frame(group = meta$Survival)
rownames(annotation_col) = rownames(meta)
normdata.gsva.1 = normdata.gsva.1[,rownames(meta)]
para = unique(c(seq(0,1,length = 100)))

pheatmap(normdata.gsva.1,
         show_colnames = F,
         cluster_rows = F,cluster_cols = F,
         annotation_col = annotation_col,
         breaks=para,
         cellwidth=2, cellheight=4,
         fontsize=2, #gaps_row = c(12,20),
         filename = 'ssgsea.pdf', 
         width = 10,
         height = 6)

data.gsva = gsva(as.matrix(datExpr1), list, method = "ssgsea", kcdf = "Gaussian", abs.ranking = TRUE)
mygroup = as.character(datMeta$Survival)
pdf("02result/01DEG/immunescore_NKcell.pdf",h=6,w=20)
draw_boxplot(data.gsva, mygroup, color=c("#F6995C","#51829B","#EADFB4"))
dev.off()

