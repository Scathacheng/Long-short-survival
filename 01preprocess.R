#long-short-survival project
#MacOS 15.4.1
#R 4.4.2
#Jiacheng Dai

setwd("Downloads/data-analysis/longshortsurvival/")

###load metadata
datMeta = read.csv("metadata.csv",header=T,row.names=1)
vars = c("gender","age","smoking","Psscore","stage","stage.norm","pathology","therapy","efficacy","efficacy.norm","PFS.norm","Survival")
datMeta = datMeta[,vars]
sapply(datMeta, class)
#肖建林 失访，家属拒绝告知
datMeta = datMeta[-121,] 

###load exp data
data = read.csv("RNA.expr.csv",header=T)

#remove duplicated gene
data1 = data[!duplicated(data$gene_id),] 
rownames(data1) = data1[,1]
data1 = data1[,-1]

expr = data[,-1]
s = unique(data$gene_id[duplicated(data$gene_id)])
data2 = sapply(s, function(i){
    colMeans(expr[data$gene_id==i,])
    })
data2 = t(data2)
test = data1[setdiff(rownames(data1),s),]   
data3 = rbind(data2, test); 
#gene = 25190

#remove empty row
data3 = data3[-1,]

#sample match
data3 = data3[,rownames(datMeta)]
data4 = data3[rowSums(data3 > 0) >= 60,]
# [1] 52477   120

save(data4,datMeta, file="01data/01ourdata/clean_data.RData")

###1
#考虑去掉II期的样本，仅看immune样本 
#25例样本
rm(list=ls())
load("01data/01ourdata/clean_data.RData")
meta=read.csv("01data/01ourdata/metadata.csv",header=T,row.names=1)
meta1 = meta[meta$stage %in% c("IV","IVA","IVB","IIIA","IIIB","IIIC"),]
meta1 = meta1[meta1$therapy %in% c("immune"),]
meta1 = meta1[meta1$Psscore %in% c(0,1),]
meta1$names = rownames(meta1)
datMeta = meta1

data5 = log2(data4+1)
data6 = data5[,rownames(datMeta)]
X = model.matrix(~gender+age+smoking, datMeta)
Y = data6
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X[,2:4]) %*% (as.matrix(beta[2:4,])))
datExpr = data6-t(to_regress)

save(datExpr, datMeta, file="01data/01ourdata/workingdata.immune.RData")

###2
#考虑去掉II和III期的样本，仅看IV期的immune+chemo样本
#71样本
rm(list=ls())
load("01data/01ourdata/clean_data.RData")
meta=read.csv("01data/01ourdata/metadata.csv",header=T,row.names=1)
meta1 = meta[meta$stage %in% c("IV","IVA","IVB"),]
meta1 = meta1[meta1$therapy %in% c("chemo+immune"),]
meta1 = meta1[meta1$Psscore %in% c(0,1),]
meta1$names = rownames(meta1)
datMeta = meta1

data5 = log2(data4+1)
data6 = data5[,rownames(datMeta)]
X = model.matrix(~gender+age+smoking, datMeta)
Y = data6
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X[,2:4]) %*% (as.matrix(beta[2:4,])))
datExpr = data6-t(to_regress)

save(datExpr, datMeta, file="01data/01ourdata/workingdata.IV_immunechemo.RData")

