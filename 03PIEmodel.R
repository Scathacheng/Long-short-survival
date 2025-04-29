#long-short-survival project
#MacOS 15.4.1
#R 4.4.2
#Jiacheng Dai

library(survival)
library(reshape2)

setwd("Downloads/data-analysis/longshortsurvival/")

#data
#load DEG1 and DEG2 data
load("01data/01ourdata/workingdata.immune.RData")
datExpr1 = datExpr
datMeta1 = datMeta
load("01data/01ourdata/workingdata.IV_immunechemo.RData")
datExpr2 = datExpr
datMeta2 = datMeta

datExpr = cbind(datExpr1, datExpr2)
datMeta = rbind(datMeta1, datMeta2)
save(datExpr,datMeta,file="01data/01ourdata/workingdata.modelinput.RData")

#load another datMeta file
vars = c("gender","age","smoking","Psscore","stage","stage.norm","pathology","therapy","efficacy","efficacy.norm","PFS.norm","PFS.end","Survival")
datMeta1 = datMeta[,vars]





######unicox
load("01data/01ourdata/workingdata.modelinput.RData")

genelist = c("KRT14","FOLR2","SLC31A2","EFCAB14")

dat = t(datExpr)
outTable = lapply(genelist, function(i){
    datMeta$gene = ifelse(dat[,i] > median(as.numeric(dat[,i])), "High", "Low")
    datMeta$gene = factor(datMeta$gene, levels = c("Low", "High"))
    cox = coxph(Surv(PFS.norm, PFS.end) ~ gene, data = datMeta)
    coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    outTable = cbind(id=i,
                z=coxSummary$coefficients[,"z"],
                HR=coxSummary$conf.int[,"exp(coef)"],
                HR.95L=coxSummary$conf.int[,"lower .95"],
                HR.95H=coxSummary$conf.int[,"upper .95"],
                pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
})
outTable = do.call(rbind, outTable)
outTable = as.data.frame(outTable)
outTable = outTable[order(outTable$pvalue),]

#save file
write.csv(outTable, "02result/03model-construct/unicox.longshort.csv")


rs_forest = outTable
rs_forest[,2:6] = sapply(rs_forest[,2:6], function(x){
    dat = as.numeric(x)
    dat = round(dat, 3)
    return(dat)
})
tabletext <- cbind(c("High.vs.Low (median)",rs_forest$id),c("HR (95% CI)",paste0(rs_forest$HR, "(", rs_forest$HR.95L, "-", rs_forest$HR.95H,")")),c("P.Value",rs_forest$pvalue))
tabletext
nrow(tabletext) + 1
library(forestplot)
#files
pdf("02result/03model-construct/unicox.forest2.pdf",onefile = F,width = 7,height = 7)

forestplot(labeltext = tabletext,
           title = "PFS (CellType Model)",
           #grid = 1,
           mean = c(NA, rs_forest$HR), #设置均值
           lower = c(NA, rs_forest$HR.95L), #设置均值的lowlimits限
           upper = c(NA, rs_forest$HR.95H), #设置均值的uplimits限
           is.summary=c(T,rep(F, 100)),
           fn.ci_norm = fpDrawCircleCI,###
           graphwidth = unit(50, 'mm'),
           #align = "l",
           zero = 1, #设置参照值
           boxsize = 0.4, #设置点估计的方形大小
           lineheight = unit(8,'mm'),#设置图形中的行距
           colgap = unit(4,'mm'),#设置图形中的列间距
           #lwd.zero = 1,#设置参考线的粗细
           lwd.ci = 1,#设置区间估计线的粗细
           lwd.xaxis=1,#设置X轴线的粗细
           lty.ci = "solid",
           xlog=F,
           xticks.digits = 1,
           #clip = c(0.1, 16),
           graph.pos = 2,
           txt_gp = fpTxtGp(label = list(gpar(fontsize= 10),
                                         gpar(fontsize = 10
                                         )),
                            ticks = gpar(cex=1),
                            xlab  = gpar(cex = 1)),
           #xticks = c(seq(0.1,1, 0.1), seq(1,16, 2)),
           #col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           #vertices = gpar(col = "blue"),
           shapes_gp = fpShapesGp(
             #default = gpar(col = "black",lty = 2),
             zero = gpar(col = "black",lty = 2,lwd = 1),
             box = gpar(fill = "#3182bd", col = "#3182bd"), # only one parameter
             lines = list( # as many parameters as CI
               gpar(lwd = 1, lty = 2, col = "black")
             )#,
             # grid = gpar(lwd = 0, lty = 2,col = "blue")
           )
)
dev.off()



###calculate effect
rm(list=ls())
load("01data/01ourdata/workingdata.modelinput.RData")
datMeta1 = datMeta
testdata = t(datExpr)
testdata = as.data.frame(testdata)
attach(testdata)
datMeta1$main_TNM = 0.076073733486049*KRT14 + -0.271058244371388*FOLR2 + -0.471524129699845*SLC31A2 + -0.405844167909755*EFCAB14
datMeta1$KRT14_mRNA = KRT14
datMeta1$FOLR2_mRNA = FOLR2
datMeta1$SLC31A2_mRNA = SLC31A2
datMeta1$EFCAB14_mRNA = EFCAB14
datMeta1$PDL1_mRNA = CD274
detach(testdata)

#construct gene model
summary(coxph(Surv(PFS.norm, PFS.end)~main_TNM + PDL1_mRNA, data=datMeta1))$coefficients
summary(coxph(Surv(PFS.norm, PFS.end)~main_TNM + age + smoking + Psscore + PDL1_mRNA, data=datMeta1))$coefficients
datMeta1$CellType_PDL1_model = 0.6170483*datMeta1$main_TNM + -0.0881275*datMeta1$PDL1_mRNA
datMeta1$new_model = 0.38361640*datMeta1$main_TNM + -0.01196873*datMeta1$age + 0.21270505*datMeta1$smoking + 0.99106911*datMeta1$Psscore + -0.15546117*datMeta1$PDL1_mRNA
write.csv(datMeta1, "02result/03model-construct/ImmuneModel.longshort.csv")

#heatmap 
library(scales)
datMeta2 = datMeta1
datMeta2$scale_TNM = rescale(datMeta2$main_TNM,to=c(-2,2))
datMeta2$scale_CA = rescale(datMeta2$main_CA,to=c(-2,2))
datMeta2$scale_PDL1 = rescale(datMeta2$PDL1_mRNA,to=c(-2,2))
datMeta2$scale_model1 = rescale(datMeta2$CellType_model,to=c(-2,2))
datMeta2$scale_model2 = rescale(datMeta2$CellType_PDL1_model,to=c(-2,2))
write.csv(datMeta2, "02result/03model-construct/Heatmapfile.csv")


#ROC1 main figure
library(timeROC)
rm(list=ls())
datMeta1 = read.csv("02result/03model-construct/ImmuneModel.longshort.csv",header=T,row.names=1)

source("00code/ROCplot.R")

model1 = coxph(Surv(PFS.norm, PFS.end) ~ PD.L1, data = datMeta1)
model2 = coxph(Surv(PFS.norm, PFS.end) ~ CellType_PDL1_model, data = datMeta1)

lp1 <- predict(model1,type="lp",newdata = datMeta1)
ROC1 <- timeROC(T=datMeta1$PFS.norm,
                 delta=datMeta1$PFS.end, marker=lp1,
                 cause=1,weighting="marginal",iid=T,
                 times=c(6,12,24))

lp2 <- predict(model2,type="lp",newdata = datMeta1)
ROC2 <- timeROC(T=datMeta1$PFS.norm,
                 delta=datMeta1$PFS.end, marker=lp2,
                 cause=1,weighting="marginal",iid=T,
                 times=c(6,12,24))

#
plotDs=ROC1
legendLab <- c(paste0('6-month: ',sprintf("%0.3f",round(plotDs$AUC[1],3)),
                      ' (',sprintf("%0.3f",round(confint(ROC1)$CI_AUC[1,1]/100,3)),', ',sprintf("%0.3f",round(confint(ROC1)$CI_AUC[1,2]/100,3)),')'),
               paste0('12-month: ',sprintf("%0.3f",round(plotDs$AUC[2],3)),
                      ' (',sprintf("%0.3f",round(confint(ROC1)$CI_AUC[1,1]/100,3)),', ',sprintf("%0.3f",round(confint(ROC1)$CI_AUC[2,2]/100,3)),')'),
               paste0('24-month: ',sprintf("%0.3f",round(plotDs$AUC[3],3)),
                      ' (',sprintf("%0.3f",round(confint(ROC1)$CI_AUC[2,1]/100,3)),', ',sprintf("%0.3f",round(confint(ROC1)$CI_AUC[3,2]/100,3)),')'))
myPlot(plotDs=plotDs,
       colorList=c('darkblue','darkgreen','darkred'),
       smooth=T,
       main="ROC for PFS prediction of long-short survival",
       legendLab=legendLab,
       legendCol=c('darkblue','darkgreen','darkred'),
       legendLty=c(1,1))
