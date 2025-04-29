#long-short-survival project
#MacOS 15.4.1
#R 4.4.2
#Jiacheng Dai

library(MCPcounter)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(GSVA)
library(openxlsx)
library(tibble)
library(tidyverse)

setwd("Downloads/data-analysis/longshortsurvival/")
rm(list=ls())
# load data
load("01data/01ourdata/workingdata.modelinput.RData")
datMeta1 = datMeta[datMeta$Survival %in% c("short","long"),]
datExpr1 = datExpr[,rownames(datMeta1)]
# load reference
datref = read.xlsx("01data/06figure/ao.pathway.xlsx", sheet = 5)
datref
datalist = lapply(datref, function(x) x[!is.na(x)])
datalist #a list

# GSVA analysis
expr = as.matrix(datExpr1)
gsva_par = ssgseaParam(expr, datalist)
scores = gsva(gsva_par)
save(scores, file="01data/06figure/gs.exp.ao.rda")


mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
identical(colnames(scores), datMeta1$names)
scores1 = t(scores)
scores1 = as.data.frame(scores1)
identical(rownames(scores1), datMeta1$names)
scores1$group = datMeta1$Survival


gene = names(scores1)[-ncol(scores1)]
gene
pval = sapply(gene, function(x){
  dat = scores1[,c(x, "group")]
  names(dat)[1] = "gene"
  test = wilcox.test(gene ~ group, dat)
  return(test$p.value)
})
pval.sig = ifelse(pval < 0.001, "***",
                  ifelse(pval < 0.01, "**",
                         ifelse(pval < 0.05, "*", ""
                         )
                  )
)

##heatmap## 
library(ggplot2)
library(ComplexHeatmap)
library(GenAnalysis)

scores1 = scores1[,-ncol(scores1)]
scores1 = t(scores1)

scores1.heatmap = scores1
range(scores1.heatmap)

bk1 <- seq(min(scores1.heatmap),0,by=0.01)
bk2 = seq(0,max(scores1.heatmap),by=0.01)
rownames(scores1.heatmap) = ifelse(pval < 0.05, 
                               paste0(rownames(scores1.heatmap), "(",pval.sig,")"),
                               rownames(scores1)
)
datMeta1[which(datMeta1$stage=="IVB"),]$stage = "IV"
datMeta1[which(datMeta1$stage=="IVA"),]$stage = "IV"
datMeta1[which(datMeta1$stage=="IIIA"),]$stage = "III"
datMeta1[which(datMeta1$stage=="IIIB"),]$stage = "III"
datMeta1[which(datMeta1$stage=="IIIC"),]$stage = "III"  

datMeta1 = datMeta1[order(datMeta1$Survival),]

top.anno = HeatmapAnnotation(
  Gender = datMeta1$gender,
  Age = datMeta1$age,
  Smoking = datMeta1$smoking,
  PS_score = datMeta1$Psscore,
  Stage = datMeta1$stage,
  Pathological_type = datMeta1$pathology,
  BOR = datMeta1$efficacy,
  Therapy = datMeta1$therapy,
  Group = datMeta1$Survival,
  annotation_name_side = "left",
  simple_anno_size = unit(1.2, "cm"),
  gap = unit(0.2, "cm"),
  annotation_name_gp = gpar(fontsize = 30,fontface = "bold"),
  border = F,
  col = list(Therapy=c("immune"="#67ADB7","chemo+immune"="#E4A6BD"),
                        BOR=c("PD"="#ECA8A9","PR"="#74AED4","SD"="#D3E2B7"),
                        Stage=c("III"="#78A040","IV"="#F3A17C"),
                        PS_score = c("0"="#CFAFD4","1"="#F7C97E"),
                        Gender = c("M"="#F5E1D8","F"="#AFACB7"),
                        Pathological_type = c("adeno"="#67ADB7","squamous"="#F5E1D8","others"="#E4A6BD"),
                        Age = circlize::colorRamp2(c(min(datMeta1$age),
                                 max(datMeta1$age)),
                               c("white", "darkred")
                        ),
                        Smoking = c("1"="#5EAAA6","2"="#D0E8E0"),
                        Group =  c("short" = "#80b1d3", "long" = "#bc80bd")
                        ),
  annotation_legend_param = list(labels_gp = gpar(fontsize = 12),
                                 title_gp = gpar(fontsize = 12),
                                 #legend_gp = gpar(fill = 1:6),
                                 row_gap = unit(2, "cm"),
                                 title_gap = unit(5, "cm"),
                                 grid_width = unit(1, "cm"),
                                 grid_height = unit(1, "cm")
  )
)

scores1.heatmap = scores1.heatmap[,datMeta1$names]
identical(datMeta1$names, colnames(scores1.heatmap))
ht_list = Heatmap(scores1.heatmap, 
	  cluster_rows = T,
    cluster_columns = T,
    show_column_names = F,
    col = c(colorRampPalette(colors = c("blue","white"))(length(bk1)),
                colorRampPalette(colors = c("white","red"))(length(bk2))
        ),
    show_row_dend = F,
    show_column_dend = F,
    row_names_max_width = max_text_width(
          rownames(scores1.heatmap), 
          gp = gpar(fontsize = 12)
        ),
    column_split = datMeta1$Survival,
	  row_names_side = "left",
    row_title = NULL,
    column_title = NULL,
    column_gap = unit(8, "mm"),
    row_names_gp = gpar(fontsize = 12),
	  heatmap_legend_param = list(
		  at = c(-1, 0, 1),
      title_gp = gpar(fontsize = 12),
      labels_gp = gpar(fontsize = 12),
      title = "Enrichment.score",
		  labels = c(-0.5,0,0.5),
		  legend_height=unit(4,"cm")
    ),
	  top_annotation = top.anno
)


pdf("./02result/06figure/ao.immune.heatmap.pdf",width = 15,height = 20)
draw(ht_list, heatmap_legend_side = "right",
     legend_gap = unit(2, "cm"), legend_grouping = "original")
dev.off()


##boxplot##
library(tibble)
library(tidyverse)

setwd("Downloads/data-analysis/longshortsurvival/")
rm(list=ls())

load("01data/01ourdata/workingdata.modelinput.RData")
datMeta1 = datMeta[datMeta$Survival %in% c("short","long"),]
datExpr1 = datExpr[,rownames(datMeta1)]
load("01data/06figure/gs.exp.ao.rda")
ls()
#[1] "datExpr"  "datExpr1" "datMeta"  "datMeta1" "scores"  

#1
mypalette <- colorRampPalette(brewer.pal(8,"Set1"))
identical(colnames(scores), datMeta1$names)
scores1 = scores[c("NonResponse_Ivy","CRMA","B_cell_Budczies"),]

scores1 = t(scores1)
scores1 = as.data.frame(scores1)
identical(rownames(scores1), datMeta1$names)
scores1$Survival = datMeta1$Survival
table(scores1$Survival)

gene = c("NonResponse_Ivy", "CRMA", "B_cell_Budczies")
names(gene) = gene
y.lim = lapply(gene, function(x){
  dat = scores1[,x]
  y.min = min(dat)
  y.max = max(dat)
  return(c(y.min,y.max))
})
y.lim = unlist(y.lim)
blank_data <- data.frame(gene = rep(gene, each = 2),
                         x = 1, y = y.lim)

#2
library(reshape2)
library(ggplot2)
library(ggpubr)
library(rstatix)

scores2 = melt(scores1,variable.name = "gene",id.vars = "Survival")

stat.test <- scores2 %>%
  group_by(gene) %>%
  wilcox_test(value ~ Survival)
stat.test
stat.test <- stat.test %>% add_y_position()

#3
library(ggbeeswarm)
library(RColorBrewer)

p = ggplot(scores2, aes(Survival,value)) + 
  geom_beeswarm(aes(color = Survival),cex = 4,size = 2) + 
  geom_blank(data = blank_data, aes(x = x, y = y)) +
  facet_wrap(.~gene,
             scales = "free_y",
             nrow = 3
  ) +
  stat_pvalue_manual(stat.test, label = "p", tip.length = 0) +
  scale_color_manual(values = c("short" ="#80b1d3", 
                    "long" = "#bc80bd")) +
  #scale_color_manual(values = c(brewer.pal(7, "Set2")[c(1,2,5,7)])) + 
  #ggsci::scale_color_jama()+
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1)))+ 
  labs(y = "GSVA Score")+
  stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median, 
               geom = 'crossbar', width = 0.4, size = 0.2, color = 'black') +
  stat_summary(fun.data = function(x) median_hilow(x, 0.5), 
               geom = 'errorbar', width = 0.25, color = 'black')+
  theme(#panel.grid = element_blank(), 
    panel.background = element_blank(), legend.position = "bottom",
    panel.border = element_rect(color = "black",fill = NA), 
    axis.line = element_line(color = 'black'),
    # axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1,
    #                            color = 'black'),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.length.x.bottom = unit(0,"mm"),
    #axis.line.x.top = element_line(color = 'black'),
    axis.text.y = element_text(size = 15,
      color = 'black'),
    axis.title.y = element_text(size = 20),
    plot.margin = unit(c(10,10,10,10), units = "mm"))
p
ggsave(p, file = "./02result/06figure/ao.boxplot.pdf",width = 4,height = 8)


##Univariate cox regression for CRMA##
library(survival)
library(reshape2)

rm(list=ls())
load("01data/01ourdata/workingdata.modelinput.RData")
genelist = c("MAGEA10","MAGEA12","MAGEA2","MAGEA4","MAGEA8","MAGEA11","MAGEA1","MAGEA5")

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

pdf("./02result/06figure/unicox.CRMA.pdf",onefile = F,width = 7,height = 7)
forestplot(labeltext = tabletext,
           mean = c(NA, rs_forest$HR), #设置均值
           lower = c(NA, rs_forest$HR.95L), #设置均值的lowlimits限
           upper = c(NA, rs_forest$HR.95H), #设置均值的uplimits限
           is.summary=c(T,rep(F, 100)),
           fn.ci_norm = fpDrawCircleCI,###
           graphwidth = unit(50, 'mm'),
           zero = 1, #设置参照值
           boxsize = 0.4, #设置点估计的方形大小
           lineheight = unit(8,'mm'),#设置图形中的行距
           colgap = unit(4,'mm'),#设置图形中的列间距
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


##survival KM plot##
library(survminer)
library(survival)
library(cowplot)

rm(list=ls())
load("01data/01ourdata/workingdata.modelinput.RData")
identical(datMeta$names, colnames(datExpr))

datMeta$gene1 = t(datExpr["MAGEA12",])
datMeta$gene2 = t(datExpr["MAGEA2",])

datMeta$score2 = ifelse(datMeta$gene2<mean(datMeta$gene2),"Low-exp","High-exp")
datMeta$score12 = ifelse(datMeta$gene1<mean(datMeta$gene1),"Low-exp","High-exp")


#MAGEA2
datMeta$score2 = as.factor(datMeta$score2)
datMeta$score2 = relevel(datMeta$score2, ref="Low-exp")
summary(coxph(Surv(PFS.norm, PFS.end)~as.factor(score2), data = datMeta))
fit2 = survfit(Surv(PFS.norm, PFS.end)~score2, data = datMeta)

a<-ggsurvplot(fit2,data = datMeta,title="MAGEA2",
              xlim=c(0,60),
              ylab="PFS endpoints",
              xlab="PFS",
              font.x=12,
              font.y=12,
              legend.title="",
              legend.labs=c("Low-exp","High-exp"),
              legend=c(0.8,0.8),
              risk.table = TRUE,
              break.x.by = 3,
              break.y.by = 0.2,
              palette = c("#31789F","#9C2351"),
              # risk.table.title=NULL,
              #risk.table.pos="out",
              #tables.height=0.15,
              ggtheme = theme(#text = element_text(family = "Arial"),
                plot.margin = unit(rep(0,4),"cm"),
                panel.spacing = unit(rep(0,4),"cm"),
                panel.border=element_blank(),
                axis.line=element_line(size=0.5),
                legend.background = element_blank(),
                legend.text =element_text(size=11),
                plot.background = element_blank(),
                panel.background = element_blank(),
                legend.box.background=element_blank(),
                legend.key = element_blank(),
                axis.text.y = element_text(color="black"),
                axis.text.x = element_text(color="black"),
                #axis.title.y = element_text(margin=unit(c(0.1,0.1,0,0.1),"cm"),size = 10),
                axis.title.x = element_text(margin=unit(c(0.1,0.1,0,0.1),"cm"))),
              tables.theme  = theme(plot.margin = unit(rep(0,4),"cm"),
                                    panel.spacing = unit(rep(0,4),"cm"),
                                    title = element_blank(),
                                    panel.border=element_blank(),
                                    legend.key = element_blank(),
                                    axis.line=element_blank(),
                                    axis.ticks = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.text.x = element_blank(),
                                    legend.background = element_blank(),
                                    plot.background = element_blank(),
                                    panel.background = element_blank(),
                                    # axis.text.y = element_text(color="black"),
                                    legend.box.background=element_blank(),
                                    axis.title.x = element_blank()
              ),
              fontsize=4,
              tables.y.text=T,
              tables.y.text.col=TRUE
)
a$plot<-a$plot+
  ggplot2::annotate("text", x =24, y = 1,
                    label = expression(paste('HR = 1.6697 ',', ',  italic(P),' = 0.0245')),colour="black", size = 3.5,parse = TRUE)

# MAGEA12
datMeta$score12 = as.factor(datMeta$score12)
datMeta$score12 = relevel(datMeta$score12, ref="Low-exp")
summary(coxph(Surv(PFS.norm, PFS.end)~as.factor(score12), data = datMeta))
fit12 = survfit(Surv(PFS.norm, PFS.end)~score12, data = datMeta)

a<-ggsurvplot(fit12,data = datMeta,title="MAGEA12",
              xlim=c(0,60),
              ylab="PFS endpoints",
              xlab="PFS",
              font.x=12,
              font.y=12,
              legend.title="",
              legend.labs=c("Low-exp","High-exp"),
              legend=c(0.8,0.8),
              risk.table = TRUE,
              break.x.by = 3,
              break.y.by = 0.2,
              palette = c("#31789F","#9C2351"),
              # risk.table.title=NULL,
              #risk.table.pos="out",
              #tables.height=0.15,
              ggtheme = theme(#text = element_text(family = "Arial"),
                plot.margin = unit(rep(0,4),"cm"),
                panel.spacing = unit(rep(0,4),"cm"),
                panel.border=element_blank(),
                axis.line=element_line(size=0.5),
                legend.background = element_blank(),
                legend.text =element_text(size=11),
                plot.background = element_blank(),
                panel.background = element_blank(),
                legend.box.background=element_blank(),
                legend.key = element_blank(),
                axis.text.y = element_text(color="black"),
                axis.text.x = element_text(color="black"),
                #axis.title.y = element_text(margin=unit(c(0.1,0.1,0,0.1),"cm"),size = 10),
                axis.title.x = element_text(margin=unit(c(0.1,0.1,0,0.1),"cm"))),
              tables.theme  = theme(plot.margin = unit(rep(0,4),"cm"),
                                    panel.spacing = unit(rep(0,4),"cm"),
                                    title = element_blank(),
                                    panel.border=element_blank(),
                                    legend.key = element_blank(),
                                    axis.line=element_blank(),
                                    axis.ticks = element_blank(),
                                    axis.text.y = element_blank(),
                                    axis.text.x = element_blank(),
                                    legend.background = element_blank(),
                                    plot.background = element_blank(),
                                    panel.background = element_blank(),
                                    # axis.text.y = element_text(color="black"),
                                    legend.box.background=element_blank(),
                                    axis.title.x = element_blank()
              ),
              fontsize=4,
              tables.y.text=T,
              tables.y.text.col=TRUE
)
a$plot<-a$plot+
  ggplot2::annotate("text", x =24, y = 1,
                    label = expression(paste('HR = 1.6660 ',', ',  italic(P),' = 0.0259')),colour="black", size = 3.5,parse = TRUE)
