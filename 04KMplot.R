#long-short-survival project
#MacOS 15.4.1
#R 4.4.2
#Jiacheng Dai

###KM plot
library(survminer)
library(survival)
library(cowplot)

datMeta1 = read.csv("02result/03model-construct/ImmuneModel.longshort.csv",header=T,row.names=1)
data = datMeta1

#main_TNM
#main_CA
#PDL1_mRNA
#CellType_model
#CellType_PDL1_model

data$SCore = ifelse(data$main_TNM<mean(data$main_TNM),"Low-risk","High-risk")
data$SCore = as.factor(data$SCore)
data$SCore = relevel(data$SCore, ref="Low-risk")
summary(coxph(Surv(PFS.norm, PFS.end)~as.factor(SCore), data = data))
fit = survfit(Surv(PFS.norm, PFS.end)~SCore, data = data)

a<-ggsurvplot(fit,data = data,title="long-short survival - CellType Model",
              xlim=c(0,60),
              ylab="PFS endpoints",
              xlab="PFS",
              font.x=12,
              font.y=12,
              legend.title="",
              legend.labs=c("Low-risk","High-risk"),
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
                    label = expression(paste('HR = 1.6530 ',', ',  italic(P),' = 0.0306')),colour="black", size = 3.5,parse = TRUE)
#####################################################


#
datMeta1$dose = cut(datMeta1$CellType_PDL1_model, breaks = quantile(datMeta1$CellType_PDL1_model,probs=c(0,0.2,0.4,0.6,0.8,0.9,1),na.rm = T),labels = 1:6)

plotFit <- survfit(coxph(Surv(PFS.norm, PFS.end)~strata(dose),data = datMeta1))
sur<-data.frame(plotFit[["time"]],plotFit[["surv"]])

library(survminer)
library(RColorBrewer)
color <- c("#31789F","#9C2351","#DBA541","#8FBC8F","#00CED1","#7B68EE")
p <- ggsurvplot(plotFit,data = datMeta1,title="Discriminative ability of Model",
                #xlim=c(0,120),
                #pval = "HR = 4.05 , P = 2.55¡Á10-12",
                #pval.coord=c(100,1),
                #ggtheme=theme_survminer(),  
                ylab="6-survival probability",
                xlab="Overall Survival Months",censor=FALSE,
                font.x=12,
                font.y=12,
                legend.title="",
                legend.labs=c("Level 1","Level 2","Level 3","Level 4","Level 5","Level 6"),
                legend="bottom",
                break.x.by = 6,
                #break.y.by = c(0.33,0.50,0.70,0.86,0.93,0.94),
                surv.median.line="v",
                surv.plot.height=0.8,
                palette = color,
                ggtheme = theme(plot.title = element_text(hjust = 0.5),
                                plot.margin = unit(rep(0,4),"cm"),
                                panel.spacing = unit(rep(0,4),"cm"),
                                panel.border=element_rect(color="black",fill=NA),
                                axis.line=element_line(size=0.5),
                                legend.background = element_blank(),
                                legend.text =element_text(size=9),
                                plot.background = element_blank(),
                                panel.background = element_blank(),
                                legend.box.background=element_blank(),
                                #legend.box = TRUE,
                                #legend.box.spacing = unit(rep(0,4),"cm"),
                                legend.key = element_blank(),
                                axis.text.y = element_text(color=c('black',rev(color)[1:2],'black',rev(color)[3:6]),size = 9),
                                axis.text.x = element_text(color="black",size = 9),
                                #axis.title.y = element_text(size = 10),
                                axis.title.x = element_text(margin=unit(c(0.1,0.1,0,0.1),"cm"))))
p$plot<-p$plot+
  ggplot2::scale_y_continuous(expand = c(0.01,0.01),sec.axis = sec_axis(~.,name="24-survival probability",breaks = c(0,0.016,0.204,0.50,0.521,0.723,0.801,0.945),labels = c(0,1.6,20.4,50,52.1,72.3,80.1,94.5)),breaks = c(0,0.079,0.331,0.50,0.715,0.885,0.952,0.986),labels = c(0,7.9,33.1,50,71.5,88.5,95.2,98.6))+
  #ggplot2::scale_y_continuous(expand = c(0.01,0.01),sec.axis = sec_axis(~.,name="5-survival probability",breaks = c(0.23,0.50,0.54,0.68,0.77,0.86),labels = c(0.23,0.50,0.54,0.68,0.77,0.86)))+    
  
  ggplot2::scale_x_continuous(expand = c(0.01,0.01))+
  #ggplot2::geom_abline(intercept=c(36,0.863),linetype=2,size=0.5)+
  #ggplot2::geom_vline(xintercept = c(36,60),linetype=2,size=0.5)+
  ggplot2::geom_hline(yintercept=c(0.5),linetype=2,size=0.5)
p

