library(patchwork)
# 第一次使用需要安装
wrap_plots(p,nrow=2, guides="collect")

setwd('d:\\A\\work\\extra work')

#富集分析

#标准富集分析
library(DOSE)
options(connectionObserver = NULL)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)

keytypes(org.Hs.eg.db) 
x1 <-c( 'KLF9','UHRF1','CDC25A', 'CCNE1','CDK2','CCNE2','TUBB','TAF11','POLE2', 'PKMYT1','KLHDC1',
        'CDC45', 'ZBTB4', 'UBE2S','CDKN2CNEK2', 'TOMM40', 'TACC3','GPR19','TCEAL5', 'FUS',
        'SIK2' ,'AP1S1SHB', 'HS6ST1','TP73', 'GATA3', 'HOXA10','CD3EAP','SLC20A1','XKR5','SOX4')
x2 <- c('NFKB2', 'NFYA', 'NR3C1', 'MAZ', 'BAX', 'TNIP2', 'DGKZ', 'PIK3R1IL4I1')
x3<- c('KDM4B', 'STAT5B', 'BCL2', 'TRIM59')
x4<- c('TNXB', 'MAFF', 'CTGF', 'KLF4', 'JUN')
x5<- c('BRCA1', 'RAD51')

x<- append(x1,x2)
x<-  append(x,x3)
x<-  append(x,x4)
x<-  append(x,x5)
x<-c('LEF1','POU5F1','GATA3','SOX2'	)
x<-c('SMAD2',	'LEF1',	'GATA3','SOX2')
x<-read.table('d:\\A\\work\\2019\\scrna\\补充的结果分析\\42.txt')
sox2<-read.table('d:\\A\\work\\2019\\scrna\\补充的结果分析\\sox2.txt')
sox2lef1<-read.table('d:\\A\\work\\2019\\scrna\\补充的结果分析\\sox2lef1.txt')
com3<-read.table('d:\\A\\work\\2019\\scrna\\补充的结果分析\\3.txt')
comb4<-read.table('d:\\A\\work\\2019\\scrna\\补充的结果分析\\41.txt')
test = bitr(M2$V1, #数据集
            fromType="SYMBOL", #输入为SYMBOL格式
            toType="ENTREZID",  # 转为ENTERZID格式
            OrgDb="org.Hs.eg.db") #人类 数据库
head(test,2)

data(geneList, package="DOSE") #富集分析的背景基因集
gene <- names(geneList)[abs(geneList) > 2]
gene.df <- bitr(gene, fromType = "ENTREZID", toType = c("ENSEMBL", "SYMBOL"), OrgDb = org.Hs.eg.db)
head(gene.df,2)

ego_ALL <- enrichGO(gene = test$ENTREZID, 
                    universe = names(geneList), #背景基因集
                    OrgDb = org.Hs.eg.db, #没有organism="human"，改为OrgDb=org.Hs.eg.db
                    #keytype = 'ENSEMBL',
                    ont = "BP", #也可以是 CC  BP  MF中的一种
                    pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                    pvalueCutoff = 1, #P值会过滤掉很多，可以全部输出
                    qvalueCutoff = 1,
                    readable = TRUE) #Gene ID 转成gene Symbol 易读
head(ego_ALL,2)

ego_MF <- enrichGO(gene = test$ENTREZID, universe = names(geneList),OrgDb = org.Hs.eg.db,ont = "BP", pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = FALSE)
ego_MF1 <- setReadable(ego_MF, OrgDb = org.Hs.eg.db)
write.csv(as.data.frame(ego_ALL),"M2-enrich.csv",row.names =FALSE)

dotplot(ego_MF,title="EnrichmentGO")#点图，按富集的数从大到小的
##可视化--条形图
barplot(ego_MF, showCategory=20,title="EnrichmentGO")#条状图，按p从小到大排，绘制前20个Term

                                                   #条状图，按p从小到大排，绘制前20个Term

ego3 <- mutate(ego_MF, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
library(tidyverse)
library(DOSE)
ggplot(ego3,aes(richFactor,fct_reorder(Description,richFactor)))+
  geom_segment(aes(xend=0,yend=Description))+
  geom_point(aes(color=p.adjust,size=Count))+
  scale_color_gradientn(colours=c("#f7ca64", "#46bac2", "#7e62a3"),
                        trans = "log10",
                        guide=guide_colorbar(reverse=TRUE, order=1)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_dose(12)+
  labs(x="Rich Factor",
       y="",
       title = "Biological Processes")



multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

M1=read.table("d:\\A\\work\\extra work\\M1.txt" )
M2=read.table("d:\\A\\work\\extra work\\M2.txt" )





