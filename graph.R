library(igraph)
cancermodule=read.csv("edge509.csv")
#edgeweight=cancermodule[,4]
network <- graph.data.frame(cancermodule[,c(2,3)],directed = F)
fc <- cluster_fast_greedy(network)
a=membership(fc)
b=as.data.frame(t(sapply(a, "[", i = 1:max(sapply(a, length)))))

sizes(fc)
plot.membership<-function(graph,membership,main=""){
  V(graph)$member<-membership
  mem.col<-rainbow(length(unique(membership)),alpha=0.6)
  V(graph)$width=100
  V(graph)$color<-mem.col[membership]
  V(graph)$size = 6
  V(graph)$label=fc$names
  plot(graph,edge.width=E(graph)$weight,vertex.color=V(graph)$color,main=main,layout=layout.lgl)
  
}
svg(filename="2.svg",width = 40,height = 40)
plot.membership(network,a)
dev.off()
plot.membership(network,a)
library(pheatmap)
#创建数据集test测试矩阵
test = matrix(rnorm(200), 20, 10)
test[1:10, seq(1, 10, 2)] = test[1:10, seq(1, 10, 2)] + 3
test[11:20, seq(2, 10, 2)] = test[11:20, seq(2, 10, 2)] + 2
test[15:20, seq(2, 10, 2)] = test[15:20, seq(2, 10, 2)] + 4
colnames(test) = paste("Test", 1:10, sep = "")
rownames(test) = paste("Gene", 1:20, sep = "")


# 用pheatmap函数画热图
pheatmap(test)

name <- "HiSeqV2_BRCA_outcome.txt"
data <- as.matrix(read.table(name,sep = "\t",header=T))

dim(data)
#17924*1218

normal=cbind()
cancer=cbind()
for (i in 1:1218) {
  if(data[1,i]==1) {
    normal=cbind(normal,data[,i])
  }
  else{
    cancer=cbind(cancer,data[,i])
  }
  
  
}
dim(normal)
dim(cancer)
#######################################
library(ggpubr)
library(ggplot2)
st=compare_means(expression~Label, dat2, method = "t.test", paired = FALSE,group.by = 'gene')
st["p"]
storder=st[order(st['p.format']),]
pvalueordergene=storder["gene"]
module1node <-c('SP1','USF1','RBL1','MIF','MAPK1','GABPA','EP300','RELB','CREB1','BAX','HRAS','NFKBIB','AR','NFKB2','RB1',"Label")
#dat=data[,module1node]
dat=data[,append(as.vector(unlist(pvalueordergene[1])),"Label")]#这里是想通过改变gene的顺序使箱线图排序
dat=data.frame(dat)
normal=t(normal)
cancer=t(cancer)
heatmapnode<-c('SP1','USF1','RBL1','MIF','MAPK1','GABPA','EP300','RELB','CREB1','BAX','HRAS','NFKBIB','AR','NFKB2','RB1')
normaldata=normal[,heatmapnode]
cancerdata=cancer[,heatmapnode]
pheatmap(cor(normaldata),cluster_row = FALSE, cluster_col = FALSE)
pheatmap(cor(cancerdata),cluster_row = FALSE, cluster_col = FALSE)

data=t(data)
library(ggplot2)
library(tidyr)
library(dplyr)
dat2 = gather(dat,key = "gene",value = "expression",-Label) ###画箱线图的数据
dat2["sort"]=1:30



################箱线图

boxfigure=ggplot(data = dat2,aes(x = as.character(Label),y = expression,color = as.character(Label)))+
  stat_compare_means(aes(group = Label), label = "p.format",label.y=5.5)+
  geom_boxplot()+
  scale_y_continuous(name = "Expression")+
  scale_x_discrete(name = "Gene") +guides(color = guide_legend(title = 'Condition'))+
  scale_color_discrete(name = "Conditions", labels = c("Disease", "Normal"))+
  
  theme_bw()
  
boxfigure+facet_wrap(~gene,nrow = 3)+
  scale_color_manual(values=c("#8B2323", "#66CDAA"), 
                  name="Condition",
                  breaks=c("0", "1"),
                  labels=c("Disease", "Normal"))
t.test(dat2["expression"]~dat2["Label"])
#差不多就到这里



p <- ggboxplot(dat2, x = "gene", y = "expression",
               color = "Label", palette = "jco", 
               add = "jitter") # palette可以按照期刊选择相应的配色，如"npg"等
dat2['Label']=as.character(dat2['Label'])
p + stat_compare_means(aes(group = Label), label = "p.signif")

data("ToothGrowth")
df1 <- ToothGrowth
head(df1)
p <- ggboxplot(dat2, x="gene", y="expression", color = "Label", 
               #palette = c("#00AFBB", "#E7B800", "#FC4E07"), 
               add = "jitter", shape="Label")#增加了jitter点，点shape由dose映射

# 增加不同组间的p-value值，可以自定义需要标注的组间比较
my_comparisons <- list(c("0.5", "1"), c("1", "2"), c("0.5", "2"))
p+stat_compare_means(comparisons = my_comparisons)+ #不同组间的比较 
  stat_compare_means(label.y = 50)

boxfigure+ 

geom_boxplot(alpha=0.7) +

  ggtitle("Boxplot of hub gene") +

# Color by group (dose)
e + geom_boxplot(aes(color = dose))+
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))
# Change fill color by group (dose)
e + geom_boxplot(aes(fill = dose)) +
  scale_fill_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))

# Choose which items to display: group "0.5" and "2"
e + geom_boxplot() + 
  scale_x_discrete(limits=c("0.5", "2"))
# Change the default order of items
e + geom_boxplot() +
  scale_x_discrete(limits=c("2", "0.5", "1"))

e2 <- e + geom_boxplot(
  aes(fill = supp),
  position = position_dodge(0.9) 
) +
  scale_fill_manual(values = c("#999999", "#E69F00"))
e2

p=list()
#list 列表将储存所有的循环绘图，以实现后面的多图组合
library(ggplot2)
for (i in 1:(ncol(dat)-1)){
  #之所以减一，是因为最后一列是组类别
  p[[i]]=ggplot(data=dat,aes_string(x="Label",y=colnames(dat)[i]))+
    geom_boxplot(aes(color=as.character(Label)))+ 
    geom_jitter(aes(color=as.character(Label)))
  } 


library(patchwork)
# 第一次使用需要安装
wrap_plots(p,nrow=2, guides="collect")



#富集分析

#标准富集分析
library(DOSE)
library(org.Hs.eg.db)
library(topGO)
library(clusterProfiler)

keytypes(org.Hs.eg.db) 
x <- c('PGR', 'MYB', 'CTNNB1', 'UHRF1', 'IGF1R')

test = bitr(x, #数据集
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
                    ont = "ALL", #也可以是 CC  BP  MF中的一种
                    pAdjustMethod = "BH", #矫正方式 holm”, “hochberg”, “hommel”, “bonferroni”, “BH”, “BY”, “fdr”, “none”中的一种
                    pvalueCutoff = 1, #P值会过滤掉很多，可以全部输出
                    qvalueCutoff = 1,
                    readable = TRUE) #Gene ID 转成gene Symbol 易读
head(ego_ALL,2)

ego_MF <- enrichGO(gene = test$ENTREZID, universe = names(geneList),OrgDb = org.Hs.eg.db,ont = "MF", pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = FALSE)
ego_MF1 <- setReadable(ego_MF, OrgDb = org.Hs.eg.db)
write.csv(as.data.frame(ego_ALL),"ALL-enrich3.csv",row.names =FALSE)
dotplot(ego_MF,title="Module3 EnrichmentGO")#点图，按富集的数从大到小的
##可视化--条形图
barplot(ego_MF, showCategory=20,title="Module3 EnrichmentGO")#条状图，按p从小到大排，绘制前20个Term
