rm(list = ls())

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


regnetwork <- as.matrix(read.table("regnetwork.txt",sep = "\t",header=F))
dim(regnetwork)
trrust<- as.matrix(read.table("trrust_rawdata.human.tsv",sep = "\t",header=F))

#从这里开始更改背景网络了
relall <- as.matrix(read.csv("relall+1.csv",header=T))
backgroundnetwork=relall[,2:3]
dim(backgroundnetwork)
backgroundnetwork=unique(backgroundnetwork)
dim(backgroundnetwork)
node1=unique(backgroundnetwork[,1])
node2=unique(backgroundnetwork[,2])
node=unique(c(node1,node2)) #20404  #118


nodewithdata=rownames(data) #17924
nodeinnetwork=list()
for(i in nodewithdata) {
  print(i)
  if(i %in% node){
    #print(i)
    nodeinnetwork=append(nodeinnetwork,i)
  }
}

nodeinnetwork=unique(nodeinnetwork)#15985
print(length(nodeinnetwork))#118


###算的太慢了
backgroundnetwork_use=rbind()
for(i in 1:dim(backgroundnetwork)[1]){
    if(backgroundnetwork[i,1] %in% nodeinnetwork){
      if(backgroundnetwork[i,2] %in% nodeinnetwork){
        backgroundnetwork_use=rbind(backgroundnetwork_use,backgroundnetwork[i,])
      }
    }
    
  }
dim(backgroundnetwork_use)    #127938 2

###########这是对图形的处理，去掉selfloop
library(igraph)

#reN <- graph.data.frame(backgroundnetwork_use,directed = F)
reN <- graph.data.frame(backgroundnetwork_use,directed = F)
reN_sim <- simplify(reN,remove.loops=T,remove.multiple = T)#去除selfloop
#plot(reN_sim)
ed_reN_sim <- get.edgelist(reN_sim) #获取边的信
dim(ed_reN_sim)    #122983 2









####  PCC matrix calculation and visulation  建立PCC矩阵  这里用的是exp，就是第一个读入的表达数据

normalexp=normal[-1,]
cancerexp=cancer[-1,]


normalexp=cancer[-1,]####这里直接改了



ex_v <- apply(as.matrix(normalexp),2,
              function(x)as.numeric(x)) 
dim(ex_v)
colnames(ex_v) <-colnames(normalexp)
rownames(ex_v) <- rownames(normalexp) #numeric之后重新赋行列名
ex_v=t(ex_v)
dim(ex_v) #121 17923 #1097 17923







#直接计算网络里的边的PCC,直接判断
matPCC=rbind()
for(i in 1:dim(ed_reN_sim)[1]){
  print(i)
  PCC=cor(ex_v[,ed_reN_sim[i,1]],ex_v[,ed_reN_sim[i,2]],method = "pearson")
  if(is.na(PCC)=="FALSE"){
    if(abs(PCC)>0.4){#这里设定的值
      one=c(ed_reN_sim[i,],PCC)
      matPCC=rbind(matPCC,one)
      }
    }
  
}

dim(matPCC)#0.4 26069 3  原来122983 0.5 14634  #2957 3



###### 0-MI  这里使用的是离散后的表达

library(infotheo)#他的包实现了基于几个熵估计的信息理论的各种度量
#ed=matPCC[,1:2]
ed_2=matPCC[,1:2]

#write.csv(matPCC,"normalmatPCC.csv")
#write.csv(matPCC,"cancermatPCC.csv")


nbins <- sqrt(NROW(ex_v))#大概是存储方格个数
dis_data <- discretize(ex_v,disc="equalfreq",nbins) 
colnames(dis_data) <- colnames(ex_v)

MI_0 <- rbind()
for(i in 1:dim(ed_2)[1])
{
  #i <- 1
  loc1 <- ed_2[i,1]
  loc2 <- ed_2[i,2]
  MI <- mutinformation(dis_data[,loc1],dis_data[,loc2],method="emp")
  MI_0 <- rbind(MI_0,MI)
}
median(MI_0)
mean(MI_0)
####### plot hist MI

hist(MI_0,main='MI_0_distribution',breaks =20)

MI_0_thre <- 0.55
MI_0_left <- rbind()
for(i in 1:length(MI_0))
{
  if(MI_0[i] > MI_0_thre)
  {
    MI_0_left <- rbind(MI_0_left,c(ed_2[i,],MI_0[i]))
  }
}
dim(MI_0_left)


MI_0_net <- graph.data.frame(MI_0_left[,c(1,2)],directed = F)
pdf('./MI_0_net.pdf')
plot(MI_0_net, vertex.color="purple",vertex.size=8,
     label.font=2,label.cex=2,label.color='black',main='MI_0_net')
dev.off()

##### Calculate 1-MI if exists
library(pracma)
MI_1 <- rbind()
MI_2 <- rbind()
MI_3 <- rbind()
MI_re <- rbind()
Node_0 <- unlist(vertex_attr(MI_0_net))
for(i in 1:dim(MI_0_left)[1])
{
  #i <-1
  loc1 <- MI_0_left[i,1]
  loc2 <- MI_0_left[i,2]
  ### if they have shared one-order neighbor
  ne1 <- setdiff(Node_0[unlist(ego(MI_0_net, order = 1, MI_0_left[i,1]))],MI_0_left[i,1])
  ne2 <- setdiff(Node_0[unlist(ego(MI_0_net, order = 1, MI_0_left[i,2]))],MI_0_left[i,1])
  nn <- unique(setdiff(intersect(ne1,ne2),c(MI_0_left[i,1],MI_0_left[i,2])))
  if(isempty(nn)==F)
  {
    # 1-order
    for(j in 1:length(nn))
    {
      loc3 <- nn[j]
      con <- dis_data[,loc3]
      mi <- condinformation(dis_data[,loc1],dis_data[,loc2],con,method="emp")
      MI_1 <- rbind(MI_1,c(MI_0_left[i,c(1,2)],nn[j],1,mi))
    }
   
    
  }
  else{
    #MI_1 <- rbind(MI_1,c(MI_0_left[i,c(1,2)],0,0)) 
    MI_re <-rbind( MI_re,c(MI_0_left[i,c(1,2)],0,'left'))
  }
}
colnames(MI_re) <- c('node1','node2','MI_order','states')
colnames(MI_1) <- c('node1','node2','neighbor1','MI_order','states')



#### select the maximum MI and CMI
thre <- 0.5
#### 1-order CMI
u1 <- unique(MI_1[,1])
MI_1_left_max <- rbind()
dim(MI_1)
for(i in 1:length(u1))
{
  loc <- which(MI_1[,1] %in% u1[i])
  k1 <- unique(MI_1[loc,2])
  can <- MI_1[loc,2]
  for(j in 1:length(k1))
  {
    z1 <- apply(as.matrix(MI_1[loc[which(can %in% k1[j])],dim(MI_1)[2]]),
                2,function(x) as.numeric(x))
    if(max(z1)> thre)
    {
      MI_1_left_max <- rbind(MI_1_left_max,c(u1[i],k1[j],max(z1)))
    }
  }
}
dim(MI_1_left_max)


MI_1_net <- graph.data.frame(MI_1_left_max[,c(1,2)],directed = F)
for(i in 1:dim(MI_1_left_max)[1])
{
  #i <-1
  loc1 <- MI_1_left_max[i,1]
  loc2 <- MI_1_left_max[i,2]
  ### if they have shared one-order neighbor
  ne1 <- setdiff(Node_0[unlist(ego(MI_1_net, order = 1, MI_1_left_max[i,1]))],MI_1_left_max[i,1])
  ne2 <- setdiff(Node_0[unlist(ego(MI_1_net, order = 1, MI_1_left_max[i,2]))],MI_1_left_max[i,1])
  nn <- unique(setdiff(intersect(ne1,ne2),c(MI_1_left_max[i,1],MI_1_left_max[i,2])))
  if(isempty(nn)==F)
  {
   
    # 2-order
    if(length(nn)>1)
    {
      for(k in 1:(length(nn)-1))
      {
        loc3 <- nn[k]
        for(b in (k+1):length(nn))
        {
          loc4 <- nn[b]
          con <- c(dis_data[,loc3],dis_data[,loc4])
          mi <- condinformation(dis_data[,loc1],dis_data[,loc2],con,method="emp")
          MI_2 <- rbind(MI_2,c(MI_0_left[i,c(1,2)],nn[k],nn[b],2,mi))
        }
      }
    }
  }
}

colnames(MI_2) <- c('node1','node2','neighbor1','neighbor2','MI_order','states')
dim(MI_2)
#### 2-order CMI
thre=0.5
u2 <- unique(MI_2[,1])
MI_2_left_max <- rbind()
for(i in 1:length(u2))
{
  loc <- which(MI_2[,1] %in% u2[i])
  k1 <- unique(MI_2[loc,2])
  can <- MI_2[loc,2]
  for(j in 1:length(k1))
  {
    z1 <- apply(as.matrix(MI_2[loc[which(can %in% k1[j])],dim(MI_2)[2]]),
                2,function(x) as.numeric(x))
    if(max(z1) > thre)
    {
      MI_2_left_max <- rbind(MI_2_left_max,c(u2[i],k1[j],max(z1)))
    }
  }
}
dim(MI_2_left_max)

MI_2_net <- graph.data.frame(MI_2_left_max[,c(1,2)],directed = F)
for(i in 1:dim(MI_2_left_max)[1])
{
  #i <-1
  loc1 <- MI_2_left_max[i,1]
  loc2 <- MI_2_left_max[i,2]
  ### if they have shared one-order neighbor
  ne1 <- setdiff(Node_0[unlist(ego(MI_2_net, order = 1, MI_1_left_max[i,1]))],MI_2_left_max[i,1])
  ne2 <- setdiff(Node_0[unlist(ego(MI_2_net, order = 1, MI_1_left_max[i,2]))],MI_2_left_max[i,1])
  nn <- unique(setdiff(intersect(ne1,ne2),c(MI_2_left_max[i,1],MI_2_left_max[i,2])))
  if(isempty(nn)==F)
  {
    #3 order
    if(length(nn)>2)
    {
      for(k in 1:(length(nn)-2))
      {
        loc3 <- nn[k]
        for(b in (k+1):(length(nn)-1))
        {
          loc4 <-nn[b]
          for(h in (b+1):length(nn))
          {
            loc5 <- nn[h]
            con <- c(dis_data[,loc3],dis_data[,loc4],dis_data[,loc5])
            mi <- condinformation(dis_data[,loc1],dis_data[,loc2],con,method="emp")
            MI_3 <- rbind(MI_3,c(MI_2_left_max[i,c(1,2)],nn[k],nn[b],nn[h],3,mi))
          }
        }
      }
    }
  }
}
colnames(MI_3) <- c('node1','node2','neighbor1','neighbor2','neighbor3','MI_order','states')
#### 3-order CMI
#max(apply(as.matrix(MI_3[,7]),2,function(x) as.numeric(x)))
u3 <- unique(MI_3[,1])
MI_3_left_max <- rbind()
dim(MI_3)
for(i in 1:length(u3))
{
  #i<-3
  loc <- which(MI_3[,1] %in% u3[i])
  k1 <- unique(MI_3[loc,2])
  can <- MI_3[loc,2]
  for(j in 1:length(k1))
  {
    #j <-4
    z1 <- apply(as.matrix(MI_3[loc[which(can %in% k1[j])],dim(MI_3)[2]]),
                2,function(x) as.numeric(x))
    if(max(z1) > thre)
    {
      MI_3_left_max <- rbind(MI_3_left_max,c(u2[i],k1[j],max(z1)))
    }
  }
}
dim( MI_3_left_max )
dim(MI_3)
MI_3[100:105,]

thre2 <- 1.4
all_net <- rbind(MI_1_left_max,rbind(MI_2_left_max,MI_3_left_max))
net_left <- rbind()
z1 <- unique(all_net[,1])
for(i in 1:length(z1))
{
  k1 <- which(all_net[,1] == z1[i])
  k2 <- unique(all_net[k1,2])
  can <- all_net[k1,2]
  val <- as.numeric(all_net[k1,3])
  for(j in 1:length(k2))
  {
    zz <- which(can %in% k2[j])
    va <- max(val[zz])
    if(va > thre2)
    {
      net_left <- rbind(net_left,c(z1[i],k2[j],va))
    }
  }
}

dim(unique(net_left))
unique(net_left[,2])
write.csv(MI_1,"cancerMI_1.csv")
write.csv(MI_2,"cancerMI_2.csv")
write.csv(MI_3,"cancerMI_3.csv")
write.csv(net_left,"cancernetleft1.4.csv")


###############################################################
##### Calculate 1-MI if exists
library(pracma)
MI_1 <- rbind()
MI_2 <- rbind()
MI_3 <- rbind()
MI_re <- rbind()
Node_0 <- unlist(vertex_attr(MI_0_net))
for(i in 1:dim(MI_0_left)[1])
{
  #i <-1
  loc1 <- MI_0_left[i,1]
  loc2 <- MI_0_left[i,2]
  ### if they have shared one-order neighbor
  ne1 <- setdiff(Node_0[unlist(ego(MI_0_net, order = 1, MI_0_left[i,1]))],MI_0_left[i,1])
  ne2 <- setdiff(Node_0[unlist(ego(MI_0_net, order = 1, MI_0_left[i,2]))],MI_0_left[i,1])
  nn <- unique(setdiff(intersect(ne1,ne2),c(MI_0_left[i,1],MI_0_left[i,2])))
  if(isempty(nn)==F)
  {
    # 1-order
    for(j in 1:length(nn))
    {
      loc3 <-  nn[j]
      con <- dis_data[,loc3]
      mi <- condinformation(dis_data[,loc1],dis_data[,loc2],con,method="emp")
      MI_1 <- rbind(MI_1,c(MI_0_left[i,c(1,2)],nn[j],1,mi))
    }
    # 2-order
    if(length(nn)>1)
    {
      for(k in 1:(length(nn)-1))
      {
        loc3 <- nn[k]
        for(b in (k+1):length(nn))
        {
          loc4 <- nn[b]
          con <- c(dis_data[,loc3],dis_data[,loc4])
          mi <- condinformation(dis_data[,loc1],dis_data[,loc2],con,method="emp")
          MI_2 <- rbind(MI_2,c(MI_0_left[i,c(1,2)],nn[k],nn[b],2,mi))
        }
      }
    }#else{
    # MI_2 <- rbind(MI_2,c(MI_0_left[i,c(1,2)],'non','non',1,'left'))
    #}
    # 3-order
    if(length(nn)>2)
    {
      for(k in 1:(length(nn)-2))
      {
        loc3 <- nn[k]
        for(b in (k+1):(length(nn)-1))
        {
          loc4 <- nn[b]
          for(h in (b+1):length(nn))
          {
            loc5 <- nn[h]
            con <- c(dis_data[,loc3],dis_data[,loc4],dis_data[,loc5])
            mi <- condinformation(dis_data[,loc1],dis_data[,loc2],con,method="emp")
            MI_3 <- rbind(MI_3,c(MI_0_left[i,c(1,2)],nn[k],nn[b],nn[h],3,mi))
          }
        }
      }
    }#else{
    # MI_3 <- rbind(MI_3,c(MI_0_left[i,c(1,2)],'non','non','non',2,'left'))
    #}
    
    
  }
  else{
    #MI_1 <- rbind(MI_1,c(MI_0_left[i,c(1,2)],0,0)) 
    MI_re <-rbind( MI_re,c(MI_0_left[i,c(1,2)],0,'left'))
  }
}
colnames(MI_re) <- c('node1','node2','MI_order','states')
colnames(MI_1) <- c('node1','node2','neighbor1','MI_order','states')
colnames(MI_2) <- c('node1','node2','neighbor1','neighbor2','MI_order','states')
colnames(MI_3) <- c('node1','node2','neighbor1','neighbor2','neighbor3','MI_order','states')
#### select the maximum MI and CMI
thre <- 0.5
#### 1-order CMI
u1 <- unique(MI_1[,1])
MI_1_left_max <- rbind()
for(i in 1:length(u1))
{
  loc <- which(MI_1[,1] %in% u1[i])
  k1 <- unique(MI_1[loc,2])
  can <- MI_1[loc,2]
  for(j in 1:length(k1))
  {
    z1 <- apply(as.matrix(MI_1[loc[which(can %in% k1[j])],dim(MI_1)[2]]),
                2,function(x) as.numeric(x))
    if(max(z1)> thre)
    {
      MI_1_left_max <- rbind(MI_1_left_max,c(u1[i],k1[j],max(z1)))
    }
  }
}

#### 2-order CMI
u2 <- unique(MI_2[,1])
MI_2_left_max <- rbind()
for(i in 1:length(u2))
{
  loc <- which(MI_2[,1] %in% u2[i])
  k1 <- unique(MI_2[loc,2])
  can <- MI_2[loc,2]
  for(j in 1:length(k1))
  {
    z1 <- apply(as.matrix(MI_2[loc[which(can %in% k1[j])],dim(MI_2)[2]]),
                2,function(x) as.numeric(x))
    if(max(z1) > thre)
    {
      MI_2_left_max <- rbind(MI_2_left_max,c(u2[i],k1[j],max(z1)))
    }
  }
}
#### 3-order CMI
#max(apply(as.matrix(MI_3[,7]),2,function(x) as.numeric(x)))
u3 <- unique(MI_3[,1])
MI_3_left_max <- rbind()
dim(MI_3)
for(i in 1:length(u3))
{
  #i<-3
  loc <- which(MI_3[,1] %in% u3[i])
  k1 <- unique(MI_3[loc,2])
  can <- MI_3[loc,2]
  for(j in 1:length(k1))
  {
    #j <-4
    z1 <- apply(as.matrix(MI_3[loc[which(can %in% k1[j])],dim(MI_3)[2]]),
                2,function(x) as.numeric(x))
    if(max(z1) > thre)
    {
      MI_3_left_max <- rbind(MI_3_left_max,c(u2[i],k1[j],max(z1)))
    }
  }
}

#######
## 1-order, 2-order, 3-order if exists, select the maximum
######
#######
##set new threshold for the 
#######
thre2 <- 0.7
all_net <- rbind(MI_1_left_max,rbind(MI_2_left_max,MI_3_left_max))
net_left <- rbind()
z1 <- unique(all_net[,1])
for(i in 1:length(z1))
{
  k1 <- which(all_net[,1] == z1[i])
  k2 <- unique(all_net[k1,2])
  can <- all_net[k1,2]
  val <- as.numeric(all_net[k1,3])
  for(j in 1:length(k2))
  {
    zz <- which(can %in% k2[j])
    va <- max(val[zz])
    if(va > thre2)
    {
      net_left <- rbind(net_left,c(z1[i],k2[j],va))
    }
  }
}
#hist(as.numeric(net_left[,3]))

MI_net <- graph.data.frame(net_left[,c(1,2)],directed = F)

plot(MI_net, vertex.color="purple",vertex.size=8,
     label.font=2,label.cex=2,label.color='black',main='MI_net')

library(igraph)
cancermodule=read.csv("differentnetwork.csv")
#edgeweight=cancermodule[,4]
network <- graph.data.frame(cancermodule[,c(2,3)],directed = F)
fc <- cluster_fast_greedy(network)
a=membership(fc)
b=as.data.frame(t(sapply(a, "[", i = 1:max(sapply(a, length)))))
write.csv(b,"module_2.csv")
sizes(fc)
plot.membership<-function(graph,membership,main=""){
  V(graph)$member<-membership
  mem.col<-rainbow(length(unique(membership)),alpha=0.6)
  V(graph)$color<-mem.col[membership]
  V(graph)$size = 5
  V(graph)$label=NA
  legend(x=-1.5,y=1.5,levels(factor(V(graph)$location)),pch=21,col="#777777",pt.bg=vcolor)
  plot(graph,edge.width=E(graph)$weight,vertex.color=V(graph)$color,main=main,layout=layout.lgl)
  
}
plot.membership(network,a)
legend(1, 95, legend=c("Line 1", "Line 2"),
       col=c("red", "blue"), lty=1:2, cex=0.8,
       title="Line types", text.font=4, bg='lightblue')


stop()



