library(KEGGgraph)
library(readr)
library(org.Hs.eg.db)
# 这个文件是用来提取调节通路中所有的基因的



filename <- paste (c('extdata/'),c('hsa05224.xml'), sep ="", collapse = NULL)
    
toyKGML <- system.file(filename,package="KEGGgraph")
    
toygraph <- parseKGML2Graph(toyKGML,expandGenes=TRUE)
    
toypathway <- parseKGML(toyKGML)
    
pathName <- getName(toypathway)
    
    
nodeData <- getKEGGnodeData(toygraph)
    
n <- length(nodeData)
    
if(n > 0){
  
  for(i in 1:n){
    x <-nodeData[[i]]
    
    nodename <- getName(x)
    
    nodeEntryID <- getEntryID(x)
    
    nodeType <- getType(x)
    
    
    
    if(nodeType=="gene"){
      if(1>0){
        
        
        p <- nodename
        q <- nodeEntryID
        
        if(require(org.Hs.eg.db)) {
          ioGeneID <- translateKEGGID2GeneID(nodeEntryID)
          nodesNames <- sapply(mget(ioGeneID, org.Hs.egSYMBOL, ifnotfound=NA), "[[",1)
          #nodesNames <- sapply(mget(ioGeneID, org.Mm.egSYMBOL, ifnotfound=NA), "[[",1)
        } else {
          nodesNames <- names(ios)
        }
        
        a_list <- list(
          NAME=nodesNames[[1]],
          TYPE=nodeType
        )
        print(1)
        write.table(a_list,"pathwaygene.txt",append=TRUE,sep=",")
      }
      
    } 
  }
}
