import pandas as pd

# expdata=pd.read_csv("HiSeqV2_BRCA_outcome.txt",header=0,index_col=0,sep="\t")
# cancernetwork=pd.read_csv()
# normalnetwork=pd.read_csv()
#
# #找相同的，然后把相同的从疾病里边删掉
# uniquenetwork=pd.DataFrame()
# for i in range(0,len(cancernetwork.iloc[:,0])):
#     if cancernetwork.iloc[i,:] not in normalnetwork:
#         cancernetwork.append(cancernetwork[i,:])
#
# #读到cytoscape 画出来看看
# #分模块
# trrust=pd.read_csv()


humandatabase=pd.read_csv("human.sourceall",sep="\t",header=None, low_memory=False)
#print(humandatabase)
genenode=pd.read_csv("pathwaygene.txt",header=None)
#print(genenode.head())
trrust=pd.read_csv("trrust_rawdata.human.tsv",sep="\t",header=None)



def rel():


    rel=[]
    for i in range(0,len(humandatabase.iloc[:,0])):
        #print(list(genenode.iloc[:,0]),humandatabase.iloc[i,0])
        if humandatabase.iloc[i,0] in list(genenode.iloc[:,0]):
            rel.append(humandatabase.iloc[i, :])

        if humandatabase.iloc[i,2] in list(genenode.iloc[:,0]):

            rel.append(humandatabase.iloc[i,:])
    print(len(rel))
    rel=pd.DataFrame(rel)
    rel.to_csv("breastgenerel.csv")

def trrustrel():
    trrel=[]
    for i in range(0,len(trrust.iloc[:,0])):
        #print(list(genenode))
        if trrust.iloc[i,0] in list(genenode.iloc[:,0]):
            trrel.append(trrust.iloc[i, :])

        if trrust.iloc[i,1] in list(genenode.iloc[:,0]):

            trrel.append(trrust.iloc[i,:])
    print(len(trrel))
    trrel=pd.DataFrame(trrel)
    trrel.to_csv("breastgenetrrel.csv")

def relall():#regnetwork and trrust combine

    reg=pd.read_csv("breastgenerel.csv",index_col=0)

    trr=pd.read_csv("breastgenetrrel.csv",index_col=0)
    reg2=reg[['0', '2']]
    print("reg",reg2.head(),len(reg2))

    trr2=trr[['0', '1']]
    trr2.columns = list(reg2)
    print("trr",trr2.head(),len(trr2))
    mrel=reg2.append(trr2)
    print("mrel",mrel.head())

    newrelunique = mrel.drop_duplicates(subset=['0', '2'], keep='first')
    print(len(newrelunique))
    newrelunique.to_csv("relall+1.csv")

if __name__=='__main__':
    # rel()
    # trrustrel()
    # relall()
    edge509=pd.read_csv("edge509.csv",index_col=0)
    print(edge509)
    node=pd.read_csv("figure/node.txt",header=None)
    print(node)
    nodediffedge=[]
    for i in range(0,len(edge509.iloc[:,0])):
        #print(list(genenode))
        if edge509.iloc[i,0] in list(node.iloc[:,0]):
            nodediffedge.append(edge509.iloc[i, :])

        if edge509.iloc[i,1] in list(node.iloc[:,0]):

            nodediffedge.append(edge509.iloc[i,:])

    print(len(nodediffedge))
    trrel=pd.DataFrame(nodediffedge)
    trrel.to_csv("nodediffedge.csv")
