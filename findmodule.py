import pandas as pd
import networkx as nx

cancer=pd.read_csv("cancer+1.csv",index_col=0,header=0)
normal=pd.read_csv("normal+1.csv",index_col=0,header=0)
cancernetwork=cancer.iloc[:,0:2]
normalnetwork=normal.iloc[:,0:2]




def node():
    cancernode=cancernetwork.iloc[:,0].append(cancernetwork.iloc[:,1])
    cancernode=cancernode.drop_duplicates( keep='first')
    print(len(cancernode))
    normalnode=normalnetwork.iloc[:,0].append(normalnetwork.iloc[:,1])
    normalnode=normalnode.drop_duplicates( keep='first')
    print(len(normalnode))
    allnode=cancernode.append(normalnode)
    print(len(allnode))
    index=allnode.duplicated()
    commonnode=allnode[index]
    print(len(commonnode))
    return commonnode


def cancer_normal_unique_edge():
    #cancerdifferentnetwork = pd.read_csv("cancerdifferentnetwork.csv", index_col=0, header=0)
    #normaldifferentnetwork = pd.read_csv("normaldifferentnetwork.csv", index_col=0, header=0)
    cancerdifferentnetwork = pd.read_csv("cancerdifferentnetwork+1.csv", index_col=0, header=0)
    normaldifferentnetwork = pd.read_csv("normaldifferentnetwork+1.csv", index_col=0, header=0)
    differentnetwork = pd.read_csv("edge509.csv", index_col=0, header=0)
    cancerdifferentnetwork = cancerdifferentnetwork.values.tolist()
    normaldifferentnetwork = normaldifferentnetwork.values.tolist()
    commonnetwork = differentnetwork.values.tolist()

    cancerleftunique = []
    normalleftunique = []

    for index,i in enumerate(commonnetwork):
        if i not in cancerdifferentnetwork:
            print(i)
            normalleftunique.append(i)
            commonnetwork[index].append(1)
        else:
            print(index)
            commonnetwork[index].append(2)

    print(len(normalleftunique))
    normalleftunique = pd.DataFrame(normalleftunique)
    #normalleftunique.to_csv("normalleftuniquel.csv")
    commonnetwork=pd.DataFrame(commonnetwork)
    commonnetwork.to_csv("commonnetworkwithmark1.csv")

    commonnetwork = differentnetwork.values.tolist()
    for index,i in enumerate(commonnetwork):
        if i not in normaldifferentnetwork:
            cancerleftunique.append(i)
            commonnetwork[index].append(0)
        else:
            commonnetwork[index].append(2)

    print(len(cancerleftunique))
    cancerleftunique = pd.DataFrame(cancerleftunique)
    #cancerleftunique.to_csv("cancerleftuniquel.csv")
    commonnetwork=pd.DataFrame(commonnetwork)
    commonnetwork.to_csv("commonnetworkwithmark2.csv")




def mergeORdifferentnetwork(commonnode): #全网络或者共同的网络或者同节点不同的网络
    mrel=cancernetwork.append(normalnetwork)
    print("mrel",mrel.head(),len(mrel))
    index=mrel.duplicated()
    common=mrel[index]
    #print(common)
    newrelunique = mrel.drop_duplicates( keep='first')
    newrelunique.to_csv("mergenetwork.csv")

    differentnetwork=pd.DataFrame()
    for i in range(0,len(newrelunique.iloc[:,0])):
        if newrelunique.iloc[i,0] in list(commonnode):
            if newrelunique.iloc[i,1] in list(commonnode):
                differentnetwork=differentnetwork.append(newrelunique.iloc[i,:])
    print(differentnetwork,len(differentnetwork))
    differentnetwork.to_csv("differentnetwork.csv")


def cancerunique(cancernetwork, normalnetwork,common):
    cancerunique = []

    cancernetwork = cancernetwork.values.tolist()
    normalnetwork = normalnetwork.values.tolist()
    commonnetwork = common.values.tolist()

    for i in cancernetwork:
        if i not in commonnetwork:
            cancerunique.append(i)

    print(len(cancerunique))

    df = pd.DataFrame(cancerunique)
    df.to_csv("cancerunique.csv")

if __name__=='__main__':
    commonnode=node()
    # cancerdifferentnetwork = pd.DataFrame()
    # for i in range(0, len(cancernetwork.iloc[:, 0])):
    #     if cancernetwork.iloc[i, 0] in list(commonnode) or cancernetwork.iloc[i, 1] in list(commonnode):
    #             cancerdifferentnetwork = cancerdifferentnetwork.append(cancernetwork.iloc[i, :])
    # cancerdifferentnetwork=cancerdifferentnetwork.drop_duplicates()
    # print(cancerdifferentnetwork, len(cancerdifferentnetwork))
    # cancerdifferentnetwork.to_csv("cancerdifferentnetwork+1.csv")
    # normaldifferentnetwork = pd.DataFrame()
    # for i in range(0, len(normalnetwork.iloc[:, 0])):
    #     if normalnetwork.iloc[i, 0] in list(commonnode) or normalnetwork.iloc[i, 1] in list(commonnode):
    #         normaldifferentnetwork = normaldifferentnetwork.append(normalnetwork.iloc[i, :])
    #
    # normaldifferentnetwork=normaldifferentnetwork.drop_duplicates()
    # print(normaldifferentnetwork, len(normaldifferentnetwork))
    # normaldifferentnetwork.to_csv("normaldifferentnetwork+1.csv")
    #
    # edge509=cancerdifferentnetwork.append(normaldifferentnetwork)
    # edge509=edge509.drop_duplicates()
    # print(edge509, len(edge509))
    # edge509.to_csv("edge509.csv")
    #
    cancer_normal_unique_edge()



