#  preperation for calculating score based on jaccard
import networkx as nx
import pandas as pd

### read edge list to networkx ###

# the format of each line: (src dst whole_duration/total_duration total_duration)

def createGraph(filename) :#, w_index

    G = nx.Graph()
    
    for line in open(filename) :
        line=line.strip('\n')
        line = line.replace('"',"")

        strlist = line.split(",")

        #print(strlist)

        n1 = strlist[1]

        n2 = strlist[2]

        #weight = float(strlist[w_index])

        #G.add_weighted_edges_from([(n1, n2, weight)])
        G.add_edges_from([(n1, n2)])
    
    return G

Gnormal=createGraph("normal+1.csv")
Gcancer=createGraph("cancer+1.csv")
Gall=createGraph("edge509.csv")
Gcancerdifferentnetwork=createGraph("cancerdifferentnetwork+1.csv")
Gnormaldifferentnetwork=createGraph("normaldifferentnetwork+1.csv")
print(len(Gnormal.node()))
print(len(Gcancer.node()))
print(len(Gall.node()))
print(len(Gcancerdifferentnetwork.node()))
print(len(Gnormaldifferentnetwork.node()))

Gallnode_normalnetwork=nx.Graph()
Gallnode_cancernetwork=nx.Graph()

Gallnode_normalnetwork.add_nodes_from(Gall.node())
Gallnode_normalnetwork.add_edges_from(Gnormaldifferentnetwork.edges())
Gallnode_cancernetwork.add_nodes_from(Gall.node())
Gallnode_cancernetwork.add_edges_from(Gcancerdifferentnetwork.edges())
print(len(Gallnode_normalnetwork.node()))
print(len(Gallnode_cancernetwork.node()))
#print(Gallnode_normalnetwork.edges)

import numpy as np
A=np.array(nx.adjacency_matrix(Gall).todense())
print(A)
Normalmatrix=np.array(nx.adjacency_matrix(Gallnode_normalnetwork).todense())
print(Normalmatrix)
Cancermatrix=np.array(nx.adjacency_matrix(Gallnode_cancernetwork).todense())
print(Cancermatrix)

df1=pd.DataFrame(Normalmatrix)
df2=pd.DataFrame(Cancermatrix)

df1.to_csv("normalmatrix.csv")
df2.to_csv("cancermatrix.csv")


