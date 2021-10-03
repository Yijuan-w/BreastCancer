import pandas as pd
import numpy as np
import scipy.spatial.distance as dist

dfnormal=pd.read_csv("normalmatrix.csv",index_col=0)
dfcancer=pd.read_csv("cancermatrix.csv",index_col=0)

distlist=[]
for i in range(0,len(dfnormal)):
    v1 = dfnormal.loc[i]
    v2 = dfcancer.loc[i]

    matv = np.array([v1, v2])
    #print(matv)
    ds = dist.pdist(matv, 'jaccard')
    distlist.append(ds[0])
print(distlist)
dfdistlist=pd.DataFrame(distlist)
dfdistlist.to_csv("distlist.csv")

# 输出
# [[1 1 0 1 0 1 0 0 1] [0 1 1 0 0 0 1 1 1]]

# [ 0.75]
