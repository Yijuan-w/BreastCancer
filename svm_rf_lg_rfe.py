import numpy as np
from sklearn.feature_selection import RFE, RFECV
from sklearn.svm import SVC, SVR
from sklearn.svm import LinearSVC
from sklearn import model_selection
import pandas as pd
from imblearn.under_sampling import RandomUnderSampler
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn import metrics
from matplotlib import pyplot as plt
from imblearn.combine import SMOTEENN
from sklearn.model_selection import train_test_split
from sklearn.metrics import confusion_matrix,classification_report
from sklearn.metrics import f1_score,precision_score,recall_score,roc_auc_score,accuracy_score,roc_curve
from sklearn.preprocessing import StandardScaler,RobustScaler,MinMaxScaler
ee = SMOTEENN()

ss = RobustScaler()


class RandomForestClassifierWithCoef(RandomForestClassifier):
    def fit(self, *args, **kwargs):
        super(RandomForestClassifierWithCoef, self).fit(*args, **kwargs)
        self.coef_ = self.feature_importances_


rf = RandomForestClassifierWithCoef(n_estimators=500, min_samples_leaf=5, n_jobs=-1)


def svmRFE(features,label):
    ee = SMOTEENN()
    features, label = ee.fit_sample(features, label)
    label.astype(np.int)
    X_train, X_test, y_train, y_test = train_test_split(features, label, test_size=0.3, random_state=66)

    y = y_train
    #print(y.head())
    x = X_train

    print(features.head())
    print('开始训练特征')
    svc = SVC(kernel="linear", C=1, probability=True)

    rfe = RFE(estimator=svc, n_features_to_select=3, step=1)
    rfe.fit(x, y)
    print('训练特征选择完成，输出特征：')
    result = features.columns[rfe.get_support()]
    print(result)
    print(rfe.ranking_)

    re = RFECV(estimator=svc, step=1, cv=5)
    re.fit(x, y)
    result2 = features.columns[re.get_support()]
    y_pred2 = re.predict(X_test)
    y_proba2 = re.predict_proba(X_test)[:, 1]
    print(classification_report(y_test, y_pred2, digits=2))
    print(result2)
    print('特征输出完毕')
    fpr1, tpr1, thresholds1 = roc_curve(y_true=y_test, y_score=y_proba2)
    return fpr1,tpr1

def rfRFE(features,label):
    ee = RandomUnderSampler(random_state=0)
    features, label = ee.fit_sample(features, label)
    label.astype(np.int)
    X_train, X_test, y_train, y_test = train_test_split(features, label, test_size=0.3, random_state=66)

    y = y_train
    #print(y.head())
    x = X_train

    print(features.head())
    print('开始训练特征')
    rf = RandomForestClassifierWithCoef(n_estimators=500, min_samples_leaf=5, n_jobs=-1)
    re = RFECV(estimator=rf, step=1, cv=10)
    re.fit(x, y)
    result2 = features.columns[re.get_support()]
    y_pred2 = re.predict(X_test)
    y_proba2 = re.predict_proba(X_test)[:, 1]
    print(classification_report(y_test, y_pred2, digits=2))
    print(result2)
    print(re.ranking_)
    print('特征输出完毕')
    fpr1, tpr1, thresholds1 = roc_curve(y_true=y_test, y_score=y_proba2)
    return fpr1,tpr1

def lgRFE(features,label):
    ee = RandomUnderSampler(random_state=6)
    features, label = ee.fit_sample(features, label)
    label.astype(np.int)
    X_train, X_test, y_train, y_test = train_test_split(features, label, test_size=0.5, random_state=6)

    y = y_train
    #print(y.head())
    x=X_train
    #x = ss.fit_transform(X_train)
    #X_test=ss.fit_transform(X_test)

    print(features.head())
    print('开始训练特征')
    rf = LogisticRegression(max_iter=10000)
    re = RFECV(estimator=rf, step=1, cv=10)
    re.fit(x, y)
    result2 = features.columns[re.get_support()]
    y_pred2 = re.predict(X_test)
    y_proba2 = re.predict_proba(X_test)[:, 1]
    print(classification_report(y_test, y_pred2, digits=2))
    print(result2)
    print(re.ranking_)
    print('特征输出完毕')
    fpr1, tpr1, thresholds1 = roc_curve(y_true=y_test, y_score=y_proba2)
    return fpr1,tpr1

df=pd.read_table("HiSeqV2_BRCA_outcome.txt",header=0,index_col=0)
df=df.T
#df = pd.read_csv('E:\Pycharm\PycharmLearning\RF\\RFEtestyu1012.csv', encoding='gbk')
print(df.head())

module1=[]
module2=[]
module3=[]
module4=[]
module5=[]
moduleinfo=pd.read_csv("module_2.csv",index_col=None,header=None)
moduleinfo=moduleinfo.T



print(moduleinfo.head())
for i in range(0,len(moduleinfo.iloc[:,0])):
    if moduleinfo.iloc[i,1]=="1":
        module1.append(moduleinfo.iloc[i,0])
    if moduleinfo.iloc[i,1]=="2":
        module2.append(moduleinfo.iloc[i,0])
    if moduleinfo.iloc[i,1]=="3":
        module3.append(moduleinfo.iloc[i,0])
    if moduleinfo.iloc[i,1]=="4":
        module4.append(moduleinfo.iloc[i,0])
    if moduleinfo.iloc[i, 1] == "5":
        module5.append(moduleinfo.iloc[i, 0])
print("module1",module4)

for i in range(1,6):
    if i == 1:
        features = df[module1]
        print(len(module1))
        label = df.Label
        fpr1,tpr1=svmRFE(features,label)
    if i == 2:
        features = df[module2]
        print(len(module2))
        label = df.Label
        fpr2,tpr2=svmRFE(features,label)
    if i == 3:
        features = df[module3]
        print(len(module3))
        label = df.Label
        fpr3,tpr3=svmRFE(features, label)
    if i == 4:
        features = df[module4]
        print(len(module4))
        label = df.Label
        fpr4, tpr4 = svmRFE(features, label)
    if i == 5:
        features = df[module5]
        print(len(module5))
        label = df.Label
        fpr5, tpr5 = svmRFE(features, label)

plt.style.use('seaborn-darkgrid')
Font={'size':18, 'family':'Times New Roman'}
Font2={'size':17, 'family':'Times New Roman'}
plt.figure(figsize=(8,6))
plt.tick_params(labelsize=15)
plt.xlim(0,1) ##设定x轴的范围
plt.ylim(0.0,1.1) ## 设定y轴的范围
plt.xlabel('False Postive Rate',Font2)
plt.ylabel('True Postive Rate',Font2)
roc_auc1 = metrics.auc(fpr1, tpr1)
roc_auc2 = metrics.auc(fpr2, tpr2)
roc_auc3 = metrics.auc(fpr3, tpr3)
roc_auc4 = metrics.auc(fpr4, tpr4)
roc_auc5 = metrics.auc(fpr5, tpr5)
plt.plot(fpr1, tpr1, 'b', label='Module1 = %0.3f' % roc_auc1 , color='Red',linewidth=2)
plt.plot(fpr2, tpr2, 'b', label='Module2 = %0.3f'% roc_auc2 , color='k',linewidth=2)
plt.plot(fpr3, tpr3, 'b', label='Module3 = %0.3f' % roc_auc3, color='RoyalBlue',linewidth=2)
plt.plot(fpr4,tpr4,label='Module4 = %0.3f' % roc_auc4, color='seagreen',linewidth=2)
plt.plot(fpr5,tpr5,label='Module5 = %0.3f' % roc_auc5, color='gold',linewidth=2)
plt.legend(loc = 'lower right', prop=Font)

plt.show()



