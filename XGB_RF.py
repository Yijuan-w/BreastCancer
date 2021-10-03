from sklearn.datasets import load_iris
import xgboost as xgb
from xgboost import plot_importance
from matplotlib import pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import confusion_matrix,classification_report
import pandas as pd
from imblearn.combine import SMOTEENN
from sklearn.metrics import f1_score,precision_score,recall_score,roc_auc_score,accuracy_score,roc_curve
from sklearn.model_selection import StratifiedShuffleSplit

#测试其他模型
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.naive_bayes import MultinomialNB
from sklearn.svm import LinearSVC


from sklearn.model_selection import GridSearchCV

from sklearn.svm import LinearSVC
from sklearn.preprocessing import StandardScaler,RobustScaler,MinMaxScaler
ee = SMOTEENN()

ss = MinMaxScaler()



#加载样本数据集
data=pd.read_table("HiSeqV2_BRCA_outcome.txt",header=0,index_col=0)

data=data.T
print(data.head())
y=data["Label"]
module1=[]
module2=[]
module3=[]
module4=[]
moduleinfo=pd.read_csv("module.csv",index_col=None,header=None)
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
print("module1",module1)

features = data[['BRCA1', 'BRCA2']]
label = y
#print(features)
#print(label)

ee = SMOTEENN()
features, label = ee.fit_sample(features, label)
X_train, X_test, y_train, y_test = train_test_split(features, label, test_size=0.2, random_state=66)

ss = StandardScaler()
X_train = ss.fit_transform(X_train)
X_test = ss.transform(X_test)


y_train = y_train.values
y_train = y_train.ravel()
# y_train.to_csv("y_train")

# 训练模型
model = xgb.XGBClassifier(max_depth=3, min_child_weight=3, learning_rate=0.1, n_estimators=200,
                          objective='binary:logistic', gamma=1e-5)

X_train = ss.fit_transform(X_train)
X_test = ss.transform(X_test)
model.fit(X_train, y_train)
print(X_train,X_test)
# 对测试集进行预测
y_pred = model.predict(X_test)
y_proba = model.predict_proba(X_test)[:,1]

# 计算准确率
accuracy = accuracy_score(y_test, y_pred)
precision = precision_score(y_test, y_pred)

print('accuracy:%2.f%%' % (accuracy * 100))
print('precision:%2.f%%' % (precision * 100))

print(confusion_matrix(y_test, y_pred))
print(classification_report(y_test, y_pred, digits=2))
print(roc_auc_score(y_test,y_pred))
fpr, tpr, thresholds=roc_curve(y_true=y_test,y_score=y_proba)
plt.plot(fpr,tpr,marker = 'o')

plt.show()
# 显示重要特征
plot_importance(model)
plt.show()

model2 = RandomForestClassifier()
model2.fit(X_train, y_train)
y_pred2 = model2.predict(X_test)
y_proba2 = model2.predict_proba(X_test)[:,1]
print(classification_report(y_test, y_pred2, digits=2))
fpr2, tpr2, thresholds2=roc_curve(y_true=y_test,y_score=y_proba2)
plt.plot(fpr2,tpr2)
plt.show()
print("##################")

