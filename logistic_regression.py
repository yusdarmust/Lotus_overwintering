#!/usr/bin/env python3
#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.linear_model import LinearRegression
import scipy.stats as stats

# Loading data
data = pd.read_csv('../data/20200506_Lj_accession_winter_timecourse_mapped_vs_Gifu.RPKM.txt',sep='\t')

def remove_duplicates(seq):
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

# Extracing line names
lines = [n.split('_', 1)[0] for n in data.columns[1:]]
lines = remove_duplicates(lines)
lines = np.array(lines) # For easy indexing

# Mask for which lines to use. Set the following indices to 0 to exclude a line
# Gifu:0, 008:1, 013:2, 016:3, 022:4, 028:5, 044:6, 066:7, 072:8, 110:9, 111:10, 113:11
mask = np.ones(len(lines))

# Extracting gene names
genes = data.iloc[:,0]

# Extracting expression levels
expressions = data.iloc[:,1:].to_numpy()
# Reshaping expression levels into gene x line x month
expressions = expressions.reshape(-1,12,4)

# Winter Survival phenotype
# Gifu, 008, 013, 016, 022, 028, 044, 066, 072, 110, 111, 113
#phenotype = np.array([1,0,1,0,0,0,1,0,1,0,0,1])

# Latitude
# Gifu, 008, 013, 016, 022, 028, 044, 066, 072, 110, 111, 113
latitude = np.array([35.717,35.018,35.210,35.135,31.192,31.159,38.387,32.132,32.852,32.282,32.860,32.863])
lat_bool = latitude >= np.mean(latitude)

# Pop3_membership
# Gifu, 008, 013, 016, 022, 028, 044, 066, 072, 110, 111, 113
pop3_membership = np.array([86,47,93,73,65,39,96,48,88,52,31,79])
phenotype = pop3_membership >= np.mean(pop3_membership)

# Helper functions for pomegranate, which requires a list of 1d arrays
def unpack_data(x):
    for i in range(x.shape[0]):
        yield x[i,:]

pct_mod = 100/expressions.shape[0]


scores = np.zeros(len(genes))
lat_scores = np.zeros(len(genes))
misclassified = np.zeros((len(genes),len(lines)))

for row in range(expressions.shape[0]):
    features = expressions[row,:,:].reshape(-1,4)
    #detrender = LinearRegression().fit(latitude.reshape(-1,1),features)
    #features = features - detrender.predict(latitude.reshape(-1,1))
    #features = np.concatenate((latitude.reshape(-1,1),features),axis=1)
    lr = LogisticRegression(solver='liblinear',max_iter=50).fit(features, phenotype)
    scores[row] = lr.score(features,phenotype)
    misclassified[row,:] = (lr.predict(features) != phenotype)

    lr = LogisticRegression(solver='liblinear',max_iter=50).fit(features, lat_bool)
    lat_scores[row] = lr.score(features,lat_bool)
    #print(row*pct_mod)

targets = ['LotjaGi1g1v0458500.1', 'LotjaGi6g1v0343100.1']
name_filter = np.array([True if g in targets else False for g in genes])

print(genes[name_filter])
print(scores[name_filter])
print(lines[np.any(misclassified[name_filter,:],axis=0)])

print(sum(scores == 1))
print(sum(scores > 0.9))
print(sum(scores > 0.7))
#%%
'''
plt.subplot(411)
plt.bar(lines,np.sum(misclassified[scores>0.9,:],axis=0))
plt.title('Misclassified lines, Genes >0.9 Accuracy')
plt.subplot(412)
plt.bar(lines,np.sum(misclassified[scores>0.7,:],axis=0))
plt.title('Misclassified lines, Genes >0.7 Accuracy')
plt.subplot(413)
plt.bar(lines,np.sum(misclassified[scores<0.7],axis=0))
plt.title('Misclassified lines, Genes <0.7 Accuracy')
plt.subplot(414)
plt.bar(lines,np.sum(misclassified,axis=0))
plt.title('Misclassified lines, All Genes')
plt.show()
'''
df = pd.DataFrame()
df['Gene'] = genes
df['Surv_Score'] = scores
df['Lat_Score'] = lat_scores
df.to_csv('../data/lr_surv_lat_scores.txt')
print(df)

#%%


plt.figure()
plt.hist2d(scores, lat_scores, density=False, cmap='plasma')
plt.colorbar(label='Number of Genes')
plt.plot(scores[name_filter],lat_scores[name_filter],'xr')
#plt.hist(scores)
#plt.vlines(scores[name_filter],ymin=0, ymax = 10000,colors='r')
plt.title('Latitude Score vs Survival Score')
plt.xlabel('Survival Score')
plt.ylabel('Latitude Score')
#plt.xlim([0,1])
#plt.ylim([0,10000])
plt.legend(['Genes of Interest','Accuracy Distribution'])
plt.show()

'''
# Filtering by gene of interest for preliminary assessment
targets = ['LotjaGi1g1v0458500.1', 'LotjaGi6g1v0343100.1']
name_filter = np.array([True if g in targets else False for g in genes])
expressions = expressions[name_filter,:,:]
'''


# %%
