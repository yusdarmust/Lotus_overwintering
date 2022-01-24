import pandas as pd
import sys

p1, p2 = sys.argv[1:3]

df = pd.read_csv(p1,
                 header=None,
                 comment='#',
                 sep='\t')
edf = pd.DataFrame(index=df.index,
                   columns=df.columns[9:],
                   dtype='str')
sdf = df.iloc[:, 9:]

for i in range(len(sdf.index)):
    for j in range(len(sdf.columns)):
        c = sdf.iloc[i, j]
        c = c.split(':')
        c_1 = c[1].split(',')
        if c[0] == '0/1':
            if c_1[1] == '0':
                if int(c_1[2]) - int(c_1[1]) <= 5:
                    edf.iloc[i, j] = '1/1' + ':' + c[1]
                elif int(c_1[0]) - int(c_1[1]) <= 5:
                    edf.iloc[i, j] = '0/0' + ':' + c[1]
                else:
                    edf.iloc[i, j] = '0/1' + ':' + c[1]
            if c_1[0] == '0':
                edf.iloc[i, j] = '0/0' + ':' + c[1]
            if c_1[2] == '0':
                edf.iloc[i, j] = '1/1' + ':' + c[1]
            elif '0' not in c_1:
                edf.iloc[i, j] = './.' + ':' + '0,0,0'
        else:
            edf.iloc[i, j] = c[0] + ':' + c[1]

edf = pd.concat((df.iloc[:, :9], edf), axis=1)
edf.to_csv(p2,
           sep='\t',
           index=False)
