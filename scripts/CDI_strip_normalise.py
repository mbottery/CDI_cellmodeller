import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import scipy.stats as stats
import pandas as pd

def segment_norm(list,sigma,size):
    if (size*2)+1>len(list):
        print("ERROR: segement size larger than list")
        return 0
    m = []
    for i in range(len(list)):
        segment = []
        centre = i
        idx_1 = i-(size)
        idx_2 = i+(size+1)
        if idx_1 < 0:
            segment=segment+list[idx_1:]
            segment=segment+list[:idx_2]
        if idx_1 >= 0 and idx_2 <= len(list):
            segment=segment+list[idx_1:idx_2]
        if idx_2 > len(list):
            idx_2 = idx_2-len(list)
            segment=segment+list[idx_1:]
            segment=segment+list[:idx_2]
        # segment of data
        segment = np.asarray(segment)
        # replace zeros with -1 so that we can multiply values by weightings
        segment[segment==0]=-1
        # calculate the weightings using normal probability density function
        weightings = scipy.stats.norm(0,sigma).pdf(range(-size,size+1))
        # if weighted mean is less than 0 it is type 0
        if np.mean(segment*weightings)<0:
            m.append(0)
        # if weighted mean is greater than 0 it is type 1
        else:
            m.append(1)
    return(m)

def find_patches(strip):
    patches =[]
    cT =[]
    patch = 1
    for i in range(len(strip)):
        if i != len(strip)-1:
            if strip[i]==strip[i+1]:
                patch+=1
            if strip[i]!=strip[i+1]:
                patches.append(patch)
                cT.append(strip[i])
                patch = 1
        if i == len(strip):
            patches.append(patch)
            cT.append(strip[i])
    return patches,cT

# read in the data
#data = pd.read_csv("cdi_cell_data.csv")

# make new column of the interaction between group, rep and radius
# so I can index each annuli of each simulation separately
data['sim']=data['group']+'.'+data['rep'].map(str)+'.'+data['radius'].map(str)

print(len(data.index))
data=data[data['density']==20000]
print(len(data.index))


block_size = 25
sigma = 3

normalised = pd.DataFrame(columns=data.columns.values)

print('Normalise')
j=0
max = len(data['sim'].unique())
for i in data['sim'].unique():
    print('{0} of {1}'.format(j,max))
    j+=1
    df = data[data['sim']==i].copy()
    cellType = df['cellType'].tolist()
    norm_cellType = segment_norm(cellType,sigma,block_size)
    df['cellType']=norm_cellType
    normalised=normalised.append(df)

normalised.to_csv('CDI_strips_normalised_3.csv',sep=',',index=False)
