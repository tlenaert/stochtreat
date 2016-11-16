#!/usr/bin/env python3

import pandas as pd
import numpy as np

import sys

filename=sys.argv[1]
data=pd.read_csv(filename)
# print(data.columns.values)
data=data[data["IM_initial_dose"]==400]

timepointnames=["M3","M6","M12"]
for timepoint in timepointnames:
    data[timepoint]=data[timepoint+"_IS"]/data["BL_IS"]

# print(data[timepoint],data[timepoint+"_IS"],data["BL_IS"])


# print(data)
scores=["LOW","INT","HIGH"]
sokalsets=[]
for score in scores:
    sokalsets.append(data[data.sokal==score])
    # print(sokalsets[-1].shape[0])

eurosets=[]
for score in scores:
    eurosets.append(data[data.euro==score])
    # print(data[data.euro==score])

print("overall", end=" ")
for timepoint in timepointnames:
    print(data[timepoint].mean(),end=" ")
print()
    
for sokaldata,score in zip(sokalsets,scores):
    print([sokaldata[timepoint].mean()  for timepoint in timepointnames],end=" ")
    print(", #sokal",score)

for eurodata,score in zip(eurosets,scores):
    print([eurodata[timepoint].mean()  for timepoint in timepointnames],end=" ")
    # print(eurodata[timepoint].mean(),end=",")
    print(", #euro",score)


# print(data[data.sokal=="LOW"])

