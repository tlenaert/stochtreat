#!/usr/bin/env python3

import pandas as pd
import numpy as np

import sys

filename=sys.argv[1]
data=pd.read_csv(filename)
# print(data.columns.values)

timepointnames=["M3","M6","M12"]
for timepoint in timepointnames:
    data[timepoint]=data[timepoint+"_IS"]/data["BL_IS"]

# print(data[timepoint],data[timepoint+"_IS"],data["BL_IS"])


# print(data)
scores=["LOW","INT","HIGH"]
sokalsets=[]
for score in scores:
    sokalsets.append(data[data.sokal==score])

eurosets=[]
for score in scores:
    eurosets.append(data[data.euro==score])
    # print(data[data.euro==score])

for sokaldata,score in zip(sokalsets,scores):
    print("sokal",score, end=" ")
    for timepoint in timepointnames:
        print(sokaldata[timepoint].mean(),end=" ")
    print()

for eurodata,score in zip(eurosets,scores):
    print("euro",score, end=" ")
    for timepoint in timepointnames:
        print(eurodata[timepoint].mean(),end=" ")
    print()


# print(data[data.sokal=="LOW"])

