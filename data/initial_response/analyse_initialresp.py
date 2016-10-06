#!/usr/bin/env python3

import numpy as np
import sys

filenames = sys.argv[1:]
t_max=20.
t_steps=100

rawdata=[np.loadtxt(filename,comments="#") for filename in filenames]

data=np.concatenate(rawdata)

# print(np.where(data[:,1]==1))
relapse=data[np.where(data[:,1]==1)]
norelapse=data[np.where(data[:,1]==0)]

np.savetxt("relapsedat.txt",relapse)
np.savetxt("norelapsedat.txt",norelapse)
# for x in relapse:
#     for y in x:
#         print(y, end=" ")
#     print()
#
# for x in norelapse:
#     for y in x:
#         print(y, end=" ")
#     print()
