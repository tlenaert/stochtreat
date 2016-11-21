#!/usr/bin/env python3

import numpy as np
import helperfunctions
import sys
import warnings
import math


##get tested epsilon values
epsilonvalues=60
epsilon_min=0.5
epsilon_max=1.0
epsilon_n=0.85
epsilon_i_range=np.linspace(0.8,epsilon_max,epsilonvalues)
epsilon_c_range=np.linspace(0.55,0.80,epsilonvalues)
eps_list=[]
for epsilon_c in epsilon_c_range:
    for epsilon_i in epsilon_i_range:
        if epsilon_i<epsilon_c:
            continue
        eps_list.append([epsilon_c,epsilon_i])
eps_list_raw=np.array(eps_list)
eps_list=[]

#process file data
filenames = helperfunctions.natural_sort(sys.argv[1:])
rawdata=[]
rawstd=[]
for filename,eps in zip(filenames,eps_list_raw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        filedat=np.loadtxt(filename,comments="#",skiprows=0)[:3]
        filestd=np.loadtxt(filename,comments="#",skiprows=0)[3:]
    if filedat.size==0:
        continue
    rawdata.append([filedat])
    rawstd.append([filestd])
    eps_list.append(eps)


for i,d,s,e in zip(range(len(rawdata)),rawdata,rawstd,eps_list):
    for eps in e:
        print(eps,end=" ")
    for med in d[0]:
        print(med,end=" ")
    for std in s[0]:
        print(std,end=" ")
    print()
