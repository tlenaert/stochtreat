#!/usr/bin/env python3

import numpy as np
import helperfunctions
import sys

# array_to_find=np.array([0.00266636,0.00119456,0.000851222])
array_to_find=np.array([3.016534688,2.210308663, 2.04162678])
delta=2.0e-0

##get tested epsilon values
epsilon_min=0.5
epsilon_max=1.0
epsilonvalues=30
epsilon_i_range=np.linspace(epsilon_min,epsilon_max,epsilonvalues)
epsilon_c_range=np.linspace(epsilon_min,epsilon_max,epsilonvalues)
eps_list=[]
for epsilon_c in epsilon_c_range:
    for epsilon_i in epsilon_i_range:
        if epsilon_i<epsilon_c:
            continue
        eps_list.append([epsilon_c,epsilon_i])
eps_list=np.array(eps_list)

#process file data
filenames = helperfunctions.natural_sort(sys.argv[1:])
rawdata=[[np.loadtxt(filename,comments="#")] for filename in filenames]
data=np.concatenate(rawdata)
eps_list=eps_list[~np.isnan(data).any(axis=1)]
data= data[~np.isnan(data).any(axis=1)]

a=data <  array_to_find+delta
b=data >  array_to_find-delta
c=np.logical_and(a,b)

resultindices=np.where(c.all(axis=1))
print(resultindices)
result=eps_list[resultindices]
resultdata=data[resultindices]
print("#finding ",array_to_find)
for y,yresult in zip(result,resultindices):
    for x in y:
        print(x,end=" ")
    print()


