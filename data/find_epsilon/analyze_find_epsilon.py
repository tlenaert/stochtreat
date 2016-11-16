#!/usr/bin/env python3

import numpy as np
import helperfunctions
import sys

# array_to_find=np.array([0.00266636,0.00119456,0.000851222])
to_find_data=[
    [0.11994494812140269, 0.031847672264213756, 0.025543728434327066] , #sokal LOW
    [0.064294621947554351, 0.041933604102809392, 0.032632945359571189] , #sokal INT
    [0.12944505110038521, 0.16560671082658152, 0.0825145962552611] , #sokal HIGH
    [0.11902869904379794, 0.031901634379805963, 0.03729351383155937] , #euro LOW
    [0.088553538937121576, 0.083289731014841065, 0.03001561769103259] , #euro INT
    [0.11597804280543321, 0.11706848735364135, 0.12163073243492417]  #euro HIGH
]
to_find_names=["sokalLOW","sokalINT","sokalHIGH","euroLOW","euroINT","euroHIGH"]
array_to_find=np.array(to_find_data[int(sys.argv[1])])
name_to_find=to_find_names[int(sys.argv[1])]

delta_start=5e-1
delta_step=1e-5

##get tested epsilon values
epsilon_min=0.5
epsilon_max=1.0
epsilonvalues=50
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
filenames = helperfunctions.natural_sort(sys.argv[2:])
rawdata=[[np.loadtxt(filename,comments="#")] for filename in filenames]
data=np.concatenate(rawdata)
eps_list=eps_list[~np.isnan(data).any(axis=1)]
data= data[~np.isnan(data).any(axis=1)]

def find_closest_values(data,array_to_find,delta):
    a=data <  array_to_find+delta
    b=data >  array_to_find-delta
    c=np.logical_and(a,b)
    
    resultindices=np.where(c.all(axis=1))[0]
    return(resultindices)

number=np.inf
d=delta_start
inx=[]
while (number > 1):
    oldinx=inx
    inx=find_closest_values(data,array_to_find,d)
    number=len(inx)
    d=d-delta_step
    if (number<1):
        inx=oldinx
resultindices=inx

result_eps=eps_list[resultindices]
resultdata=data[resultindices]

print("#finding",name_to_find,array_to_find,"precision="+str(d))
for epslist,indice in zip(result_eps,resultindices):
    for x in epslist:
        print(x,end=" ")
    for y in data[indice]:
        print(y,end=" ")
    print()


