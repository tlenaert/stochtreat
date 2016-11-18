#!/usr/bin/env python3

import numpy as np
import helperfunctions
import sys
import warnings
import math

# array_to_find=np.array([0.00266636,0.00119456,0.000851222])
# to_find_data=[
#     [0.11994494812140269, 0.031847672264213756, 0.025543728434327066] , #sokal LOW
#     [0.064294621947554351, 0.041933604102809392, 0.032632945359571189] , #sokal INT
#     [0.12944505110038521, 0.16560671082658152, 0.0825145962552611] , #sokal HIGH
#     [0.11902869904379794, 0.031901634379805963, 0.03729351383155937] , #euro LOW
#     [0.088553538937121576, 0.083289731014841065, 0.03001561769103259] , #euro INT
#     [0.11597804280543321, 0.11706848735364135, 0.12163073243492417]  #euro HIGH
# ]
to_find_data=[
[0.013335852927678644, 0.003168316831683168, 0.001331085278196823] , #sokal LOW
[0.016415868673050615, 0.0029578380531110836, 0.001119067968492862] , #sokal INT
[0.03188643788635723, 0.005482346453259462, 0.00285670035325997] , #sokal HIGH
[0.012729061043539474, 0.003363537748437772, 0.001529934963971043] , #euro LOW
[0.021978650945826055, 0.003143673416564866, 0.0012283644102662355] , #euro INT
[0.060362173038229376, 0.009255533199195172, 0.006306495690561278] , #euro HIGH
[0.01688620945460927, 0.003220529598200593, 0.0015214083894805476] # overall
]

to_find_names=["sokalLOW","sokalINT","sokalHIGH","euroLOW","euroINT","euroHIGH"]
array_to_find=np.array(to_find_data[int(sys.argv[1])])
name_to_find=to_find_names[int(sys.argv[1])]

delta_start=5e-1
delta_step=1e-5

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
filenames = helperfunctions.natural_sort(sys.argv[2:])
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


# for i,d,s,e in zip(range(len(rawdata)),rawdata,rawstd,eps_list):
#     if math.isclose(e[1],epsilon_i_range[32],abs_tol=1e-5):
#         print(i,d[0],s[0],e)
# exit(0)

data=np.concatenate(rawdata)
data_std=np.concatenate(rawstd)
eps_list=np.array(eps_list)
eps_list=eps_list[~np.isnan(data).any(axis=1)]
data_std= data_std[~np.isnan(data).any(axis=1)]
data= data[~np.isnan(data).any(axis=1)]

def find_all_inrange(data,array_to_find,delta):
    a=data <  array_to_find+delta
    b=data >  array_to_find-delta
    c=np.logical_and(a,b)
    
    resultindices=np.where(c.all(axis=1))[0]
    return(resultindices)

def found_close_enough(inarray,array_to_find,data_std):
    for x,a,s in zip (inarray,array_to_find,data_std):
        if x-s > a or x+s < a:
            return False
    return True


number=np.inf
d=delta_start
inx=[]
while (number > 5):
    oldinx=inx
    inx=find_all_inrange(data,array_to_find,d)
    number=len(inx)
    d=d-delta_step
    if (number<1):
        inx=oldinx
resultindices=inx
# resultindices=[]
# for i,x in enumerate(data):
#     if found_close_enough(x,array_to_find,data_std[i]):
#         resultindices.append(i)
#         d=np.sum(data_std[i])
#         # print(i,d,data_std[i],x)

result_eps=eps_list[resultindices]
resultdata=data[resultindices]

print("#finding",name_to_find,array_to_find,"precision="+str(d))
for epslist,indice in zip(result_eps,resultindices):
    for x in epslist:
        print(x,end=" ")
    for y in data[indice]:
        print(y,end=" ")
    print()


