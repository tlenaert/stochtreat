#!/usr/bin/env python3

import numpy as np
import helperfunctions
import sys
import warnings
import math



#process file data
filenames = helperfunctions.natural_sort(sys.argv[2:])
rawdata=[]
for filename,eps in zip(filenames,eps_list_raw):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        filedat=np.loadtxt(filename,comments="#",skiprows=0)
    if filedat.size==0:
        continue
    rawdata.append([filedat])


