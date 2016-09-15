#!/usr/bin/env python3

import numpy as np
import helperfunctions
import sys

filenames = helperfunctions.natural_sort(sys.argv[1:])
treattime_range=np.linspace(0.,5.,100.)

rawdata=[np.loadtxt(filename,comments="#") for filename in filenames]

data=np.concatenate(rawdata)
probs=data[:,0]
# print(probs)
# exit(0)


for t,p in zip(treattime_range,probs):
    print(t,p)


