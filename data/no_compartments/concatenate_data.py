#!/usr/bin/env python3

import numpy as np
import sys

filenames = sys.argv[1:]

rawdata=[np.loadtxt(filename,comments="#") for filename in filenames]

data=np.concatenate(rawdata)

for x in data:
    for y in x:
        print(y, end=" ")
    print()


