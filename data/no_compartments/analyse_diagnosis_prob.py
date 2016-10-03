#!/usr/bin/env python3

import numpy as np
import helperfunctions
import sys
import os

dirnames = helperfunctions.natural_sort(sys.argv[1:])

for i,dirname in enumerate(dirnames):
    diagfracs=[]
    filenames = [fn for fn in os.listdir(dirname) if fn.endswith(".dat") ]
    for filename in filenames:
        with open(dirname+"/"+filename) as f:
            for line in f:
                if line.startswith("#Fraction diagnosed"):
                    diagfracs.append(float(line.split()[-1]))
                    break

    print(i+1,np.mean(diagfracs))

