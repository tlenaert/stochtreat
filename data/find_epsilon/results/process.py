#!/usr/bin/env python3
import numpy as np
import pandas as pd
import helperfunctions
import sys

filename=sys.argv[1]
data=pd.read_csv(filename,delim_whitespace=True)

fitted=data[["#case","epsc","epsb","z"]].as_matrix().T
print(fitted)


