#!/usr/bin/env python3
# from __future__ import with_statement
# from __future__ import print_function
import subprocess
import sys 
import os
import glob
import shutil
import datetime
from math import sqrt,pow
import numpy as np
import time

folder_name_prefix="find_epsilon"
no_patients=5000
treattest=0
treattime=2
output="3timepointaverage"

epsilonvalues=60
if len(sys.argv[:]) > 1 :
    epislonstep=sys.argv[1]
epsilon_min=0.5
epsilon_max=1.0
epsilon_n=0.85
epsilon_i_range=np.linspace(0.80,epsilon_max,epsilonvalues)
epsilon_c_range=np.linspace(0.55,0.80,epsilonvalues)

run_script_name="run.sh"

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def copy_tools(dest):
    scriptname=__file__
    shutil.copy2(scriptname,dest)
    shutil.copy2("./stochtreat",dest)
    shutil.copy2("./analyze_find_epsilon.py",dest)
    shutil.copy2(run_script_name,dest)



folder=folder_name_prefix+"_"+datetime.datetime.now().strftime("%m_%d_%H_%M")+"_data"
if os.path.exists(folder):
    print("folder",folder,"already exists")
    exit(0)

logfilefoldername="logfiles/"
subprocess.call(["mkdir","-p",folder])
copy_tools(folder)
    
id=0



with cd(folder):
    subprocess.call(["mkdir","-p",logfilefoldername])#create logfile foldername

    for epsilon_c in epsilon_c_range:
        for epsilon_i in epsilon_i_range:
            if epsilon_i<epsilon_c:
                continue
            parameters= "id="+str(id)+",patients="+str(no_patients)+",output="+str(output)+",treattest="+str(treattest)+",treattime="+str(treattime)+",epsc="+str(epsilon_c)+",epsb="+str(epsilon_i)+",epsn="+str(epsilon_n)
            id=id+1
            jobname= parameters
            walltime_parameter="10:30:00"#"walltime="+
            # with cd(dirname):
            print( "submitting script with",parameters)
            process=subprocess.Popen(["sbatch","--export="+parameters,"--job-name="+jobname,"-t",walltime_parameter,"-o",logfilefoldername+jobname+".out",run_script_name])
            time.sleep(0.05)
            with open("parameter_values.txt","a") as pfile:
                pfile.write(parameters+"\n")
                # status = process.wait()
                #~ print (status)


print("submission successfull, data target folder: ", folder)

#~ shutil.rmtree(folder)

# os.remove("run_program.sh")



