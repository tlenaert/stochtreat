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


treattime=2
if len(sys.argv[:]) > 1 :
    treattime=sys.argv[1]
folder_name_prefix="initialresp_"+str(treattime)
no_patients=1000
no_runs=100
treattest=1
output="nolsctime.diagtime.initresponse.fullburden.nooverview"

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

    for id in range(no_runs):
        parameters= "id="+str(id)+",patients="+str(no_patients)+",output="+str(output)+",treattest="+str(treattest)+",treattime="+str(treattime)
        jobname= parameters
        walltime_parameter="01:30:00"#"walltime="+
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



