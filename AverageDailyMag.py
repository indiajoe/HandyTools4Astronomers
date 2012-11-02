#!/usr/bin/env python
# This is to average the magnitude of the same day in the second column for a single JD in 1st column
# Sorted JD is important in input.
# Usage: ./AverageDailyMag.py inputfile.txt
#...................................indiajoe@gmail.com
import numpy as np
import sys
if len(sys.argv) < 2 :
    print("Wrong input \n Usage: "+sys.argv[0]+" inputfile.txt" )
    exit(1)

inpfile=open(sys.argv[1],'r')
Workarray=[]  #Array to work on
for line in inpfile.readlines() :
    line=line.rstrip().split()
    JD=int(float(line[0]))
    if len(Workarray) == 0 or Workarray[-1][0] != JD :   #New JD night. so opening the night's entry
        Workarray.append([JD,[]])
    Workarray[-1][1].append(float(line[1]))
    
for i in range(len(Workarray)):
#    print(str(Workarray[i][0])+' '+str(Workarray[i][1])+' '+str(np.mean(Workarray[i][1]))+' '+str(np.std(Workarray[i][1])))
#    print(str(Workarray[i][0])+' '+str(round(np.mean(Workarray[i][1]),3))+' '+str(np.std(Workarray[i][1])))
    print(str(Workarray[i][0])+' '+str(round(np.mean(Workarray[i][1]),3)))
