#!/usr/bin/env python
# This is to generate the column of two data matching the nearest UT of observation
# Usage: ./ColorCalc.py Ref1stfile.txt 2ndfile.txt
#...................................indiajoe@gmail.com
import numpy as np
import operator
import sys
if len(sys.argv) < 3 :
    print("Wrong input \n Usage: "+sys.argv[0]+" Ref1stfile.txt 2ndfile.txt" )
    exit(1)

REFfile=open(sys.argv[1],'r')
File2match=open(sys.argv[2],'r')
RefArray=[inp.rstrip().split() for inp in REFfile.readlines()]
InpArray=[inp.rstrip().split() for inp in File2match.readlines()]
REFfile.close()
File2match.close()
inpJD=[eval(inp[0]) for inp in InpArray]
len_of2ndfile=len(inpJD) # List of the JDs in 2nd file
for i in RefArray :
    refJD=np.ones(len_of2ndfile)*eval(i[0])  #Creating an array of RefJD value
    min_index, min_value = min(enumerate(abs(inpJD-refJD)), key=operator.itemgetter(1))
    print(' '.join([i[0],i[1],InpArray[min_index][0],InpArray[min_index][1],str(eval(i[1])-eval(InpArray[min_index][1]))]))
