#!/usr/bin/env python2
# Creates a directory for the PID, if it does not exists already
# imports resources?
# Creates a directory for the experiment ID within the corresponding PID directory
# Experiment ID is hexidecimal

import os
import sys

bigPrime = 2305843009213693951
def nextId(oldGen):
    return ((oldGen+672186240682106308) % bigPrime) # mod the 9th Mersenne prime 2**61-1
 
# first experiment id
os.chdir(os.path.dirname(os.path.realpath(__file__))+'/experimentData')
if(not os.path.isfile('generator.txt')):
    # add file to hold experiment id generator
    initGenFile = open('generator.txt','w')
    initGenFile.write('1231495648941321981')
    initGenFile.close()

# Assuming PIDs are always 4 characters
inPID = sys.argv[1][0:4] #raw_input('Input PID: ')[0:4]
pidPath=inPID+'/'

# check user input #

# create PID directory if it does not exist
if not os.path.exists(pidPath):
    os.mkdir(pidPath)

# grab last id used to generate new id
# (in integer form)
oldGenFile = open('generator.txt','r')
expId=int(oldGenFile.read())
oldGenFile.close()

# drop in generator function
expId=nextId(expId)

# replace generator with new generator
genFile = open('generator.txt','w')
genFile.write(str(expId))
genFile.close()

# use hexidecimal for shorter user ids
hexExpId = hex(expId)[2:] # drop the 0x part
# 0 pad hexExpId's so they look nice
hexExpId='0'*(len(hex(bigPrime)[2:])-len(hexExpId))+hexExpId
os.mkdir(pidPath+hexExpId)
os.chdir(pidPath+hexExpId)
print hexExpId
# Now can mess with experiment directory 


    


