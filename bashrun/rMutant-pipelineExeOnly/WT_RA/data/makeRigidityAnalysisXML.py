# import system calls
import sys

# get PDBID and chainID
pdbID = sys.argv[1]
chainID = sys.argv[2]

# generate output name
fileOutName = pdbID + "." + chainID + ".ra.xml"

# file input and output
fileIn = open("rigidityAnalysisTemplate.txt", "r")
fileOut = open(fileOutName, "w")

# parse input file
for line in fileIn :

    if "{PDBID}" in line :
        line = line.replace("{PDBID}", pdbID)

    if "{CHAINID}" in line :
        line = line.replace("{CHAINID}", chainID)

    fileOut.write(line)

# close input and output files
fileIn.close()
fileOut.close()
