# import system calls
import sys

# get PDBID
pdbID = sys.argv[1]
chainID = sys.argv[2]

# generate output name
if (chainID == "all"):
    fileOutName = pdbID + ".all.curate.xml"
else:
    fileOutName = pdbID + "." + chainID + ".curate.xml"

# file input and output
if (chainID == "all") :
    fileIn = open("curateTemplateAllChains.txt", "r")    
else:
    fileIn = open("curateTemplateSpecificChain.txt", "r")    

fileOut = open(fileOutName, "w")

# parse input file, pipe to output
for line in fileIn :

    if "{PDBID}" in line :
        line = line.replace("{PDBID}", pdbID)

    if "{CHAINID}" in line :
        line = line.replace("{CHAINID}", chainID)

    fileOut.write(line)

# close input and output files
fileIn.close()
fileOut.close()
