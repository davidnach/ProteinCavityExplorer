import requests
import sys
import os.path
pid = sys.argv[1]
# assumes that it will find at least one cavity
insurfFile= open(pid + ".txt", "r")
outsurfFile = open(pid + ".surfReport", "w")
cavityfound = 0
lastresidue = 0
cavitycount = 0

for line in insurfFile:
	if line [0:6] == "CAVITY":
		if cavityfound == 1:
			outsurfFile.write("}\n")
		cavityfound = 1
#area currently is a dummy number until we figure out how to calculate it
		outsurfFile.write("Pocket " + str(cavitycount) +" (area: 1.0): {")
		cavitycount += 1
		lastresidue = 0
	elif cavityfound == 1:
		if lastresidue != line.split()[5]:
			outsurfFile.write(line.split()[5]+",")
			lastresidue = line.split()[5]
outsurfFile.write("}\n")

outsurfFile.close()

