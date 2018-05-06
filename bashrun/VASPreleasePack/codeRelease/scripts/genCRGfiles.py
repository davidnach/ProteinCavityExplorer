##
## VASP: Volumetric Analysis of Surface Properties
## Copyright (c) 2014 Brian Y. Chen
## All rights reserved.
##
## This file is part of VASP
## 
## VASP is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
## 
## VASP is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with VASP.  If not, see <http://www.gnu.org/licenses/>.  
##
## File: genCRGfiles.py
##       A CRG file generator
##
## Written by 
##       Brian Y. Chen <chen@lehigh.edu>
##
## WWW URL: http://cse.lehigh.edu/~chen/
## Email: chen@lehigh.edu
## Documentation can be found here: 
## http://www.ploscompbiol.org/article/info:doi/10.1371/journal.pcbi.1000881
##
##



import sys
import getopt
import string
import os
import re
import glob

###########################################################################
### Credit
###########################################################################
def credit():
	print "genCRGfiles.py, Brian Chen, Lehigh University, 2011."

###########################################################################
### Usage
###########################################################################
def usage():
	credit();
	print("========================================================="   );
	print(""                                                          );
	print("python genCRGfiles.py"        );
	print("   Running with no input returns this usage info."  );
	print(""                                                          );
	print("python genCRGfiles.py [crgFile] [pdbFile] [outputDir]"        );
	print("   Makes a new copy of the CRG file for each amino acid in pdbFile");
	print("   where that amino acid's charge is set to zero..");
	print(""                                                          );
	print("========================================================="   );

###########################################################################
###########################################################################
### File gobble
###########################################################################
###########################################################################
def gobble(fileName):
        #print("Gobbling: "+ fileName);
        
        inputFile = open(fileName, 'rt');
        inputLine = inputFile.readline();
        
        outputList = "".split();  #### list declaration
        
        while inputLine:
                outputList.append(inputLine);
                inputLine = inputFile.readline();

        inputFile.close();
        
#       for l in outputList:
#               print l;
        
        return outputList;


###########################################################################
###########################################################################
### outputCRGfile
###########################################################################
###########################################################################
def outputCRGfile( stumpList, chargeLineList, newCRG ):
	print "newCRG: " + newCRG;
	newCRGfile = open(newCRG, 'w');
	for line in stumpList:
		newCRGfile.write(line);
	
	for line in chargeLineList:
		newCRGfile.write(line);
	
	newCRGfile.close();
	
	
	


###########################################################################
###########################################################################
### fixPDB - reads a PDB file, eliminates garbage.  Takes NO chain or A chains only.
#NOTE - this does not fix the HETATM methionines.
###########################################################################
###########################################################################
##COLUMNS      DATA TYPE        FIELD      DEFINITION
##------------------------------------------------------
## 1 -  6      Record name      "ATOM    "
## 7 - 11      Integer          serial     Atom serial number.
##13 - 16      Atom             name       Atom name.
##17           Character        altLoc     Alternate location indicator.
##18 - 20      Residue name     resName    Residue name.
##22           Character        chainID    Chain identifier.
##23 - 26      Integer          resSeq     Residue sequence number.
##27           AChar            iCode      Code for insertion of residues.
##31 - 38      Real(8.3)        x          Orthogonal coordinates for X in 
##                                         Angstroms
##39 - 46      Real(8.3)        y          Orthogonal coordinates for Y in 
##                                         Angstroms
##47 - 54      Real(8.3)        z          Orthogonal coordinates for Z in 
##                                         Angstroms
##55 - 60      Real(6.2)        occupancy  Occupancy.
##61 - 66      Real(6.2)        tempFactor Temperature factor.
##77 - 78      LString(2)       element    Element symbol, right-justified.
##79 - 80      LString(2)       charge     Charge on the atom.
##

###########################################################################
###########################################################################
### genCRGfiles
###########################################################################
###########################################################################
def genCRGfiles(crgFileName, pdbFileName, outputDir):
	crgFileList = gobble(crgFileName);
	pdbFileList = gobble(pdbFileName);
	pdbHeader = pdbFileName.split("/")[-1].split(".")[0];

	##first filter the crg file for empty lines, creating a new list.
	stumpList = [];
	for line in crgFileList:
		if len( line.split() ) != 0:
			stumpList.append(line);

	##now begin looping through the PDB file.
	currentResNum = -9999
	chargeLineList = [];
	for line in pdbFileList:
		if not line.startswith("ATOM"):
			continue;

		resNum = int(line[22:26]);
		if( currentResNum == -9999 ):
			currentResNum = resNum;	
			
		if( resNum != currentResNum ):
			##output stuff
			
			counterString = '%(cnt)04i' % {'cnt' : currentResNum};
			outputCRGfile( stumpList, chargeLineList, outputDir+"/"+pdbHeader+"-"+counterString+".crg" );
			currentResNum = resNum;	
			chargeLineList = [];
			
		else:
			atom_str = "%(atom)-6s" % {'atom': line[12:16].split()[0] };
			resName_str = "%(rName)-3s" % {'rName': line[17:20].split()[0] };
			resNum_str = "%(rNum) 4i" % {'rNum': int( line[22:26].split()[0] ) };
			chain_str = line[21];
			charge_str = " 0.0";

			crgFileLine = atom_str + resName_str + resNum_str + chain_str + charge_str + "\n";
			chargeLineList.append(crgFileLine);

	##outputStuff.
	counterString = '%(cnt)04i' % {'cnt' : resNum};
	outputCRGfile( stumpList, chargeLineList, outputDir+"/"+pdbHeader+"-"+counterString+".crg" );






###########################################################################
###########################################################################
### Main Method
###########################################################################
###########################################################################
def main():
	
#	testDebug();
	
	#print sys.argv[1]
	#print len(sys.argv)
	#print "Hello, World!"

	################################################ Usage Statement
	if len(sys.argv) == 1:
		usage()
		raise SystemExit, 5

	#print("python genCRGfiles.py [crgFile] [pdbFile] [outputDir]"        );
	################################################ 
	if len(sys.argv) == 4:
		crgFileName = sys.argv[1];
		pdbFileName = sys.argv[2];
		outputDir = sys.argv[3];
		
		genCRGfiles(crgFileName, pdbFileName, outputDir);
		raise SystemExit, 5


if __name__ == "__main__":
	main()
	


