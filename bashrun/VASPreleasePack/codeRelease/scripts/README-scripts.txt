==============================================================================
This document describes data and source code used as input for the following 
publication:

VASP-E: Specificity Annotation with a Volumetric Analysis of
Electrostatic Isopotentials

by Brian Y. Chen, 
Dept. Computer Science and Engineering, Lehigh University

Please direct any questions to <chen@cse.lehigh.edu>
==============================================================================

License:

The source code provided here is free software: you can redistribute it and/or 
modify it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or (at your 
option) any later version.

The source code provided is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License 
for more details.

A copy of the GNU General Public License is available in each software 
directory described below.  If you have not received this license, see 
<http://www.gnu.org/licenses/>.  



Compilation:

Scripts in this folder are interpreted; they do not require compilation.


Overall Usage:

Delphi uses .CRG files to specify a range of biophysical settings used when
solving the Poisson-Boltzmann equation.  CRG files are used to nullify the 
contribution of an individual amino acid; an individual CRG file is generated
for every variation of the given protein with a different nullified residue.

genCRGfiles.py generates CRG files necessary for nullifying every individual
amino acid once.  Once the CRG files are generated the the delphi_customCRG.pl
script is used to generate GRD files that reflect the nullification in the
electrostatic potential field.


genCRGfiles.py:

This script generates CRG files using the following commandline:

python genCRGfiles.py [crgFile] [pdbFile] [outputDir]

[crgFile] is a path to a default CRG file with default charges.  The 
amber.crg default charge file is provided here for this purpose.  [pdbFile]
is the pdb file for which you want to generate CRG files.  [outputDir]
is a path where you would like the output files generated.  For example:

mkdir test;
python ./genCRGfiles.py charmm.crg ./1a0j.pdb ./test

will generate CRG files for 1a0j.pdb (provided as a demonstrate) in the test
folder.

* You should only perform this nullification on pdb files with a single chain.



delphi_customCRG.pl:

This script calls DelPhi to solve the Poisson-Boltzmann equation for a given
PDB file and CRG file.  The commandline for running it is:

perl delphi_customCRG.pl [pdbFile] [crgFile]

*  This script requires setup for your system before it can be run.  In
particular, it must know the location of delphi, in order to call it.  To set
this path, change the following line to reflect the location of your delphi
executable:

my $DELPHI = '/infolab/exchange/software/delphi_new/executable/delphi';

Also, change the following line to reflect the location of the SIZ file, which
specifies the Van der Waals radii for all atoms in use:

my $SIZ = '/infolab/exchange/software/delphi_new/parameters/amber.siz

A copy of this file is provided for you, but you may use your own as well.

* it is crucial to note that Delphi generates temporary files in the current
working directory, and it can be confused if two instances of Delphi run at the
same time - they will overright each other's temp files and receive
contradictory data.  For this reason it is strongly recommended that you run
each instance of delphi from an individual folder from which no other instance
will be run.







