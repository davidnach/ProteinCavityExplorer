#!/bin/bash
# usage: runrmutantandvasp.sh [pid] [expId]

clear

echo "The script starts now."
pdbID=${1}
expID=$2
echo "running pdbfetch with ${pdbID}"
python downloadopdb.py ${pdbID}
echo "copying pdb file for VASP later."
expPath=../experimentData/${pdbID}/${expID}
cp -r pdbs/${pdbID}.pdb VASPreleasePack/codeRelease/surfProcessingRelease/debug
cp -r pdbs/${pdbID}.pdb ${expPath}
echo "running rmutant with ${pdbID} and chain A."
cd rMutant-pipelineExeOnly
./invokePipeline_WT.sh ${pdbID} A
cd ..
cp -r rMutant-pipelineExeOnly/WT_RA/data/${pdbID}.A.ra.out/${pdbID}.A.processed.pdb.knr ${expPath}
cp -r rMutant-pipelineExeOnly/WT_RA/data/${pdbID}.A.ra.xml ${expPath}
cd VASPreleasePack/codeRelease/surfProcessingRelease/debug
echo "running VASP with ${pdbID} into surf file with radius 5.0 and resolution .5."
./surfProcessing -surfProbeGen ${pdbID}.pdb ${pdbID}_Envelope.SURF 5.0 .5
# Try r=1.4 Angstrom probe
./surfProcessing -surfProbeGen ${pdbID}.pdb ${pdbID}.SURF 1.4 .5
cd ../../../..
cp -r VASPreleasePack/codeRelease/surfProcessingRelease/debug/${pdbID}.SURF ${expPath} 
cp -r VASPreleasePack/codeRelease/surfProcessingRelease/debug/${pdbID}_Envelope.SURF ${expPath}

# clean rMutant
cd  rMutant-pipelineExeOnly
./cleanPipeline.sh
cd ..


#### surfracer does not detect the correct type of cavity
#### it detects internal, solvent inaccessible cavities

#(echo 1;echo ${pdbID}.pdb;echo 1.4;echo 1;echo ; echo;) | ./surfrace5_0_linux_64bit >> /dev/null
#
#if [ "$?" -eq "139" ]; then surfFailed=true;
#else
#    # copy stuff
#    cp ${pid}.txt -r ${expPath}
#    cp ${pid}_residue.txt -r ${expPath}
#    cp -r results.txt ${expPath}
#    rm results.txt;
#    # do pracpa
#fi









