#!/bin/bash

#File : invokePipeline_WT_allChains.sh
#Last Revision : 11 May 2016
#Authors : Rebecca Hsieh, Filip Jagodzinski, WWU

#This file contains the command-line programs/arguments that you'd
#issue for invoking curation and rigidity analysis for the wild type
#for any PDB

# retrieve the number of arguments received
numArgs=$#

# check that exactly 1 argument has been received
if [ $numArgs -ne "1" ]; then
    echo " Expecting exactly 1 argument."
    echo " Have received ${numArgs}"
    echo " Quitting"
    exit
fi

# retrieve the 1 argument and save into variables
pdbID=${1}
chainID="all"

# set variable to remember the top-most directory
topDir=`pwd`

# curate of WT
cd WT_C/data
./runCurate.sh ${pdbID} ${chainID}

printf '\n'
echo " WILD TYPE CURATION COMPLETED "
printf '\n'

# perform rigidity analysis of curated WT
cd ${topDir}
cd WT_RA/data

`cp -r ../../WT_C/data/${pdbID}.${chainID}.cur.out .`
./runRigidityAnalysis.sh ${pdbID} ${chainID}

printf '\n'
echo " WILD TYPE RIGIDITY ANALYSIS COMPLETED "
printf '\n'


cd ${topDir}
printf '\n'
echo " CONFIRMING RUN-TIME SUCCESS."
echo " CHECKING FOR EXPECTED FILES."

allFilesExist=true
listOfFiles=( WT_C/data/${pdbID}.${chainID}.cur.out/${pdbID}.${chainID}.processed.pdb.knr
	      WT_RA/data/${pdbID}.${chainID}.ra.out/${pdbID}.${chainID}.processed.pdb/user/${pdbID}.${chainID}.processed.pdb_postPG_MetricsBBH.txt
	    )

for file in "${listOfFiles[@]}"
do
    if [ -f $file ]
    then
	echo " Confirmed. file $file exists"
    else
	allFilesExist=false
	echo " Error! file $file does not exist"
    fi
done

if [ ${allFilesExist} = false ]
then
    echo " Programs were not successful!"
else
    echo " Programs ran successfully."
    echo " First 28 lines of the rigidity data"
    head -28 WT_RA/data/${pdbID}.${chainID}.ra.out/${pdbID}.${chainID}.processed.pdb/user/${pdbID}.${chainID}.processed.pdb_postPG_MetricsBBH.txt
    
fi

