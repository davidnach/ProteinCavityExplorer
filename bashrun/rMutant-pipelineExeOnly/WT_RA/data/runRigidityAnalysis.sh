#!/bin/bash

numArgs=$#
pdbID=""
chainID=""

# check for that we are using python3 
# if no python3, which python3 will return nothing
corrPython=`which python3`
sizeofcorrPython=${#corrPython}

if [ $numArgs -ne "2" ]; then

    echo " Expecting 2 arguments"
    echo " quitting"
    exit
 
elif [ $numArgs == "2" ]; then

    # if corrPython returns nothing, python3 file doesn't exist
    if [ $sizeofcorrPython -lt 1 ]; then

        echo " Python3 does not exist"
        exit

    else

        echo " Python3 exists"
    fi

    pdbID=${1}
    chainID=${2}

    echo " Received arguments $pdbID and $chainID"

    # make config file
    echo " Making config file for ${pdbID}"
    /usr/bin/python3 makeRigidityAnalysisXML.py ${pdbID} ${chainID}

    # performing rigidity analysis
    echo " Performing Rigidity Analysis on ${pdbID}"
    ../bin/kinariMolWithBondInputExe ${pdbID}.${chainID}.ra.xml

fi


