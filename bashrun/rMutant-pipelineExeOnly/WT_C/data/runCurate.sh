#!/bin/bash

# executable list
declare -a wtcArr=("../bin/curatePDBFileAllStepsExe" "../external/bndlst/bndlst" "../external/hbplus/hbplus" "../external/reduce/reduce_src/reduce")

# retrieve number of command line arguments
numArgs=$#

# check for python 3.5 or 3.4
python35=`which python3.5`
python34=`which python3.4`
if [ ${#python35} -gt 0 ]; then
    python=`which python3.5`
elif [ ${#python34} -gt 0 ]; then
     python=`which python3.4`   
else
    echo " need python 3.5 or 3.4"
    echo " quitting"
    exit
fi

# check that correct number of arguments have been received
if [ $numArgs -ne "2" ]; then
    echo " Expecting 2 arguments"
    echo " quitting"
    exit

elif [ $numArgs == "2" ]; then
    
    # check list of required programs/executables
    for i in "${wtcArr[@]}"
    do
        if [ ! -f ${i} ]; then
            #compile here from some other file call
            #move folder 
            echo "$i does not exist"
            exit
        fi
    done
    
    pdbID=${1}
    retainedChains=${2}

    # make the xml config file
    echo " Making config file"
    ${python} makeCurateXML.py ${pdbID} ${retainedChains}

    # curate the file
    echo " Curating"
    ../bin/curatePDBFileAllStepsExe ${pdbID}.${retainedChains}.curate.xml

fi
