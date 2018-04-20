#!/bin/bash

rm -rf *~

cd WT_C/data
./cleanStage.sh
cd ../..

cd WT_RA/data
./cleanStage.sh
cd ../..
