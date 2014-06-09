#! /bin/bash

# Runs a local CHARMM optimization over a single result from LoopModelTest
# Usage:  opt_loop.sh <filename-xray> <filename-model> <aa#-first-in-loop>
# <aa#-last-in-loop> <filename-for-temporary-without-.pdb>

scwrl -i $2 -s $6 -o $5_tmp.pdb > /dev/null
grep "OXT" $2 >> $5_tmp.pdb

./do_min.sh $5_tmp $3 $4 $1
echo '------------------------------------------------------------------'
