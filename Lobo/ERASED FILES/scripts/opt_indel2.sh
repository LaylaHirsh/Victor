#! /bin/bash

# Runs a local CHARMM optimization over a single result from LoopModelTest
# Usage:  opt_loop.sh <filename-xray> <filename-model> <aa#-first-in-loop>
# <aa#-last-in-loop> <filename-for-temporary-without-.pdb>

scwrl -i $2 -o $5_tmp2.pdb > /dev/null
grep "OXT" $2 >> $5_tmp2.pdb

./do_min2.sh $5_tmp2 $3 $4
echo '------------------------------------------------------------------'
