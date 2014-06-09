#! /bin/bash

# Runs a local CHARMM optimization over a single result from LoopModelTest
# Usage:  opt_loop.sh <filename-xray> <filename-model> <aa#-first-in-loop>
# <aa#-last-in-loop>

echo '---> ' $2
echo '-------------------------------'
./do_fit.pl $1 $2 $3 $4 
echo '-------------------------------'
./do_min.pl $2 $3 $4
echo '-------------------------------'
./do_fit.pl $1 etmp.pdb $3 $4
echo '-------------------------------'
#rm etmp.pdb
