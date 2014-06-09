#! /bin/bash 
#echo "Usage: do_fit.sh <filename-xray> <filename-model> <aa#-first-in-loop> 
# <aa#-last-in-loop>"

# build missing side chains with SCWRL and fix OXT bug in output 
scwrl -i $2 -o temp.pdb > /dev/null
grep "OXT" $2 >> temp.pdb

let st=$3-1
let en=$4+1

# customize input file for ProFit with parameters
sed -e "s/zzz/$1/g
        s/aaa/$st/g
        s/bbb/$en/g
        s/ccc/$3/g
        s/ddd/$4/g" < fit.in2 > fit.this2
sed -e "s/zzz/$1/g
        s/aaa/$st/g
        s/bbb/$en/g
        s/ccc/$3/g
        s/ddd/$4/g" < fit.in > fit.this

# calculate RMSD with ProFit
echo ">>>>> $1 vs. $2 UNOPTIMIZED"
echo ">>>  LOCAL:" 
profit < fit.this2 | grep "RMS" | tail -1
echo ">>>  GLOBAL:" 
profit < fit.this | grep "RMS" | tail -1
rm fit.this2 fit.this
