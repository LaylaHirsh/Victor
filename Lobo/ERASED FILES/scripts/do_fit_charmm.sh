
#echo "Usage: do_fit_charmm.sh <filename-xray> <filename-model>
# <aa#-first-in-loop> <aa#-last-in-loop>"

# fix OT1/2 to O/OXT compatibility between CHARMM and rest of the world
sed -e "s/OT2/OXT/g
        s/OT1/O  /g" < $2 > temp2.pdb

let st=$3-1
let en=$4+1

# customize input file for ProFit with parameters
sed -e "s/zzz/$1/g
        s/aaa/$st/g
        s/bbb/$en/g
        s/ccc/$3/g
        s/ddd/$4/g" < fit1.in2 > fit1.this2
sed -e "s/zzz/$1/g
        s/aaa/$st/g
        s/bbb/$en/g
        s/ccc/$3/g
        s/ddd/$4/g" < fit1.in > fit1.this

# calculate RMSD with ProFit
echo ">>>>> $1 vs. $2 MINIMIZED"
echo ">>>  LOCAL:"
profit < fit1.this2 | grep "RMS" | tail -1
echo ">>>  GLOBAL:"
profit < fit1.this | grep "RMS" | tail -1
