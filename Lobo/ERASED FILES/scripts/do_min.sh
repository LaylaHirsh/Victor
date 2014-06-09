# !/bin/csh

# fix segment id for CHARMM input
cat /dev/null > tmp1.pdb
./add_tag.pl $1.pdb TMP1 > tmp1.pdb

# now build input files for CHARMM minimization
charmm pdbfile=tmp1 toppar=/usr/local/src/c27b3/toppar/ < build.inp > build.out

let st=$2-1
let en=$3+1

# minimization proper
charmm pdbfile=tmp1 indstart=$st indend=$en < min-loop.inp > min.out

grep '^EEEEE' min.out
echo '------===========------'
echo 'RAPDF + SOLV'

cat /dev/null > etmp2.pdb
sed -e "s/OT2/OXT/g
        s/OT1/O  /g
        s/HSD/HIS/g
        s/CD  ILE/CD1 ILE/g" < tmp1_min.pdb > etmp2.pdb

#./LoopModelEval -i etmp2.pdb -m etmp2.pdb -s $st -e $en
pdb3energy -i etmp2.pdb
echo '-=-=-=-=-=-=-=-'
