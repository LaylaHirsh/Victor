# !/bin/csh

# fix segment id for CHARMM input
cat /dev/null > tmp1.pdb
./add_tag.pl $1.pdb TMP2 > tmp2.pdb

# now build input files for CHARMM minimization
charmm pdbfile=tmp2 toppar=/usr/local/src/c27b3/toppar/ < build.inp > build.out

let st=$2-1
let en=$3+1

# minimization proper
charmm pdbfile=tmp2 indstart=$st indend=$en < min-loop.inp > min.out

grep '^EEEEE' min.out
echo 'RAPDF + SOLV'

cat /dev/null > etmp3.pdb
sed -e "s/OT2/OXT/g
        s/OT1/O  /g
        s/HSD/HIS/g
        s/CD  ILE/CD1 ILE/g" < tmp2_min.pdb > etmp3.pdb

./LoopModelEval -i etmp3.pdb -m etmp3.pdb -s $st -e $en
echo '-=-=-=-=-=-=-=-'
