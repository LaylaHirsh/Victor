# !/bin/csh

# fix segment id for CHARMM input
cat /dev/null > tmp1.pdb
./add_tag.pl $1.pdb TMP1 > tmp1.pdb

# now build input files for CHARMM minimization
charmm pdbfile=tmp1 toppar=/usr/local/src/c27b3/toppar/ < build.inp > build.out

let st=$2-1
let en=$3+1
let sttwo=$4-1
let entwo=$5+1

# minimization proper
charmm pdbfile=tmp1 indstart=$st indend=$en indstarttwo=$sttwo indendtwo=$entwo < min-2loops.inp > min.out

grep '^EEEEE' min.out
echo 'RAPDF + SOLV'

cat /dev/null > two_loops.pdb
sed -e "s/OT2/OXT/g
        s/OT1/O  /g
        s/HSD/HIS/g
        s/CD  ILE/CD1 ILE/g" < tmp1_min.pdb > two_loops.pdb

echo '-=-=-=-=-=-=-=-'
