#! /bin/csh
##### Build first subunit ############
cp $1 protein.pdb
echo 1 | pdb2gmx_d -f protein.pdb -o protein_gmx.pdb -p protein.top -ignh 
echo 1 | pdb2gmx_d -f protein.pdb -o protein.gro -p protein.top  -ignh 
rm *.itp.* *.top.* *.gro.* \#*.pdb.*\#
grep "HISA" protein.gro | grep HD1 | cut -b 1-5 > hsd.txt
grep "HISB" protein.gro | grep HE2 | cut -b 1-5 > hse.txt
grep "HISH" protein.gro | grep HE2 | cut -b 1-5 > hsp.txt
if (! -z hsp.txt) then
awk -f get_hsp.awk > tmp.pdb
mv tmp.pdb protein_gmx.pdb
endif
if (! -z hsd.txt) then
awk -f get_hsd.awk > tmp.pdb
mv tmp.pdb protein_gmx.pdb
endif
if (! -z hse.txt) then
awk -f get_hse.awk > tmp.pdb
mv tmp.pdb protein_gmx.pdb
endif
grep SG protein_gmx.pdb > sg.txt
./ssbonds sg.txt disu.str
set n = `(wc -l protein_gmx.pdb | awk '{print $1 - 6}')`
echo "NATOM  $n" > prot.pdb
grep ^ATOM protein_gmx.pdb | sed -f pdb2cmm.sed | cut -b 1-66 | awk '{printf("%66s%10s\n", $0,"PROT")}' >> prot.pdb
echo "TER" >> prot.pdb
charmm < build.inp > build.out 
cp protout.pdb $1.out
exit 0
