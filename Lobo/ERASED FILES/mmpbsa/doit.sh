#!/bin/csh
date
@ j = 0
@ n = 10
####### Cleaning up first... ##################
echo "Removing files from previous runs"
set old_file_list = ( sasa.nrg \
phinrg_w_m_a.nrg \
mm.nrg \
rmsd.dat\
report.out ) 
foreach file ($old_file_list[*])
if (-e $file) then
rm $file
else
echo "$file not found"
endif
end
####### Cycling over proteins ##################
while ($j <= $n)
  echo "STO FACENDO LA PROTEINA:" protein$j.pdb >> report.out
  echo "STO FACENDO LA PROTEINA:" protein$j.pdb 
  ./build.sh protein$j.pdb
  echo "protein built"
  ./calc_cont.sh
  echo "protein calc'd" 
  ./analyse.sh 
  echo "output analysed" 
@ j = $j + 1
end
