#! /bin/csh
########################################
# 1 - Remove files from previous runs
########################################
cat protout.pdb | grep -v END | grep -v REMARK > uhbd/protein.pdb
cd uhbd
echo "Removing files from previous runs"
rm aa*.pdb
set old_file_list =  ( grid_w.grd \
grid_a.grd \
sasa.out \
rmsd.out \
phinrg_a.out \
phinrg_w.out  )
foreach file ($old_file_list[*])
if (-e $file) then
rm $file 
else
echo "$file not found"
endif
end
############################################
# 2 - get the structure and calculate 
############################################
./prepuhbd.csh
set n = `(./prot2res protein.pdb)`
uhbd < uhbd_w.uinp > uhbd_w.uout
uhbd < uhbd_a.uinp > uhbd_a.uout
cd ..
exit 0
