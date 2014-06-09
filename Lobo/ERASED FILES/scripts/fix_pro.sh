# this script removes CG and CD atoms from PRO residues which would 
# otherwise cause erroneous 3D coordinate shifts in LoopModelTest
grep 'CD  PRO' -v $1 | grep 'CG  PRO' -v > $2

