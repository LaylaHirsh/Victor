#! /bin/csh -f
#
set n = `(./prot2res protein.pdb)` 
echo There are $n residues 
#
# copy the setup file in uhbd.uinp
#
sed "s/pdie pdie/pdie 1.0/" setup_header.protoinp | sed "s/sdie sdie/sdie 80.0/" \
     | sed "s/ions ions/ions 100.0/" | sed "s/grid_x/grid_w/" > uhbd_w.uinp 
sed "s/pdie pdie/pdie 1.0/" setup_header.protoinp | sed "s/sdie sdie/sdie 1.0/"  \
     | sed "s/ions ions/ions 0.0/" | sed "s/grid_x/grid_a/" > uhbd_a.uinp 
#
# cycle thorugh residues 
#
set i = 1
while( $i <= $n )
    sed "s/pdie pdie/pdie 1.0/" rescalc.protoinp | sed "s/sdie sdie/sdie 80.0/" | sed "s/ n_aa / $i /" \
      | sed "s/ions ions/ions 100.0/" | sed "s/grid_x/grid_w/" | sed "s/n_aa.pdb/aa$i.pdb/" >> uhbd_w.uinp 
    sed "s/pdie pdie/pdie 1.0/" rescalc.protoinp | sed "s/sdie sdie/sdie 1.0/" | sed "s/ n_aa / $i /" \
      | sed "s/ions ions/ions 0.0/" | sed "s/grid_x/grid_a/" | sed "s/n_aa.pdb/aa$i.pdb/" >> uhbd_a.uinp 
    @ i = $i + 1
end
    cat close.protoinp >> uhbd_w.uinp
    sed "s/print sasa/\!print sasa/" close.protoinp | grep -v saforc >> uhbd_a.uinp
