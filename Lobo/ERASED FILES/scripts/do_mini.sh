#
# params: <x-ray-file> <basename-pdb&seq> <scwrl-basename> <AA-start> <AA-end>
#
time ~/Victor/Nazgul/lobo -i $1 -s $4 -e $5 -o $2.pdb --sol1 1000 --maxWrite 20 --cluster --simLimit 0.3 --scwrl $2.seq --verbose | tee $2.res.out
cd ./scripts/
rm etmp.*
./scwrl_all.pl $2.pdb $3 $1.seq $1
time ./opt_all_loop.sh $1 $3 $4 $5 | tee $2.opt.out
