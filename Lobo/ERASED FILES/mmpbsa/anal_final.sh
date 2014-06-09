#! /bin/sh
grep "Total acces" report.out | awk '{print $5*0.02}' > sasa.nrg
grep "TOTAL PHISITE ENERGY" report.out > tmp.out
./printdispari.awk tmp.out > phinrg_w.dat
./printpari.awk tmp.out > phinrg_a.dat
grep "ENER>" report.out > ener.out
###### PER FARE I TEST 
grep "ENER INTERN>" report.out | awk '{print $3 +$4 +$5 +$6+$7}' > tmp_intern.nrg
grep "ENER EXTERN>" report.out  | awk '{print $3 +$4 * 0.5}' > tmp_extern.nrg
pr -m -t tmp_intern.nrg tmp_extern.nrg | awk '{print $1 + $2}' > mm.nrg
pr -m -t -w 240 phinrg_w.dat phinrg_a.dat | awk '{print ($6 - $12)*0.5}' > phinrg_w_m_a.nrg
pr -m -t -w 240 mm.nrg phinrg_w_m_a.nrg | awk '{print ($1 + $2)}' > mmpb.nrg
pr -m -t -w 240 mmpb.nrg sasa.nrg | awk '{print ($1 + $2)}' > mmpbsa.nrg

