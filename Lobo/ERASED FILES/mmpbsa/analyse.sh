#! /bin/sh
grep "Total acces" uhbd/uhbd_w.uout > sasa.out
grep "TOTAL PHISITE ENERGY" uhbd/uhbd_w.uout > phinrg_w.out
grep "TOTAL PHISITE ENERGY" uhbd/uhbd_a.uout > phinrg_a.out
grep "ENER>" build.out | tail -1 > ener.out
grep "ENER INTERN>" build.out | tail -1 > ener_intern.out
grep "ENER EXTERN>" build.out | tail -1 > ener_extern.out
grep "ENER CONSTR>" build.out | tail -1 > ener_constr.out
grep "ENER POLAR>" build.out | tail -1 > polar.out
pr -m -t -w 240 phinrg_w.out phinrg_a.out | awk '{print ($6 - $12)}' >> phinrg_w_m_a.nrg
awk '{print $3}' ener.out >> mm.nrg
echo $1 >> report.out
cat phinrg_w.out phinrg_a.out ener.out ener_intern.out ener_extern.out ener_constr.out sasa.out >> report.out
