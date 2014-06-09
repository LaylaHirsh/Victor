#!/bin/awk -f


BEGIN{ i=0 
        while((getline < "protein_gmx.pdb" ) >0) {
            nohis = 1
            line = $0
            resnum = $5
            f1 = $1
            f2 = $2
            f3 = $3
            f4 = $4
            f5 = $5
            f6 = $6
            f7 = $7
            f8 = $8
            f9 = $9
            f10 = $10
            
    while((getline  < "hsp.txt" ) >0){
            if($1==resnum && f1 == "ATOM") { 
        printf("ATOM %6i%5s HSP %5i    %8.3f%8.3f%8.3f %5.2f %5.2f\n",
                f2, f3, f5, f6, f7, f8, f9, f10)
                nohis = 0
                }
              }
                close("hsp.txt")

            if(nohis) 
            print line 
        }
}
