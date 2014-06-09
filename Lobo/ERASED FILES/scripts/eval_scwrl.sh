#! /bin/csh

foreach i ( $2.* )
    ls $i > tmp.str
    cut -d. -f3 tmp.str
    rm tmp.str
    scwrl -i $i -o temp.pdb > /dev/null
    grep "OXT" $i >> temp.pdb
    ./LoopModelEval -i $1 -m temp.pdb -s $3 -e $4
end

