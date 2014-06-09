#! /bin/csh

foreach i ( $2.* )
    ls $i > tmp.str
    cut -d. -f3 tmp.str
    ./LoopModelEval -i $1 -m $i -s $3 -e $4
end

