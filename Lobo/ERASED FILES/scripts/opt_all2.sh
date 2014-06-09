#! /bin/csh

echo '------------------------------------------------------------------'
foreach in ( $2* )
echo "Energy2 for            $in "
./opt_indel2.sh $1 $in $3 $4 $in
end
