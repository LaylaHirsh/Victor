#! /bin/csh

echo '------------------------------------------------------------------'
foreach in ( $2* )
echo "Energy for            $in "
./opt_indel.sh $1 $in $3 $4 $in $5
end
