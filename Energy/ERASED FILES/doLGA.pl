#!/usr/bin/perl -w

if (@ARGV < 1)
{
    print "Usage:\n";
    print "doLGA.pl <datafile>\n";

    die "Not enough arguments.\n";
}

open (FILE1, $ARGV[0]);

while ($data = <FILE1>)
{
    chomp($data);
    
    $out = "res-" . "$data" . ".rmsd";
    `cat /dev/null > $out`;

    print "-" x 50 . "\n$data\n";

    $base = substr($data, 0, 5) . "TS\*";

    @list = `ls $base`;
    chomp(@list);

    for ($i = 0; $i <= $#list; $i++)
    {
	print "--> $list[$i]\n";
	`echo $list[$i] >> $out`;
	`runlga.mol_mol.pl $list[$i] $data -3 -sia -d:5 -o1`;
	`cat RESULTS/$list[$i].$data.res | grep 'SUMMARY' >> $out`;
    }


}

close(FILE1);

