#!/usr/bin/perl -w

if (@ARGV < 1)
{
    print "Usage:\n";
    print "doscwrl.pl <datafile> \n";

    die "Not enough arguments.\n";
}

open (FILE1, $ARGV[0]);

while ($data = <FILE1>)
{
    chomp($data);
    
    print ">>> $data\n";

    $tmp = substr($data,0,4) . "_scwrl." . substr($data,9,12) . ".seq";
    $tmp2 = $data . ".pdb";
    $tmp3 = $data . ".scwrl.pdb";
    
    for ($i = 0; $i < 50; $i++)
    {
	if ($i < 10)
	{
	    $ext = "00";
	}
	elsif ($i < 100)
	{
	    $ext = "0";
	}
	else
	{	 
	    $ext = "";
	}

	$ext = $ext . $i;

	print "$data $ext\n";

	`scwrl -i $tmp2.$ext -o $tmp3.$ext -s $tmp`;
#	print "scwrl -i $tmp2.$ext -o $tmp3.$ext -s $tmp \n";

    }
}


