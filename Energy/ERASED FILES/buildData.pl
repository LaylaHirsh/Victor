#!/usr/bin/perl -w

if (@ARGV < 1)
{
    print "Usage:\n";
    print "buildData.pl <datafile> <outputfile>\n";

    die "Not enough arguments.\n";
}

open (FILE1, $ARGV[0]);

$outfile = $ARGV[1];

`cat /dev/null > $outfile`;

while ($data = <FILE1>)
{
    chomp($data);
    
    $rmsdFile = "res-" . "$data" . ".rmsd";
    $energyFile = "res-" . substr($data,0,5) . ".energy";

    print "-" x 50 . "\n$data\n";

    open(FILE2, $energyFile);
    @energyData = <FILE2>;
    chomp(@energyData);
    close(FILE2);
    
    open(FILE3, $rmsdFile);
    @rmsdData = <FILE3>;
    chomp(@rmsdData);
    close(FILE3);

    open(OUT, ">>$outfile");

    print OUT "DECOY\n";

    for ($i = 0; $i <= $#energyData; $i++)
    {
	$id = substr($energyData[$i], 0, 13);
	if (substr($id, 12, 1) eq "\t")
	{
	    $id = substr($id, 0, 12);
	}
	$data = substr($energyData[$i], 13, 99);
	chomp($data);
	
#	print "--> ;;$id;;$data\n";
	
	# search for corresponding GDT_TS score, error if not found:

	for ($j = 0; $j <= $#rmsdData; $j++)
	{
	    $idRMSD = substr($rmsdData[$j], 0, 13); 

#	    print "#-#-#-> ;;$idRMSD;;\n";

	    if ($id eq $idRMSD)
	    {
		if (substr($rmsdData[$j+1], 0, 12) ne "SUMMARY(GDT)")
		{
		    print "Error: No RMSD data found for energy value of id = $id\n";
		}
		else
		{
		    $dataRMSD = substr($rmsdData[$j+1], 48, 6);
		    chomp($dataRMSD);
		    print OUT "$dataRMSD\t$data\n";
		};
		break;
	    }
	}
    }
}

close(FILE1);
close(OUT);

#
# END of main program
#
