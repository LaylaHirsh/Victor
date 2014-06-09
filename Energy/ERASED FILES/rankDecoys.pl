#!/usr/bin/perl -w

## -*- Perl -*-----------------------------------------------------------------
##  $Id: rankDecoys.pl,v 1.2 2007-12-17 09:23:12 biocomp Exp $
##
##  Author:             Silvio Tosatto
##
##  Project Name:       FRST
##
##  Date:               04/05
##
##  Description:
##    This script generates rank1, z-score, c.c., f.e., logP_B1, logP_B10  
##    output for decoy set data.
##
## ---------------------------------------------------------------------------
use Getopt::Std;

$| = 1;
my ($script) = ($0 =~ m|([^/]*)$|);

$Use = "$script - Decoy rank data generator
This script generates rank1, z-score, c.c., f.e., logP_B1, logP_B10  output 
for decoy set data.
 Usage: 
\t $script -i <filelist> 
\n";

#
# START of main program
#

getopts("hi:", \%args) or die "Aborting. (-h for help)\n"; 

if ($args{h}) 
{
   die "$Use";
}

if (!$args{i})
{
    die "No input filelist specified. See source code for info.\n";
}

for ($i = 1; $i <= 5; $i++)
{
    `cat /dev/null > zscore.frst.r$i`;
}

$count = 0;
$prev = "zzzz";
open (FILE1, $args{i});
while ($data = <FILE1>)
{
    chomp($data);

#    if ((substr($data, 0,5) eq "_") || (substr($data, 0,5) eq "-"))
#    {
#	$dir = substr($data, 0,6);
#    }
#    else
#    {
#	$dir = substr($data, 0,4);
#    }

    $dir = $data;
    print ";;;$dir;;;\n";

    if ($dir ne $prev)   # get rid of duplicates
    {
#	print "$data;;;\t;;;$dir;;;\n";
	$count++;
	$native = "$dir.pdb";

	`rm -f $dir/$dir.frst`;

        # generate FRST for all PDB -- assuming this has been already done!

#	`ls $dir/ | grep '.pdb' > $dir/filelist`;

	`awk '{gsub("is", "", \$0); gsub("cRMSD between $native and ","",\$0); print \$0}' $dir/rmsds > $dir/$dir.rmsd`;

	for ($i = 1; $i <= 5; $i++)
	{
	    $j = $i + 1;
	    `sort -n -k$j $dir/$dir.energy.frst > $dir/$dir.frst.s$i`;
	    `energy2zscore -s -i $dir/$dir.frst.s$i -r $dir/$dir.rmsd -n $native -c $i > $dir/$dir.frst.r$i`;
	    `cat  $dir/$dir.frst.r$i >> zscore.frst.r$i`;
	}

	$prev = $dir;
    }
}

print "Ranking data generated for a total of $count directories.\n";
    
