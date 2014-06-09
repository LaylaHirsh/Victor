#! /usr/bin/perl

if (@ARGV < 4)
{
    print "Usage:\n";
    print "scwrl_all.pl <basename-model> <basename-scwrl> <scwrl-seqfile> <xray-filename\n";
    
    die "Not enough arguments.\n";
}

$base = $ARGV[0];
$scwrl = $ARGV[1];
$seq = $ARGV[2];
$oxtfile = $ARGV[3];

$id = "$$"; # process id

@list = `ls $base*`;

for ($i = 0; $i <= $#list; $i++)
{
    chomp($list[$i]);
    print "$list[$i]\n";
    $out = $scwrl . ".";
    if ($i < 10)
    {
	$out .= "00";
    }
    else
    {
	$out .= "0";
    }
    $out .= $i;
    `scwrl -i $list[$i] -o $out -s $seq`;
    `grep "TER" -v $out > scw.$id`;
    `grep "OXT" $oxtfile >> scw.$id`;
    `grep "TER" $oxtfile >> scw.$id`;
    `mv scw.$id $out`;
}
