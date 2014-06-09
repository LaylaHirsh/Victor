#! /usr/bin/perl

if (@ARGV < 4)
{
    print "Usage:\n";
    print "do_fit.pl <filename-xray> <filename-model> <1stAA-loop> <lastAA-loop> [<id>]\n";
    
    die "Not enough arguments.\n";
}

$st = $ARGV[2] - 1;
$en = $ARGV[3] + 1;

if (@ARGV < 5)
{
    $id = "$$"; # this is the current process id!!!
}
else
{
    $id = $ARGV[4];
};

# customize input files for ProFit with parameters:

open (FIT1, ">fit.$id");
print FIT1 "quiet
reference $ARGV[0]
mobile $ARGV[1]
atoms n,ca,c,o
zone $ARGV[2]-$ARGV[3]:$ARGV[2]-$ARGV[3]
fit
ratoms n,ca,c
rzone $ARGV[2]-$ARGV[3]\n";
close FIT1;

open (FIT2, ">fit2.$id");
print FIT2 "quiet
reference $ARGV[0]
mobile $ARGV[1]
atoms n,ca,c,o
zone -$st:-$st
zone $en-:$en-
fit
ratoms n,ca,c
rzone $ARGV[2]-$ARGV[3]\n";
close FIT2;

# calculate RMSD with ProFit

$out1 = `profit < fit.$id | grep \"RMS\"| tail -1`;
chomp($out1);
@res1 = split(/:/,$out1);

$out2 = `profit < fit2.$id | grep \"RMS\" | tail -1`;
chomp($out2);
@res2 = split(/:/,$out2);

print "RMSD: \t\t LOCAL: $res1[1] \t\t GLOBAL: $res2[1] \n"; 

# delete temporary files:
#`rm fit.$id fit2.$id`;

