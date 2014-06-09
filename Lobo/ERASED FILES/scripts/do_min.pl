#! /usr/bin/perl

if (@ARGV < 3)
{
    print "Usage:\n";
    print "do_min.pl <filename-model> <1stAA-loop> <lastAA-loop> [<id>]\n";

    die "Not enough arguments.\n";
}

if (@ARGV < 4)
{
    $id = "$$"; # this is the current process id!!!
}
else
{
    $id = $ARGV[3];
}

$ext = substr($ARGV[0], -3, 3);

# fix segment id for CHARMM input
`cat /dev/null > tmp1.pdb`;
`./add_tag.pl $ARGV[0] TMP1 > tmp1.pdb`;

# now build input files for CHARMM minimization
`charmm pdbfile=tmp1 toppar=/usr/local/src/c27b3/toppar/ < build.inp > build.$id.out`;

$st = $ARGV[1] - 1;
$en = $ARGV[2] + 1;

# minimization proper
`charmm pdbfile=tmp1 indstart=$st indend=$en < min-loop.inp > min.$id.out`;

$out1 = `grep '^EEEEE' min.$id.out`;
chomp($out1);
@resC = split(/  /,$out1);
print "CHARMM: \t $resC[1] \t $resC[2] \t $resC[3] \n";

print `cat /dev/null > etmp.$ext.pdb`;
print `sed -e "s/OT2/OXT/g
        s/OT1/O  /g
        s/HSD/HIS/g
        s/CD  ILE/CD1 ILE/g" < tmp1_min.pdb > etmp.$ext.pdb`;

#./LoopModelEval -i etmp2.pdb -m etmp2.pdb -s $st -e $en
$out2 =  `../../Energy/pdb2energy -i etmp.$ext.pdb`;
chomp($out2);
@res1 = split(/\n/,$out2);
@res11 = split(/\t/,$res1[0]);
@res21 = split(/\t/,$res1[1]);

print "RAPDF: \t\t $res21[0] \n SOLV: \t\t $res11[0]\n";

# clean up temporary files, except "etmp" model:

#`rm build.$id.out`; 

`cp etmp.$ext.pdb etmp.pdb`;


