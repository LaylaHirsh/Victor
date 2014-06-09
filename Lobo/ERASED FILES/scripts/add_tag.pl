#! /usr/bin/perl -w

if (!defined(@ARGV)) {
    die "you must specify an input file\n"; 
}

$tag = $ARGV[1]; #uc(substr($ARGV[0], 0, 4));

open(FILE,"<$ARGV[0]") or die "file not found\n";

while (defined($item = <FILE>)) {

    if (substr($item, 0, 4) eq "ATOM") {
	$tmp = substr($item, 0, 66) . "      " . $tag;
	print "$tmp\n"; 
    }
    else  {
	print "$item\n"; 
    }

}




