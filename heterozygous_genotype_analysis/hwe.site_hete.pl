use strict;
use warnings;

my ($hwe, $out) = @ARGV;

open FL, "tail -n +2 $hwe |";
open FLS, ">$out";
while (<FL>) {
    chomp;
    my @tmp = split;
    my @num = split /\//, $tmp[2];
    my $het = $num[1]/($num[0]+$num[1]+$num[2]);
    print FLS "$tmp[0]\t$tmp[1]\t$het\n";
};
close FL;
close FLS;
