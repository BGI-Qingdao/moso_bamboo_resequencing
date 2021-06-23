use strict;
use warnings;

my ($anc, $vcf, $out) = @ARGV;

my %anc;
open FL, $anc;
while (<FL>) {
    chomp;
    my @tmp = split;
    $anc{$tmp[0]}{$tmp[1]} = $tmp[2];
};
close FL;

if ($vcf =~ /\.gz$/) {
    open FL, "gzip -dc $vcf |";
} else {
    open FL, $vcf;
};

open FLS, ">$out";
print FLS "physPos\tgenPos\tx\tn\n";
while (<FL>) {
    chomp;
    my @tmp = split;
    next if (/^#/);
    
    my $anc_alle;
    if (exists $anc{$tmp[0]}{$tmp[1]}) {
        $anc_alle = $anc{$tmp[0]}{$tmp[1]};
    } else {
        next;
    };
    
    my @alle = split /,/, $tmp[4];
    unshift(@alle, $tmp[3]);
    
    my $anc_p;
    foreach my $p (0 .. $#alle) {
        $anc_p = $p if ($alle[$p] eq $anc_alle);
    };

    my ($x, $n) = (0,0);
    foreach my $p (10 .. $#tmp) {
        my ($a1, $a2) = ($1, $2) if ($tmp[$p] =~ /^(\S)\/(\S)/);
        next if ($a1 eq '.' or $a2 eq '.');
        $n += 2;
        if ($a1 ne $anc_p) {
            $x++;
        };
        if ($a2 ne $anc_p) {
            $x++;
        };
    };
    print FLS "$tmp[0]\t$tmp[1]\t",$tmp[1]/1000000, "\t$x\t$n\n";
};
close FLS;
close FL;
