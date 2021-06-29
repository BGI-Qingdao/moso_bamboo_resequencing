use strict;
use warnings;

my ($in, $bin, $outdir, $outlen) = @ARGV;

my $len = $bin;

$/ = '>';
open FL,$in;
<FL>;

my $i = 1;
while (my $seq = <FL>){
    chomp $seq;
    my $id = $1 if ($seq =~ /^(\S+)/);
    $seq =~ s/.*?\n//;
    $seq =~ s/\n//g;
    my $x = length($seq);
    if ($x >= $len){
        $seq =~ s/(\w{60})/$1\n/;
        $seq =~ s/\n$//;
        
        open FLS,">$outdir/$id.fa";
        print FLS ">$id\n$seq";
        close FLS;

        open FLS, ">$outlen"
        print FLS "$i\t$id\t$x\n";
        close FLS;

        $i++;
    };
};
close FL;
