use strict;
use warnings;

my ($hwe, $window, $bed, $out) = @ARGV;

open FL, $bed;
my %len;
while (<FL>) {
    chomp;
    my @tmp = split;
    $len{$tmp[0]} = $tmp[2];
};
close FL;

my %stat;
open FL, "tail -n +2 $hwe |";
while (<FL>) {
    chomp;
    my @tmp = split;
    my @num = split /\//, $tmp[2];
    my $het = $num[1]/($num[0]+$num[1]+$num[2]);
    my $index = int(($tmp[1]-1)/$window);
    $stat{$tmp[0]}{$index}{'num'}++;
    #$stat{$tmp[0]}{$index}{'max'} = $tmp[1];
    $stat{$tmp[0]}{$index}{'het'} += $het;
};
close FL;

open FLS, ">$out";
foreach my $chr (sort {$a<=>$b} keys %stat) {
    my $chr_len = $len{$chr};
    my $max_index = int(($chr_len-1)/$window);
    foreach my $index (0 .. $max_index) {
    #foreach my $index (sort {$a<=>$b} keys %{$stat{$chr}}) {
        my $start = $index * $window + 1;
        my $end = ($index + 1 ) * $window;
        if (exists $stat{$chr}{$index}) {
            my $avg_ho = sprintf("%.6f", $stat{$chr}{$index}{'het'}/$stat{$chr}{$index}{'num'});
            print FLS "$chr\t$start\t$end\t$stat{$chr}{$index}{'num'}\t$avg_ho\n";
        } else {
            print FLS "$chr\t$start\t$end\t0\t0\n";
        };
    };
};
close FLS;

