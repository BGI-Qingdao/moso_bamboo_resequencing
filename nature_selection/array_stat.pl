#!/bin/env perl
use strict;
use warnings;
use Getopt::Long;

my $usage = <<END;
    perl $0 [Option] --in <in.file> --out <out.stat>
 
 [Option]
  --min     min value [null]
  --max     max value [null]
  --col     which column to use [1]
  --help    help info
END
die "$usage\n" if (@ARGV == 0);

my ($in, $out, $col, $help);
my ($min_cutoff, $max_cutoff);

GetOptions(
    'in=s'      => \$in,
    'out=s'     => \$out,
    'min:s'     => \$min_cutoff,
    'max:s'     => \$max_cutoff,
    'col:i'     => \$col,
    'help!'     =>\ $help,
) or die;
die "$usage\n" if ($help);

$min_cutoff = 'null' if (!defined $min_cutoff);
$max_cutoff = 'null' if (!defined $max_cutoff);
$col = 1 if (!defined $col);
$col = $col - 1;

# store all data
my @array;
open FL, $in;
while (<FL>) {
    chomp;
    next if (/^\s*$/);
    next if (/^#/);
    my @tmp = split;
    my $num = $tmp[$col];
    next if ($num =~ /[^\d\-\+\.eE]/);
    
    # value filter
    if ($min_cutoff ne 'null') {
        next if ($tmp[$col] < $min_cutoff);
    };
    if ($max_cutoff ne 'null') {
        next if ($tmp[$col] > $max_cutoff);
    };

    push (@array, $num);
};
close FL;

# sort
my @array_sort = sort {$a<=>$b} @array;

#
my $item_num = $#array_sort + 1;   die "the item number is to small (<4)\n" if ($item_num < 4);
my $min   = $array_sort[0];
my $max   = $array_sort[-1];
my $avg   = &calc_avg(\@array_sort);
my $sd    = &calc_sd(\@array_sort);
my $Q1    = &calc_per(\@array_sort, 0.25);
my $Q2    = &calc_per(\@array_sort, 0.50);
my $Q3    = &calc_per(\@array_sort, 0.75);
my $low01  = &calc_per(\@array_sort, 0.001);
my $low05  = &calc_per(\@array_sort, 0.005);
my $low1  = &calc_per(\@array_sort, 0.01);
my $low2  = &calc_per(\@array_sort, 0.02);
my $low5  = &calc_per(\@array_sort, 0.05);
my $low10 = &calc_per(\@array_sort, 0.10);

my $top01  = &calc_per(\@array_sort, 0.999);
my $top05  = &calc_per(\@array_sort, 0.995);
my $top1  = &calc_per(\@array_sort, 0.99);
my $top2  = &calc_per(\@array_sort, 0.98);
my $top5  = &calc_per(\@array_sort, 0.95);
my $top10 = &calc_per(\@array_sort, 0.90);

# ouput
open FLS, ">$out";
my $stat = <<END;
> BASE part
value number: $item_num
> RANGE part
min: $min  max: $max
avg: $avg  sd:  $sd
> Quartile part
Q1: $Q1  Q2: $Q2  Q3: $Q3
> TOP LOW (<=) part
low0.1: $low01 low0.5: $low05 low1: $low1  low2: $low2  low5: $low5  low10: $low10
> TOP HIGH (>=) part
top0.1: $top01 top0.5: $top05 top1: $top1  top2: $top2  top5: $top5  top10: $top10
END
print FLS $stat;

#
sub calc_avg {
    my $array = $_[0];
    
    my $total = 0;
    my $num = 0;
    foreach my $item (@$array) {
        $total += $item;
        $num ++;
    };
    my $avg = sprintf("%.6f", $total/$num);
    $avg =~ s/0+$//;
    $avg =~ s/\.$//;

    return $avg;
};

sub calc_sd {
    my $array = $_[0];
    
    my $total = 0;
    my $num = 0;
    #avg
    foreach my $item (@$array) {
        $total += $item;
        $num ++;
    };
    my $avg = $total/$num;
    # sd
    my $q = 0;
    foreach my $item (@$array) {
        $q += ($item - $avg) ** 2;
    };

    my $sd = sprintf( "%.6f", ($q/($num-1)) ** 0.5 );
    $sd =~ s/0+$//;
    $sd =~ s/\.$//;
    
    return $sd;   
};

# 
sub calc_per {
    my ($array, $per) = @_;
    my @array = @$array;
    
    my $point = ($#array + 1 + 1) * $per;

    my $int = int($point);
    my $float = $point - $int;
    
    my $result;
    if ($int == 0) {
        $result = sprintf( "%.6f", $array[0]);
    } elsif ($int > $#array) {
        $result = sprintf( "%.6f", $array[-1]);  
    } else {
        $result = sprintf( "%.6f", $array[$int-1] + ($array[$int] - $array[$int-1])*$float );
    };
    $result =~ s/0+$//;
    $result =~ s/\.$//;

    return $result;
};

