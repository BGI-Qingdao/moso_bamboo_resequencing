#!/bin/env perl
#================================================================
# Version: 2.0
# Author: Sun Shuai                E-mail: sunshuai@genomics.cn
# Last modified: 2020-02-21 17:13
# Filename: bed_region_ann.pl
#================================================================
use strict;
use warnings;
use Getopt::Long;
use Data::Dumper;

my $usage = <<END;

    $0 -- Annotate span gene region for bed file

 [usage]
  perl $0 [option] --gff <in.gff> --bed <in.bed> --out <out>

 [option]
   --col    three cols to define one region [1,2,3]
   --extend extend length around defined BED region [0]
   --help   help information 

 [version]
  2.0 -- ouput all genes for each region with details (including: mRNA, CDS, and non-CDS).
  1.0 -- only output one gene for each region, should be not used any more.
 
 [note]
  gff and bed is 1-based in this script
END
die "$usage\n" if (@ARGV == 0);
my ($gff,$bed,$out);
my ($col, $extend);
my ($help);
GetOptions(
    'gff=s'     => \$gff,
    'bed=s'     => \$bed,
    'out=s'     => \$out,
    'col:s'     => \$col,
    'extend:i'  => \$extend,
    'help!'     => \$help,
) or die;
$col = '1,2,3' unless (defined $col);
$extend = 0 unless (defined $extend);
die "$usage\n" if ($help);
my ($c1,$c2,$c3) = split /,/,$col;
$c1--; $c2--; $c3--;

my $bin = 1000000;  # bin length: 1 Mb


# deal with gff info
open FL,$gff;
#my %gff;

my %length; # feature 'length' 

my %mRNA;
my %CDS;

my %non_CDS_content;
$non_CDS_content{'intron'} = 1;  # for most genes have intron, so set as one default part of non-CDS.

while (<FL>){
    chomp;
    next if (/^\s*$/);
    next if (/^#/);
    my @tmp = split;
    if ($tmp[2] eq 'mRNA'){
        my $id = $1 if ($tmp[8] =~ /ID=([^\s;]+)/);
        my $bin_start = int($tmp[3]/$bin);
        my $bin_end   = int($tmp[4]/$bin);
        foreach my $p ($bin_start .. $bin_end) {
            push(@{$mRNA{$tmp[0]}{$p}}, [$id, $tmp[3], $tmp[4]]); # id, start, end, length
        };
        $length{$tmp[0]}{$id}{'mRNA'} = $tmp[4] - $tmp[3] + 1;
        #$gff{$tmp[0]}{'mRNA'}{$id}[0] = $tmp[3];
        #$gff{$tmp[0]}{'mRNA'}{$id}[1] = $tmp[4];
    } elsif ($tmp[2] eq 'CDS'){
        my $id = $1 if ($tmp[8] =~ /Parent=([^\s;]+)/);
        push(@{$CDS{$tmp[0]}{$id}}, [$tmp[3],$tmp[4]]);
        $length{$tmp[0]}{$id}{'CDS'} += ($tmp[4] - $tmp[3] + 1);
    } else {
        $non_CDS_content{$tmp[2]} = 1;
    }
};
close FL;

my $non_CDS_content = join(",", sort keys %non_CDS_content);

open FL,$bed;
open FLS,">$out";
my $time = `date "+%F %T"`;
print FLS "# START at $time";
print FLS "## USED cols: ", $c1+1, ",", $c2+1, ",", $c3+1,"\n";
print FLS "## EXTEND distance: $extend\n";
print FLS "## non_CDS include: $non_CDS_content [just skip tags which do not belong to gene parts]\n";
print FLS "## ADDED col: mRNA_ID:mRNA[length,covered_length,covered_rate],CDS[length,covered_length,covered_rate],non-CDS[length,covered_length,covered_rate]\n";

while (my $line = <FL>){
    chomp $line;
    next if ($line =~ /^\s*$/);
    if ($line =~ /^#/){
        print FLS "$line\n";   
        next;
    };
    my @tmp = split /\s+/,$line;
    
    # bed region
    my ($chr,$start,$end) = ($tmp[$c1],$tmp[$c2],$tmp[$c3]);
    $start = $start - $extend;
    $end   = $end + $extend;
    #foreach my $id (keys %{$gff{$chr}{'mRNA'}}){
    
    my $bin_start = int($start/$bin);
    my $bin_end   = int($end/$bin);
   
    my $ann_result = '';

    # all cover genes
    my %mRNA_cov;
    foreach my $p ($bin_start .. $bin_end) {
        next unless (exists $mRNA{$chr}{$p});
        foreach my $i (@{$mRNA{$chr}{$p}}) {
            my $m_id    = $i->[0];
            my $m_start = $i->[1];
            my $m_end   = $i->[2];
            if ($m_start <= $end and $m_end >= $start){
                my @four = sort {$a <=> $b} ($start, $end, $m_start, $m_end);
                $mRNA_cov{$m_start}{$m_end}{'id'}  = $m_id;
                $mRNA_cov{$m_start}{$m_end}{'cov'} = $four[2] - $four[1] + 1;
            };
        }
    };
    #print Dumper(\%mRNA_cov),"\n";
    # for other features of mRNA
    foreach my $m_start (sort {$a<=>$b} keys %mRNA_cov) {
        foreach my $m_end (sort {$a<=>$b} keys %{$mRNA_cov{$m_start}}) {
            my $m_id = $mRNA_cov{$m_start}{$m_end}{'id'};
            
            # covered
            my $CDS_cov = 0;
            foreach my $i (@{$CDS{$chr}{$m_id}}){
                my $c_start = $i->[0];
                my $c_end   = $i->[1];
                
                #print "$m_id\t$c_start\t$c_end\n";

                if ($c_start <= $end and $c_end >= $start) {
                    my @four = sort {$a <=> $b} ($start, $end, $c_start, $c_end);
                    $CDS_cov += ($four[2] - $four[1] + 1);
                    #print "$CDS_cov\n";
                };
            };
            
            my $mRNA_len = $length{$chr}{$m_id}{'mRNA'};
            my $mRNA_cov = $mRNA_cov{$m_start}{$m_end}{'cov'};
            my $mRNA_rate = sprintf("%.2f", $mRNA_cov/$mRNA_len);
            
            my $CDS_len  = $mRNA_len;
            $CDS_len = $length{$chr}{$m_id}{'CDS'} if (exists $length{$chr}{$m_id}{'CDS'});
            my $CDS_rate = 0;
            if ($CDS_len > 0) {
                $CDS_rate = sprintf("%.2f", $CDS_cov/$CDS_len);
            };

            my $non_CDS_len = $mRNA_len - $CDS_len;
            my $non_CDS_cov = $mRNA_cov - $CDS_cov;
            my $non_CDS_rate = 0;
            if ($non_CDS_len > 0) {
                $non_CDS_rate = sprintf("%.2f", $non_CDS_cov/$non_CDS_len);
            };

            $ann_result .= "$m_id:mRNA[$mRNA_len,$mRNA_cov,$mRNA_rate],CDS[$CDS_len,$CDS_cov,$CDS_rate],non-CDS[$non_CDS_len,$non_CDS_cov,$non_CDS_rate];";
            #print "$ann_result\n";
        };
    };
    print FLS "$line\t$ann_result\n";
};
close FLS;
close FL;

