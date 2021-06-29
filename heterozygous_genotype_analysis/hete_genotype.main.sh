genome=$1
vcf=$2
window_len=$3
outdir=$4
basename=$5

Bin=`cd \`dirname $0\`;pwd`

help(){
    cat <<HELP
     hete_genotype.main.sh -- detect high and low heterozygous genotype regions.
 
 [usage]
  hete_genotype.main.sh <genome.fa> <in.population.vcf.gz> <window length> <outdir> <basename>
 
 [note]
  (1) window length in bp.
  Needed software: VCFTools, BEDTtools (I have uploaded pre-complied software in this directory. If it don't not work, you should complie these tools yourself.)
 
 [example]
  sh hete_genotype.main.sh ../01.variation/01.SNP/Moso.SNP.final.ex_OG.vcf.gz 200000 ./window_200k/ Moso.SNP.final.ex_OG 
HELP
}

[ -z $1 ] && { help; exit 1;}

# Step1. count the number of heterozygous and homozygous genotypes for each SNP site using VCFTools --hwe
$Bin/vcftools --gzvcf $vcf --hwe --out $outdir/$basename

# Step2. caluate the ratio of heterozygous genotypes for each SNP site
perl $Bin/hwe.site_hete.pl $outdir/$basename.hwe $outdir/$basename.site_hete.txt

# step3. calulate the averge ratio of heterozygous genotypes in windows (I used 200kb there to intend to detect longer heterozygous regions)
perl $Bin/fa_bed.pl $genome $outdir/$basename.ref.len
perl $Bin/hwe.window_hete.pl $outdir/$basename.hwe 200000 $outdir/$basename.ref.len $outdir/$basename.window_hete.txt

# step4. calculate average and S.D. value
## cutoff of high and low
perl $Bin/avg_sd.pl $outdir/$basename.window_hete.txt $outdir/$basename.window_hete.stats
avg=`cut -f 1 $outdir/$basename.window_hete.stats`
sd=`cut -f 2 $outdir/$basename.window_hete.stats`
high_cutoff=`expr $avg + $sd`
low_cutoff=`expr $avg - $sd`
## long regions large than 5 Mb
awk -v cut=$high_cutoff '$5>cut' $outdir/$basename.window_hete.txt | $Bin/bedtools merge -d 900000 -i /dev/stdin | awk '$3-$2+1>=5000000' > $outdir/$base.high.SD_avg.bed 
awk -v cut=$low_cutoff  '$5<cut' $outdir/$basename.window_hete.txt | $Bin/bedtools merge -d 900000 -i /dev/stdin | awk '$3-$2+1>=5000000' > $outdir/$base.low.SD_avg.bed
