snp_file=$1
outgroup=$2
chr_num=$3
prefix=$4

Bin=`cd \`dirname $0\`;pwd`

help(){
    cat <<HELP
     Balance_selection.pipeline.sh -- detect regions under balance selection using BalLeRMix

 [usage]
  Balance_selection.pipeline.sh <snp file prefix> <outgroup.id.list> <chromosome number> <out_prefix>

 [note]
  <snp file prefix>: must have .vcf.gz .bim .fam .bed
  <outgroup.id.list>: only one sample allowed, format: family_id individua_id. It must be in snp file
  <chromosome number>: the chromosome must be numbered and the number of chromosome should not be more than 99.
  
HELP
}

[ -z $1 ] && { help; exit 1;}

# step1. get ancestral allele from the outgroup SNP information, only the allele of homozygous SNPs was considered. 
plink --bfile $snp_file --keep $outgroup --freq --out ${prefix}.outgroup --chr-set $chr_num
awk '$5==0 {print $2,$4}' ${prefix}.outgroup.frq | sed 's/\s\+/\t/g' | sed 's/:/\t/' > ${prefix}.anc.allele.txt &

# step2. convert VCF format to the file that can be used for BalLeRMix
perl $Bin/vcf2input.pl ${prefix}.anc.allele.txt $.vcf.gz ${prefix}.all.input.txt &
  ## remove sites with no variants 
awk '$4>0' ${prefix}.all.input.txt > ${prefix}.tmp
mv ${prefix}.tmp ${prefix}.all.input.txt

# step3. run BalLeRMix to calculate B2 statistics for each chromosome one by one
for i in $(seq 1 $chr_num); do
    head -n 1 ${prefix}.all.input.txt > ${prefix}.$i.input.txt
    tail -n +2 ${prefix}.all.input.txt | awk -v id=$i '$1==id' | cut -f 2- >> ${prefix}.$i.input.txt
    python2 $Bin/BalLeRMix_v2.2.py -i ${prefix}.$i.input.txt --getSpect --spect ${prefix}.$i.input.spect
    python2 $dir/BalLeRMix_v2.2.py -i ${prefix}.$i.input.txt --spect ${prefix}.$i.input.spect --physPos -o ${prefix}.$i.B2_stat.txt
done
head -n 1 ${prefix}.1.B2_stat.txt > ${prefix}.whole_genome.B2_stat.txt
for i in $(seq 1 $chr_num); do tail -n +2 ${prefix}.$i.B2_stat.txt | sed "s/^/$i\t/";done >> ${prefix}.whole_genome.B2_stat.txt

# step4. statistics
perl $Bin/array_stat.pl --col 4 --in all.B2_stat.txt --out all.B2_stat.stats

# this part is for moso bamboo population resequecing: I use the cutoff of top 0.5%, you may choose the appropriate cutoff according to your data 
tail -n +2 all.B2_stat.txt | awk '$4>270.883424' | perl -ne '@tmp=split; print "$tmp[0]\t", $tmp[1]-1, "\t$tmp[1]\n"' |  bedtools merge -i /dev/stdin -n -d 10000 | awk '$4>1'> top0.05.peak.bed
$Bin/../find_gene/bed_region_ann.pl --gff ../../00.genome/moso.hic_merge.gff --bed top0.05.peak.bed --extend 5000 --out top0.005.peak.ext_5k.gene_ann.bed
grep -v '^#' top0.005.peak.ext_5k.gene_ann.bed | cut -f 5 | sed 's/;/\n/g' | grep mRNA | cut -d ':' -f 1 | sort | uniq > top0.005.peak.ext_5k.gene.list
