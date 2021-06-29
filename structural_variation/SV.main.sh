genome=$1
bam_list=$2
outdir=$3
basename=$4


Bin=`cd \`dirname $0\`;pwd`

help(){
    cat <<HELP
     SV.main.sh -- detect population-scaled SV regions using Breakdancer, Manta, and SUVOVOR
 
 [usage]
  SV.main.sh <genome path> <bam.list> <outdir> <basename>
 
 [note]
  <bam.list> format (with two columns): sample_id bam_path 
  
HELP
}

[ -z $1 ] && { help; exit 1;}

# Step1. run Breakdancer and Manta for each sample, seprately
sed -i 's/\s\+/\t/g' $bam_list
sample_num=`cat $bam_list | wc -l`

for line in $(seq 1 $sample_num); do
  id=`head -n $line $bam_list | tail -n 1 | cut -f 1`
  bam=`head -n $line $bam_list | tail -n 1 | cut -f 2`
  mkdir -p $outdir/Manta_BreakDancer/$id
  # run Manta 
  $Bin/configManta.py --bam $bam --referenceFasta $genome --runDir=$outdir/Manta_BreakDancer/$id
  python2 $outdir/Manta_BreakDancer/$id/runWorkflow.py
  
  ## convert INV, and only retain INS, DEL and INV 
  python2 $Bin/convertInversion.py $Bin/samtools $genome $outdir/Manta_BreakDancer/$id/results/variants/candidateSV.vcf.gz > $outdir/Manta_BreakDancer/$id/candidateSV.vconvertINV.vcf
  
  # split by SV type: INS, DEL, and INV
  grep -E 'SVTYPE=INS' $outdir/Manta_BreakDancer/$id/candidateSV.vconvertINV.vcf > $outdir/Manta_BreakDancer/$id/INS.vcf
  grep -E 'SVTYPE=DEL' $outdir/Manta_BreakDancer/$id/candidateSV.vconvertINV.vcf > $outdir/Manta_BreakDancer/$id/DEL.vcf
  grep -E 'SVTYPE=INV' $outdir/Manta_BreakDancer/$id/candidateSV.vconvertINV.vcf > $outdir/Manta_BreakDancer/$id/INV.vcf
  
  # run BreakDaner
  ## Generating configuration file
  perl $Bin/bam2cfg.pl -v 100 -t $Bin/samtools -q 20 -o $outdir/Manta_BreakDancer/$id/breakdancer.cfg -h $bam
  $Bin/breakdancer-max -m 10000000 -q 20 -y 20 -r 1 $outdir/Manta_BreakDancer/$id/breakdancer.cfg 1> $outdir/Manta_BreakDancer/$id/breakdancer.ctx 2> $outdir/Manta_BreakDancer/$id/breakdancer.log
 
  ## only retain ITX, at least three reads to support this ITX and CTX
  grep -E "^#|ITX|CTX" $outdir/Manta_BreakDancer/$id/breakdancer.ctx | perl -ne '$num=$1 if (/bam\|(\d+)/); print if ($num>=3)' > $outdir/Manta_BreakDancer/$id/breakdancer.filt.ctx
  perl $Bin/ctx_2_vcf.pl $outdir/Manta_BreakDancer/$id/breakdancer.filt.ctx $outdir/Manta_BreakDancer/$id/breakdancer.filt.vcf
  grep -E 'SVTYPE=ITX' $outdir/Manta_BreakDancer/$id/breakdancer.filt.vcf > $outdir/Manta_BreakDancer/$id/ITX.vcf
  #grep -E 'SVTYPE=CTX' $outdir/Manta_BreakDancer/$id/breakdancer.filt.vcf > $outdir/Manta_BreakDancer/$id/CTX.vcf # not consider CTX
done

# Step2. merge SV using SURVIVOR, and filter MAF (>0.01)
for type in INS DEL INV ITX; do
  vcf_list=''
  for line in $(seq 1 $sample_num); do
    id=`head -n $line $bam_list | tail -n 1 | cut -f 1`
    vcf_list="$vcf_list $outdir/Manta_BreakDancer/$id/$type.vcf"
  done
  
  # merge
  $Bin/SURVIVOR merge $vcf_list 1000 1 1 1 0 50 $outdir/$type.merged.vcf
  
  # filter MAF 0.01
  $Bin/vcftools --vcf $outdir/$type.merged.vcf --maf 0.01 --max-maf 0.99 --recode --recode-INFO-all --out $outdir/$type.final.vcf
 done
