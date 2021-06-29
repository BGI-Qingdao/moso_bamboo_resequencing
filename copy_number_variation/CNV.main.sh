genome=$1
bam_list=$2
outdir=$3
basename=$4

bin=1000  # Reserved interface for bin length to detect CNV. The default length is currently 1 kb.

Bin=`cd \`dirname $0\`;pwd`

help(){
    cat <<HELP
     CNV.main.sh -- detect population-scaled CNV regions using Control-FREEC
 
 [usage]
  CNV.main.sh <genome path> <bam.list> <outdir> <basename>
 
 [note]
  <bam.list> format (with two columns): sample_id bam_path 
  
HELP
}

[ -z $1 ] && { help; exit 1;}

# step1. calculate GC profile using ggcount
## genome length
mkdir -p $outdir/genome/
perl $Bin/freec.split_fa.pl $genome $bin $outdir/genome $outdir/genome/ref.len

## write config file for ggcount
echo "[general]
chrLenFile = $outdir/genome/ref.len
ploidy = 2
GCcontentProfile = $outdir/genome/GC_profile.cnp
window = $bin
step = $bin
chrFiles = $outdir/genome
samtools = $Bin/samtools
bedtools = $Bin/bedtools

[sample]
inputFormat = BAM
mateOrientation = FR" > $ourdir/gccount.conf

cd $outdir   # have to move into the workdir
$Bin/gccount -conf $outdir/gccount.conf
mv $outdir/GC_profile.cnp $outdir/genome/GC_profile.cnp

# Step2. run Control-FREEC for each sample
sed -i 's/\s\+/\t/g' $bam_list
sample_num=`cat $bam_list | wc -l`

for line in $(seq 1 $sample_num); do
  id=`head -n $line $bam_list | tail -n 1 | cut -f 1`
  bam=`head -n $line $bam_list | tail -n 1 | cut -f 2`
  mkdir -p $outdir/$id
  
  # write config file for each sample
  echo "[general]
chrLenFile = $outdir/genome/ref.len
ploidy = 2
GCcontentProfile = $outdir/genome/GC_profile.cnp
window = $bin
step = $bin
chrFiles = $outdir/genome
outputDir = $outdir/$id

samtools = $Bin/samtools
bedtools = $Bin/bedtools

[sample]
mateFile = $bam
inputFormat = BAM
mateOrientation = FR" > $ourdir/$id/freec.conf
  # run control-freec
  $Bin/freec -conf $ourdir/$id/freec.conf
done

