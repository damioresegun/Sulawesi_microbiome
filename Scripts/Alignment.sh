#!/bin/bash
# Script for aligning and sorting SAM files to BAM files for all sequences within the input folder and save to the output folder
# The script aligns an input isolate fastq against the human reference genome, then extracts unmapped reads
# Unmapped reads are then saved and converted to a fastq file and checked for metric statistics
# intermediate files will be deleted i.e. SAM files and unsorted BAM files. Sorted BAM files will be renamed to just .bam rather than _sorted.bam
##NOTE: MIGHT HAVE TO RUN THIS IN A SAMTOOLS CONDA ENVIRONMENT
set -e
#set input and output folders
input=$1	# input isolate fastq to be aligned
outSave=$2	# output directory to hold initial alignment files
fastqSave=$3	# output directory to save generated unmapped fastq file
THREADS=$4	# number of threads
statsFol=$5	# output directory to save statistical assessments
reference=$6	# path to the reference to check for contimination

mkdir -p $fastqSave
mkdir -p $statsFol
#make an unmapped and mapped folder
mkdir -p $outSave/Unmapped
unmpd=$outSave/Unmapped
mpd=$outSave/Mapped
mkdir -p $mpd

#get the filename
file=$(basename $input)
echo $file
#remove the extension
fname="${file%.*.*}"
echo $fname
echo
#set the alignment name to the path and the file
alname=${outSave}/${fname}
echo $alname
check=${reference}.mmi
echo $check
if [[ -f "$check" ]]; then
	# index the reference
	echo "Index file exists. Index file will be used"
else
	echo "Index does not exist"
	echo "indexing the reference genome"
	indy="minimap2 -x map-ont -d ${reference}.mmi $reference"
	echo $indy
	eval $indy
fi
#align to reference, convert to bam, and sort
aly="minimap2 -ax map-ont ${reference}.mmi ${input} -t $THREADS | samtools view -@ $THREADS -b - | samtools sort -@ $THREADS -o ${alname}VsRef.bam -"
echo $aly
eval $aly
echo "alignment done"
echo
echo "bam conversion done"
echo ${fname} >> ${statsFol}/${fname}_FlagstatMappedVsRef_stats.txt
samtools flagstat --threads $THREADS ${alname}VsRef.bam >> ${statsFol}/${fname}_FlagstatMappedVsRef_stats.txt
echo "Converted to fastq"
# extract unmapped reads
samtools view --threads $THREADS -f 4 -b ${alname}VsRef.bam > ${unmpd}/${fname}VsRef_unmapped.bam
echo "unmapped extraction done"
#convert unmapped reads to fastq
bedtools bamtofastq -i ${unmpd}/${fname}VsRef_unmapped.bam -fq ${unmpd}/${fname}VsRef_unmapped.fastq
echo
#check the stats
assembly-stats ${unmpd}/${fname}VsRef_unmapped.fastq > ${statsFol}/${fname}_VsRef_FastQ_UnmappedStats.txt
#move fastqs to their own folder
mv ${unmpd}/${fname}VsRef_unmapped.fastq $fastqSave
echo "all done"