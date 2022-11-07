#!/bin/bash
# Bash script to call the NanoMetaPipe python script
#
# enter variable
NanoMetaPipe=/home/doresegu/scratch/private/JCS_MetaGenome_Project/NanoMetaPipe.py
reference=/home/doresegu/scratch/private/JCS_MetaGenome_Project/Index/Macaca_nemestrina_reference.fna.gz
gff=/home/doresegu/scratch/private/JCS_MetaGenome_Project/Index/Macaca_nemestrina_genomic.gff.gz
basecalled=/home/doresegu/scratch/private/JCS_MetaGenome_Project/MFMRCFS1622_BiM004/Basecalled/pass
barcodes="barcode07 barcode08"
isolates="MFMRCFS1622_DNA MFMRCFS1622_dscDNA"
output=/home/doresegu/scratch/private/JCS_MetaGenome_Project/OctoberOutputs/MFMRCFS1622_BiM004
threads=24
demultiplexer="qcat"
cdna_reads="/home/doresegu/scratch/private/JCS_MetaGenome_Project/Index/SRR10248516_1.fastq
            /home/doresegu/scratch/private/JCS_MetaGenome_Project/Index/SRR10248516_2.fastq"
krakenDB=/home/doresegu/scratch/private/Kraken_DB
cdna_adapters=/home/doresegu/scratch/private/JCS_MetaGenome_Project/Index/TruSeq3.fa
cref=/home/doresegu/scratch/private/JCS_MetaGenome_Project/CDNA/Assembly/trinity_out/Trinity-GG.fasta # path to already made transcriptome
maxMem="150G" # maximum memory
kraken="kraken2"
krakThres=2
sequence_type="both"
DnaFilt=500
CdnaFilt=100
bracken="bracken"
brackenThres=10
ncbi_db="/mnt/shared/apps/databases/ncbi"
kraken_mode="both"

mkdir -p $output
touch ${output}/logfile.txt
# both command and transcriptome path command
run="python3 $NanoMetaPipe -b $basecalled -c $barcodes -e $isolates -r $reference -o $output -kb $krakenDB -kr $kraken -w -d $demultiplexer -dfl $DnaFilt -cfl $CdnaFilt -t $threads -g $gff -s $sequence_type -cr $cref -mxm $maxMem -br $bracken -bt $brackenThres -nd $ncbi_db -km $kraken_mode -kt $krakThres 2>&1 | tee ${output}/logfile.txt"

echo $run
eval $run
