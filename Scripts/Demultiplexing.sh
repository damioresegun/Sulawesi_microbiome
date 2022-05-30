#!/bin/bash
#$ -cwd

# Script to take in raw nanopore files and take them through demultiplexing.
# Set up for guppy demultiplexer and qcat demultiplexer. However qcat is recommended and left as default.
# Uncomment for the ones where necessary. Ensure that this is run in environments where the necessary tools are. 
# Set up for Qcat version 1.1.0

set -e
#Set Paths
conEnv=$1
Base_IN=$2   # input folder containing basecalled reads
Save_OUT=$3    # output folder for the demultiplexing
THREADS=$4   # number of threads
#KIT=$4  # sequencing kit used. this version is for guppy e.g. SQK-RBK004
KIT=$5    # sequencing kit used. this version is for qcat e.g. RBK004
#ONT=$6	# path to the ONT Guppy e.g. ont-guppy/bin

Dem_OUT=$Save_OUT/Demultiplexed
Log_Ouut=$Save_OUT/LogFiles
mkdir -p $Dem_OUT
mkdir -p $Log_Ouut
source activate $conEnv
#............................................
#                 QCAT
#............................................
#echo "The files are in $Base_IN"
qdem="cat $Base_IN/pass/* | qcat -b ${Dem_OUT} --detect-middle -t $THREADS --trim -k '$KIT' --guppy 2>&1 | tee $Log_Ouut/Demultiplexing.txt"
echo $qdem
eval $qdem
pigz --best $Dem_OUT/*.fastq
conda deactivate
echo "Demultiplexing done"