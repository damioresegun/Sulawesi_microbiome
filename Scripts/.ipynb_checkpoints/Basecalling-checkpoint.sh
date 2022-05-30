#!/bin/bash
#$ -cwd
# Script to take in raw FAST5 nanopore files and take them through basecalling. 
# Currently set up for basecalling using guppy and then qcat demultiplexer using the demultiplexing script.
# Set up based on Guppy version 4.0.15
#set -e
# Set Paths
conEnv=$1
Raw_IN=$2 # input folder for raw FAST5 files
echo "Raw files are in $Raw_IN"
Save_OUT=$3   # output folder for the basecalling
FLOWCELL=$4	# state the flowcell used e.g. FLO-MIN106
KIT=$5  # sequencing kit used. e.g. SQK-RBK004
# path to the ONT Guppy
ONT=$6 # path to the nanopore basecalling package. e.g. ont-guppy/bin
# Begin
Base_OUT=$Save_OUT/Basecalled
Log_Ouut=$Save_OUT/LogFiles
mkdir -p $Base_OUT
mkdir -p $Log_Ouut
#mkdir -p $Dem_OUT
#srsh --x11 --partition=gpu --gpus=1
#source activate $conEnv
#guppy basecaller
gupbase=`$ONT/guppy_basecaller -i ${Raw_IN} --save_path ${Base_OUT} --flowcell $FLOWCELL --kit $KIT --disable_pings -r -v -q 0 -x auto --trim_adapters --gpu_runners_per_device 18 --chunks_per_runner 512 2>&1 | tee $Log_Ouut/Basecalling.txt`
#chk=`qcat --help 2>&1 | tee $Log_Ouut/Basecalling.txt`
echo $gupbase
#echo $chk
conda deactivate
exit
echo "All done"