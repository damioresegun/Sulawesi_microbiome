#!/bin/bash
# Script to carry out quality assessment of de novo assembled genomes using QUAST
# Set up with multiple functions for different types of assembled data.
# Input data can be raw assemblies (rawQuas); decontaminated assemblies (cleanQuas); medaka corrected assemblies (meddyQuas);
# pilon corrected assemblies (pilQuas); repeat masked assemblies (maskQuas or maskQuasScript) or complete annotated genomes (genoQuas)
set -e
#exec 2>&1 | tee VsPk_Flye_Quastlog.txt
Assem=$1 # path to folder holding directories of isolate assemblies i.e path/to/all/isolate/assemblies that contains directories of isolate assemblies. Can be raw, medaka, decontaminated 
# input for masked assemblies and annotated genomes differ. These simply need the path to the folder containing isolate fasta files to be assessed
output=$2 # path to output directory
Refre=$3 # path to the reference genome to be compared against
RefGf=$4 # path to the reference genome's gff file
THREADS=$5 # number of threads

cd $Assem
# function to assess the raw de novo genomes.
rawQuas() {
for folder in ${Assem}/*
do
    echo "The folder you are working on is $folder"
	echo $folder
    fna=$(basename "$folder")
	echo $fna
	fname="${fna%*.*}"
	#fname="${fnamer%*_*.*}"
    echo $fname
    #cd $folder
	mkdir -p $output/${fname}/Raw
    pwd
    echo
    echo "Running Quast on Lapp reference..."
    echo     
	#flye raw
	run="quast -t $THREADS ${Assem}/${fname}/assembly.fasta -o ${output}/${fname}/Raw/ -r ${Refre} -g ${RefGf} --large -f --circos"	
	echo $run
	eval $run
    cd $Assem
	echo "QUAST for ${fna} done successfully"
done
}
# function to assess the cleaned and decontaminated genomes.
cleanQuas() {
for folder in ${Assem}/*
do
    echo "The folder you are working on is $folder"
	echo $folder
    fna=$(basename "$folder")
	echo $fna
	fname="${fna%*.*}"
	#fname="${fnamer%*_*.*}"
    echo $fname
    #cd $folder
	mkdir -p $output/${fname}/Clean
    pwd
    #cd ${Assem}/${fname}
    echo
    echo "Running Quast on Lapp reference..."
    echo 
	run="quast -t $THREADS ${Assem}/${fname}/clean_assembly.fasta -o ${output}/${fname}/Clean/ -r ${Refre} -g ${RefGf} --large -f --circos"
	echo $run
	eval $run
	cd $Assem
done
}
# function to assess the completed and annotated genomes.
genoQuas() {
for f in ${Assem}/*
do
    echo "The isolate you are working on is $f"
	echo $f
    fna=$(basename "$f")
	echo $fna
	fname="${fna%*.*}"
	#fname="${fnamer%*_*.*}"
    echo $fname
    #cd $folder
	#mkdir -p $output/${fname}/Raw
	#mkdir -p $output/${fname}/Clean
    pwd
    #cd ${Assem}/${fname}
    #echo "Running quast Vs Ernest Reference..."
    #/shelf/apps/dro/conda/envs/quast/lib/python3.6/site-packages/quast-5.0.2-py3.6.egg-info/scripts/quast.py -t 12 /storage/home/users/dro/All_De_Novo/newAssem/sks339/sks339.contigs.fasta -o /storage/home/users/dro/Assembly_QC/newAssemblyRun/${fname}/quast/VsErnest/ -r /storage/home/users/dro/Index/PkErnestRefSeq.fasta -g /storage/home/users/dro/Index/PkErnestRefSeq.gff3 --large -f
    #echo "Quast vs Ernest Reference done..."
    echo
    echo "Running Quast on annotated genomes..."
    echo 
    mkdir -p $output/${fname}/Annotated
    #/shelf/apps/dro/conda/envs/quast/lib/python3.6/site-packages/quast-5.0.2-py3.6.egg-info/scripts/quast.py -t 12 /storage/home/users/dro/All_De_Novo/newAssem/sks339/sks339.contigs.fasta -o /storage/home/users/dro/Assembly_QC/newAssemblyRun/${fname}/quast/VsLapp/ -r /storage/home/users/dro/Index/PkRefreSeq.fna -g /storage/home/users/dro/Index/PkErnestRefSeq.gff3 --large -f
    
	##real quast command
	##generic
	quast -t $THREADS ${Assem}/${fname}.fasta -o ${output}/${fname}/Annotated -r ${Refre} -g ${RefGf} --large -f --circos
	###clean
	#quast -t $THREADS ${Assem}/${fname}/clean_assembly.fasta -o ${output}/${fname}/Clean/ -r ${Refre} -g ${RefGf} --large -f --circos
	
	
    cd $Assem
	
	### quast command for companion outputs
	#quast -t 16 $Assem/${fname}*.fasta -o $output/${fname} -r $Refre -g $RefGf --large -f 
	done
}
# function to assess the medaka corrected genomes.
meddyQuas() {
for folder in ${Assem}/*
do
    echo "The folder you are working on is $folder"
	echo $folder
    fna=$(basename "$folder")
	echo $fna
	fname="${fna%*.*}"
	#fname="${fnamer%*_*.*}"
    echo $fname
    #cd $folder
	mkdir -p $output/${fname}/Medaka
    pwd
    #cd ${Assem}/${fname}
    echo
    echo "Running Quast on Lapp reference..."
    echo 
	run="quast -t $THREADS ${Assem}/${fname}/${fname}_consensus.fasta -o ${output}/${fname}/Medaka/ -r ${Refre} -g ${RefGf} --large -f --circos"
	echo $run
	eval $run
	cd $Assem
done
}
# function to assess the pilon correct genomes
pilQuas() {
for folder in ${Assem}/*
do
    echo "The folder you are working on is $folder"
	echo $folder
    fna=$(basename "$folder")
	echo $fna
	fname="${fna%*_*}"
	#fname="${fnamer%*_*.*}"
    echo $fname
    #cd $folder
	mkdir -p $output/${fname}/Pilon
    pwd
    #cd ${Assem}/${fname}
    echo
    echo "Running Quast on Lapp reference..."
    echo 
	run="quast -t $THREADS ${Assem}/${fname}_iter2/${fname}.fasta -o ${output}/${fname}/Pilon/ -r ${Refre} -g ${RefGf} --large -f --circos"
	echo $run
	eval $run
	cd $Assem
done
}
# function to assess the masked genomes. This function is tuned to be called from a specified quality assesment script
maskQuasScript() {
for folder in ${Assem}/*
do
    echo "The folder you are working on is $folder"
	echo $folder
    fna=$(basename "$folder")
	echo $fna
	fname="${fna%*_*}"
	#fname="${fnamer%*_*.*}"
    echo $fname
    #cd $folder
	outy=$output/$fname
	mkdir -p $outy
    pwd
    #cd ${Assem}/${fname}
    echo
    echo "Running Quast on Lapp reference..."
    echo 
	run="quast -t $THREADS $folder -o $outy -r ${Refre} -g ${RefGf} --large -f --circos"
	echo $run
	eval $run
	cd $Assem
done
}
# function to assess the masked genomes. This function is tuned for this script and to be called from the command line
maskQuas() {
for folder in ${Assem}/*
do
    echo "The folder you are working on is $folder"
	echo $folder
    fna=$(basename "$folder")
	echo $fna
	fname="${fna%*_*}"
	#fname="${fnamer%*_*.*}"
    echo $fname
    #cd $folder
	mkdir -p $output/${fname}/Masked
	outy=$output/${fname}/Masked/
    pwd
    #cd ${Assem}/${fname}
    echo
    echo "Running Quast on Lapp reference..."
    echo 
	run="quast -t $THREADS $folder -o $outy -r ${Refre} -g ${RefGf} --large -f --circos"
	echo $run
	eval $run
	cd $Assem
done
}
racoQuasScript() {
for folder in ${Assem}/*
do
    echo "The folder you are working on is $folder"
	echo $folder
    fna=$(basename "$folder")
	echo $fna
	#fname="${fnamer%*_*.*}"
    #echo $fname
    #cd $folder
	outy=$output/$fname
	mkdir -p $outy
    pwd
    #cd ${Assem}/${fname}
    echo
    echo "Running Quast on Lapp reference..."
    echo 
	run="quast -t $THREADS $f/${fna}_iter4.fasta -o ${output}/${fna}/RaconPolished -r ${Refre} -g ${RefGf} --large -f --circos"
	echo $run
	eval $run
	cd $Assem
done
}