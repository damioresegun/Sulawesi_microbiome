#!/usr/bin/env python3
'''
 Name: DNA_Processing
 Author: Damilola R Oresegun
 MacKenzie Institute for Early Diagnostics
 April 2022

 Rationale: Acts as the sub-pipeline for DNA sequences for the
           NanoMetaPipe package
'''
# import modules
import os
import subprocess
#################### State functions #####################################
##########################################################################
''' Align reads against reference genome
 Function is to use minimap2 to carry out alignment against the chosen
 reference genome provided. The function checks for an index file and 
 if it does not exist, creates one and uses that. Note: the index file
 must be made using minimap2. Hence non-minimap2 indexes will not be
 detected. Alignment is carried out using the map-ont preset of minimap2
 and bam files are automatically created and reads which do not align
 are extracted, saved and converted into fastq files to use downstream.
 minimap2, samtools and bedtools must be in PATH to work '''
##########################################################################
def align(isolate,file,temp_save_dir,fastq_dir_out,THREADS,stats,REFERENCE):
    reference_mmi = REFERENCE + ".mmi"
    if os.path.isfile(reference_mmi):
        print('Reference index exists. Index file will be used')
    else:
        print('Reference index does not exist')
        print('A index will be generated and used')
        indy = ("minimap2 -x map-ont -d", reference_mmi, REFERENCE)
        runIndy = ' '.join(indy)
        print(runIndy)
        subprocess.call(runIndy, shell=True)
    print('Alignment starting...')
    aln = temp_save_dir + "/" + isolate + "_VsRef.bam"
    aly = ("minimap2 -ax map-ont", reference_mmi, file, "-t", str(THREADS),
           "| samtools view -@", str(THREADS), "-b - | samtools sort -@",
           str(THREADS), "-o", aln, "-")
    runAly = ' '.join(aly)
    print(runAly)
    subprocess.call(runAly, shell=True)
    print('Alignment done')
    stt = stats + "/" + isolate + "_FlagstatMappedVsRef_stats.txt"
    if os.path.exists(stt):
        os.remove(stt)
        os.mknod(stt)
        pass
    else:
        os.mknod(stt)
    samSt = ("samtools flagstat --threads", str(THREADS), aln, ">", stt)
    alnEx = fastq_dir_out + "/" + isolate + "VsRef_unmapped.bam"
    alnOut = fastq_dir_out + "/" + isolate + "VsRef_unmapped.fastq"
    samEx = ("samtools view --threads", str(THREADS), "-f 4 -b", aln, ">", alnEx)
    bamFq = ("bedtools bamtofastq -i", alnEx, "-fq", alnOut)
    runSamSt = ' '.join(samSt)
    runSamEx = ' '.join(samEx)
    runBamFq = ' '.join(bamFq)
    print(runSamSt)
    subprocess.call(runSamSt, shell=True)
    print(runSamEx)
    subprocess.call(runSamEx, shell=True)
    print(runBamFq)
    subprocess.call(runBamFq, shell=True)
    return alnOut
