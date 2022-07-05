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
import shutil
import subprocess

from Scripts.Tools import filter_fastq_file, makeDirectory
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
def align(isolate, file, temp_save_dir, fastq_dir_out, THREADS, stats, REFERENCE):
    reference_mmi = REFERENCE + ".mmi"
    if os.path.isfile(reference_mmi):
        print('Reference index exists. Index file will be used')
    else:
        print('Reference index does not exist')
        print('A index will be generated and used')
        runIndy = ' '.join(["minimap2 -x map-ont -d", reference_mmi, REFERENCE])
        print(runIndy)
        subprocess.call(runIndy, shell=True)
    print('Alignment starting...')
    aln = os.path.join(temp_save_dir, isolate + "_VsRef.bam")
    runAly = ' '.join(["minimap2 -ax map-ont", reference_mmi, file, "-t", str(THREADS),
           "| samtools view -@", str(THREADS), "-b - | samtools sort -@",
           str(THREADS), "-o", aln, "-"])
    print(runAly)
    subprocess.call(runAly, shell=True)
    print('Alignment done')
    stt = os.path.join(stats, isolate + "_FlagstatMappedVsRef_stats.txt")
    if os.path.exists(stt):
        os.remove(stt)
        os.mknod(stt)
        pass
    else:
        os.mknod(stt)
    runSamSt = ' '.join(["samtools flagstat --threads", str(THREADS), aln, ">", stt])
    alnEx = fastq_dir_out + "/" + isolate + "VsRef_unmapped.bam"
    alnOut = fastq_dir_out + "/" + isolate + "VsRef_unmapped.fastq"
    runSamEx = ' '.join(["samtools view --threads", str(THREADS), "-f 4 -b", aln, ">", alnEx])
    runBamFq = ' '.join(["bedtools bamtofastq -i", alnEx, "-fq", alnOut])
    print(runSamSt)
    subprocess.call(runSamSt, shell=True)
    print(runSamEx)
    subprocess.call(runSamEx, shell=True)
    print(runBamFq)
    subprocess.call(runBamFq, shell=True)
    return alnOut


def DNA_align(ready_path, dnaIsolate, statsDir, temp_align_out, aligned_out,
                threads, reference):
    for isolate in dnaIsolate:
        # set the right filtered/demultiplexed file
        file = os.path.join(ready_path, isolate + "fastq.gz")
        # set/make the stats directory
        statss = os.path.join(statsDir, "Alignment_Vs_Reference", isolate)
        makeDirectory(statss)
        # run the alignment of DNA against the reference
        hs_reads = align(isolate, file, temp_align_out, aligned_out, threads,
                            statss, reference)
        # Run de-duplication
        print('Running de-deuplicaton on ' + isolate)
        alnRename = aligned_out + "/" + isolate + "VsRef_unmapped_renamed.fastq"
        ''' filter_fastq_file function carries out deduplication '''
        filter_fastq_file(hs_reads, alnRename)
        print('De-duplication complete for ' + isolate)
        print("The reads have been renamed and saved as: " + alnRename)


def cDNA_align(ready_path, cdnaIsolate, creads, outDir, scripts, makeCref, 
                adapters, reference, creference, gff, threads, memory, aligned_out):
    for isolate in cdnaIsolate:
        # set the right filtered/demultiplexed file
        file = os.path.join(ready_path, isolate + "fastq.gz")
        # set/make the stats directory
        cdOut = os.path.join(outDir, "cDNA_Processing")
        makeDirectory(cdOut)
        # set the path to the cdna processing script
        cdProc = os.path.join(scripts, "cDNA_Processing.py")
        # check if the user chose to make a transcriptome assembly
        if makeCref is True:
            # make the command
            runCdProcy = ' '.join([cdProc, "-s reads -o", cdOut, "-r", creads[0], creads[1], 
                            "-a", adapters, "-hr", reference, "-hg", gff, 
                            "-p", str(threads), "-m", memory, "-ur", file])
            print(runCdProcy)
            # run the command
            subprocess.call(runCdProcy, shell=True)
            print('cDNA transcriptome generated')
        else:
            ''' if the user does not want to make a transcriptome, then align their reads
                against the given transcriptome using the cdna_processing script '''
            runCdProcy = ' '.join([cdProc, "-s genome -o", cdOut, "-p", str(threads), "-ur", file, 
                        "-t", creference])
            print(runCdProcy)
            # run the command
            subprocess.call(runCdProcy, shell=True)
        # get just the isolate name without the 'DNA or cDNA'
        iaw = isolate.split("_")[0]
        fileOut = os.path.join(cdOut, "Alignment", iaw + "_cDNA_VsTranscriptome.fastq")
        hfFileOut = os.path.join(aligned_out, isolate + "VsRef_unmapped.fastq")
        # copy and rename in the new location
        shutil.copy2(fileOut, hfFileOut)
        print('cDNA reads aligned against transcriptome and host reads separated')
        print('Host-free cDNA reads are now in ' + aligned_out)
        ''' filter_fastq_file function carries out deduplication'''
        print('Running de-deuplication on ' + isolate)
        # set the path for the deduplicated fastq
        hfRnm = os.path.join(aligned_out, isolate + "VsRef_unmapped_renamed.fastq")
        filter_fastq_file(hfFileOut, hfRnm)
        print('De-duplication complete for ' + isolate)
        print("The reads have been renamed and saved as: " + hfRnm)
