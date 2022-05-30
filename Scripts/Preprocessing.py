#!/usr/bin/env python3
''' Name: Preprocessing
 Author: Damilola R Oresegun
 MacKenzie Institute for Early Diagnostics
 April 2022
 Rationale: Acts as the first step of the NanoMetaPipe package.
 Carries out demultiplexing, and filtering of reads
 Requires basecalled nanopore long reads
'''
# import modules
import os
import subprocess
#################### State functions #####################################
##########################################################################
"""Demultiplexing
Purpose: Function takes in the path to the newly created demultiplexed
folder in the given output folder. The function checks if the user chose
to use either guppy or qcat as a demultiplexer. Qcat is the default as
Guppy would require the use of a GPU. It is advised that is working on a
shared cluster, the qcat option is chosen.
"""
def demultip(INP_DIR, dem_dir, DEMULP_CHOICE, THREADS, KIT):
    try:
        # check if the user selected qcat for demultiplexing
        if DEMULP_CHOICE == "qcat":
            # sets the input for qcat to the basecalled reads
            dem_INP = INP_DIR + "/*"
            # constructs the qcat command
            qdem = ("cat", dem_INP, "| qcat -b", dem_dir, 
            "--detect-middle -t", str(THREADS), "--trim -k", KIT)
        # joins the individual words into a single command string
            runQdem = ' '.join(qdem)
            # print for logging
            print(runQdem)
            # run qcat command
            #subprocess.call(runQdem, shell=True)
            print('Demultiplexing complete')
        # check if the user selected guppy for demultiplexing
        elif DEMULP_CHOICE == "guppy":
            # constructs the quppy command
            gdem = ("guppy_barcoder -i", INP_DIR, "-s", dem_dir, 
                    "--barcode_kits", KIT,
                    "-r -q 0 -t", str(THREADS), "--compress_fastq -x auto",
                    "--detect_mid_strand_barcodes --trim_barcodes", 
                    "--trim_adapters")
            runGdem = ' '.join(gdem)
            print(runGdem)
            #subprocess.call(runGdem, shell=True)
            print('Demultiplexing complete')
    except OSError as error:
        print(error)
##########################################################################

##########################################################################
"""Filtering
Filtering using Nanofilt to filter for length and quality based upon
user preference. The user will have chosen to carry out filtering 
and if decided, what length and quality to filter by. Outputs are
saved in a newly generated directory name 'Filtered_Demultiplexed_Reads'
in the output directory indicated by the user.
"""
def filt_qc(dem_dir,barcode,isolate,OUT_DIR,FILT_LENGTH,FILT_QUAL):
    """ temp = barcode + "_" + isolate
    stats_dir = os.path.join(stats, "Filtered_Demultiplexed_Reads", temp)
    if os.path.exists(stats_dir):
        pass
    else:
        os.makedirs(stats_dir) """
    print('Starting NanoFilt')
    # making a temporary variable to hold the barcode name
    temp = barcode + ".fastq.gz"
    # making the path to the demultiplexed fastq
    file_in = os.path.join(dem_dir, temp)
    # make the output folder
    filt_out = os.path.join(OUT_DIR, "Filtered_Demultiplexed_Reads")
    if os.path.exists(filt_out):
        pass
    else:
        os.mkdir(filt_out)
    # overwrite the temporary variable with the renamed fastq
    temp = isolate + ".fastq.gz"
    filt_file_out = os.path.join(filt_out, temp)
    # construct the nanofilt command
    filtSt = ("gunzip -c", file_in, "|NanoFilt -l", str(FILT_LENGTH), "-q",
              str(FILT_QUAL), "| gzip >", filt_file_out)
    runFiltSt = ' '.join(filtSt)
    print(runFiltSt)
    # run nanofilt command
    #subprocess.call(runFiltSt, shell=True)
    ''' returns the path where the filtered reads are saved and the stats
     directory where the results are saved '''
    return filt_out
##########################################################################

##########################################################################
"""Pre-Processing QC for raw and filtered
To carry out some basic QC of the reads prior to carrying out alignment.
The pre and post filtering reads can be put through this function
Function uses NanoQC, NanoStat and FastQC to carry out QC checks
Outputs will be placed in the Stats folder and named appropriately"""
def run_QC(file,barcode,stats,ofile,THREADS):
    # check if the provided folder exists already
    if os.path.exists(stats):
        pass
    else:
        os.makedirs(stats)
    # make an internal variable to not affect global
    file_in = file
    # construct the nanostat command
    nanSt = ("NanoStat", "--fastq", file_in, "--outdir", stats, "-n", ofile)
    runNanSt = ' '.join(nanSt)
    print(runNanSt)
    # run the nanostat command
    #subprocess.call(runNanSt, shell=True)
    print('nanoStat complete for ' + barcode)
    print('Proceeding to nanoQC')
    # construct the nanoQC command
    nanQ = ("nanoQC", "-o", stats, file_in)
    runNanQ = ' '.join(nanQ)
    print(runNanQ)
    # run the nanoQC command
    #subprocess.call(runNanQ, shell=True)
    print('nanoQC complete for ' + barcode)
    print('Proceeding to FastQC')
    # construct the fastqc command
    fatq = ("fastqc", "-t", str(THREADS), "-o", stats, file_in)
    runFatq = ' '.join(fatq)
    print(runFatq)
    # run the fastqc command
    #subprocess.call(runFatq, shell=True)
    print('FastQC complete')
##########################################################################