#!/usr/bin/env python3
#
# Name of pipeline: NanoMetaPipe
#
# Pipeline to take basecalled metagenomic nanopore reads through
# preprocessing to final assembled bins for analysis.
#
# MacKenzie Institute for Early Diagnosis (2022)
# Author: Damilola Oresegun, Peter Thorpe
#
import sys
import os
import configparser
import subprocess
import argparse
import mappy as mp
import logging.handlers
import time
from Bio import SeqIO
import shutil
from pathlib import Path
from Scripts.AssemblyQC import run_AssemStats, raw_Quast
from Scripts.Racon_Medaka import *
# from pycits.tools import convert_fq_to_fa
# from pycits.metapy_tools import make_folder
##########################################################################
def get_args():
    parser = argparse.ArgumentParser(description="Pipeline for metagenomics " +
                                     "data pre-processing and analysis for " +
                                     "nanopore sequenced data." +
                                     "Ensure that the tools required are in " +
                                     "in the PATH!", add_help=False)
    file_directory = os.path.realpath(__file__).split("NanoMetaPipe")[0]
    if not os.path.isfile(os.path.join(file_directory, "NanoMetaPipe.py")):
        file_directory = os.path.join(file_directory, "NanoMetaPipe")
    if not os.path.isfile(os.path.join(file_directory, "NanoMetaPipe.py")):
        print("Can't locate the correct path to the NanoMetaPipe pipeline script")
    ##########################################################################
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument("-b", "--basecalled",
                               dest="Basecalled_reads",
                               action="store",
                               type=str,
                               help="Folder of basecalled reads " +
                               "in FASTQ format. Example: path/to/" +
                               "exp_folder/pass. " +
                               "Does not work with raw reads. " +
                               "Default is for tests only",
                               required=True)
    required_args.add_argument("-c", "--barcodes",
                               dest='Barcodes',
                               action="store",
                               nargs='+',
                               help="List of barcodes used. " +
                               "e.g. -c barcode01 barcode02",
                               required=True)
    required_args.add_argument("-e", "--isolates",
                               dest='Isolates',
                               action="store",
                               nargs='+',
                               help="List of isolates used. e.g " +
                               "-e isolate1 isolate2 " +
                               "Must correspond to the order "
                               "of barcodes given " +
                               "i.e. barcode01=isolate1.",
                               required=True)
    required_args.add_argument("-r", "--reference",
                               dest="Reference_Genome",
                               action="store",
                               type=str,
                               help="Path to the reference genome to align " +
                               "against. Default is Macaca nemestrina",
                               required=True)
    ##########################################################################
    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument("-d", "--demultiplexer",
                               dest='Demultiplexer',
                               action="store",
                               type=str,
                               default="qcat",
                               choices=["qcat", "guppy"],
                               help="Choice of " +
                               "demultiplexer to use. Guppy or Qcat. " +
                               "Note: If guppy, you " +
                               "will need to have access to a " +
                               "GPU on your system to use. " +
                               "Default is qcat")
    optional_args.add_argument("-fl", "--filter_length",
                               dest='Filter_length',
                               action="store",
                               type=int,
                               default="50",
                               help="Filter length threshold. " +
                               "Default is 50")
    optional_args.add_argument("-flw", "--flowcell",
                               dest='Flowcell',
                               action="store",
                               type=str,
                               default="FLO-MIN106",
                               help="The flowcell used for this " +
                               "experiment. Default is FLO-MIN106")
    optional_args.add_argument("-fq", "--filter_quality",
                               dest='Filter_Quality',
                               action="store",
                               type=int,
                               default="10",
                               help="Phred quality threshold to filter " +
                               "demultiplexed reads. Default is 10")
    optional_args.add_argument("-ft", "--filter",
                               dest='Filter_option',
                               action="store_false",
                               help="Option to carry " +
                               "out filtering. Default is True. " +
                               "Turn off with just -ft")
    optional_args.add_argument("-g", "--gff",
                               dest='GFF',
                               action="store",
                               type=str,
                               help="The GFF file of your reference " +
                               "genome.")
    optional_args.add_argument("-k", "--kit",
                               dest='Sequencing_kit',
                               action="store",
                               type=str,
                               default='LSK109',
                               help="The sequencing kit used, without SQK" +
                               ". The default is LSK109")
    optional_args.add_argument("-n", "--name",
                               dest='Experiment_Name',
                               action="store",
                               type=str,
                               help="Give the name of " +
                               "the experiment.")
    optional_args.add_argument("-o", "--out_dir",
                               dest="Output_folder", action="store",
                               default=os.path.join(file_directory,
                                                    "Pipeline_Output"),
                               type=str,
                               help="Path to the output directory. " +
                               "Default will create the output file " +
                               "in your current working directory.")
    optional_args.add_argument("-p", "--expansion",
                               dest='Expansion_kit',
                               action="store",
                               type=str,
                               default="NBD104",
                               help="The expansion kit used for " +
                               "barcoding the isolates sequenced " +
                               "without 'EXP'. Default is NBD104")
    optional_args.add_argument("-rd", "--redo_demulp",
                               dest='Re_demultiplex',
                               action="store_true",
                               help="Option to carry out " +
                               "demultiplexing alone. To be used " +
                               "to re-generate demultiplexed reads " +
                               "if previously deleted. Off by default. " +
                               "Turn on with -rd")
    optional_args.add_argument("-t", "--threads",
                               dest='Threads',
                               action="store",
                               type=int,
                               default="24",
                               help="Number of threads. Default is 24")
    optional_args.add_argument("-w", "--cleanup",
                               dest='Cleanup',
                               action="store_true",
                               help="Option to clean up files as you go " +
                               "Note: this means that some generated " +
                               "will be deleted e.g. Basecalled outputs, " +
                               "demultiplexed and reads before alignment " +
                               "to the reference will be deleted; among " +
                               "other downstream temp files and folders. " +
                               "Read the log files to see what data are " +
                               "deleted. Cleanup off by default. Turn on " +
                               "with just -w")
    optional_args.add_argument("-h", "--help",
                               action="help",
                               default=argparse.SUPPRESS,
                               help="Displays this help message")
    ##########################################################################
    args = parser.parse_args()
    return args, file_directory


##########################################################################
# set global variables
args, FILE_DIRECTORY = get_args()
#WK_DIR = os.getcwd()
INP_DIR = args.Basecalled_reads
OUT_DIR = args.Output_folder
EXP_NAME = args.Experiment_Name
REFERENCE = args.Reference_Genome
BARCODES = args.Barcodes
ISOLATES = args.Isolates
DEMULP_CHOICE = args.Demultiplexer
FLOWCELL = args.Flowcell
FILTER_PASS = args.Filter_option
THREADS = args.Threads
CLEAN = args.Cleanup
REDEMULP = args.Re_demultiplex
REF_GFF = args.GFF
#
# set conditional variables
#if EXP_NAME is None:
#    pass
#    PREFIX = os.path.split(OUT_DIR)[-1].split()[0]
#else:
#    OUT_DIR = os.path.join(OUT_DIR, EXP_NAME)
#    PREFIX = EXP_NAME
#
if DEMULP_CHOICE == "qcat":
    G_KIT = args.Sequencing_kit
    Q_KIT = "NBD103/" + args.Expansion_kit
else:
    G_KIT = "SQK-" + args.Sequencing_kit
    Q_KIT = "EXP-" + args.Expansion_kit
#
if FILTER_PASS is True:
    FILT_LENGTH = args.Filter_length
    FILT_QUAL = args.Filter_Quality
else:
    pass
if CLEAN is True:
    print("You have chosen to do clean up. Large files and directories" +
          "will be deleted as the pipelines progress. The most " +
          "consequential of these include the deletion of the " +
          "basecalled reads as these can simply be remade with " +
          "the basecalling script and demultiplexing can be " +
          "with the -rd parameter")
else:
    print("No clean up will be done as the pipeline progresses. " +
          "Please be mindful of the space this will take up.")
##########################################################################
# Making global directories
##########################################################################
dem_dir = os.path.join(OUT_DIR, "Demultiplexed")
if os.path.exists(dem_dir):
    pass
else:
    os.makedirs(dem_dir)
stats_dir = os.path.join(OUT_DIR, "Stats")
if os.path.exists(stats_dir):
    pass
else:
    os.mkdir(stats_dir)
temp_align_out = os.path.join(OUT_DIR, "Isolate_Reads_Aligned_Vs_Reference")
if os.path.exists(temp_align_out):
    pass
else:
    os.makedirs(temp_align_out)
aligned_out = os.path.join(OUT_DIR, "Host_Free_Reads")
if os.path.exists(aligned_out):
    pass
else:
    os.makedirs(aligned_out)
if FILTER_PASS is False:
    rnm_path = os.path.join(OUT_DIR, "Isolate_Demultiplexed_Reads")
    if os.path.exists(rnm_path):
        pass
    else:
        os.makedirs(rnm_path)
##########################################################################
# Demultiplexing
# Purpose: Function takes in the path to the newly created demultiplexed
# folder in the given output folder. The function checks if the user
# chose to use either guppy or qcat as a demultiplexer. Qcat is the
# default as Guppy would require the use of a GPU. It is advised that 
# if working on a shared cluster, the qcat option is chosen.
##########################################################################
def demultip(dem_dir, DEMULP_CHOICE):
    try:
        # check if the user selected qcat for demultiplexing
        if DEMULP_CHOICE == "qcat":
            # sets the input for qcat to the basecalled reads
            dem_INP = INP_DIR + "/*"
            # constructs the qcat command
            qdem = ("cat", dem_INP, "| qcat -b", dem_dir, "--detect-middle",
                    "-t", str(THREADS), "--trim -k", Q_KIT)
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
                    "--barcode_kits", Q_KIT,
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
# Filtering
# Filtering using Nanofilt to filter for length and quality based upon
# user preference. The user will have chosen to carry out filtering 
# and if decided, what length and quality to filter by. Outputs are
# saved in a newly generated directory name 'Filtered_Demultiplexed_Reads'
# in the output directory indicated by the user.
##########################################################################
def filt_qc(dem_dir,barcode,isolate,stats):
    #temp = barcode + "_" + isolate
    #stats_dir = os.path.join(stats, "Filtered_Demultiplexed_Reads", temp)
    #if os.path.exists(stats_dir):
    #    pass
    #else:
    #    os.makedirs(stats_dir)
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
    # returns the path where the filtered reads are saved and the stats
    # directory where the results are saved
    return filt_out,stats_dir
##########################################################################
# Pre-Processing QC for raw and filtered
# To carry out some basic QC of the reads prior to carrying out alignment.
# The pre and post filtering reads can be put through this function
# Function uses NanoQC, NanoStat and FastQC to carry out QC checks
# Outputs will be placed in the Stats folder and named appropriately
##########################################################################
def run_QC(file,barcode,stats,ofile):
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
# Align reads against reference genome
# Function is to use minimap2 to carry out alignment against the chosen
# reference genome provided. The function checks for an index file and 
# if it does not exist, creates one and uses that. Note: the index file
# must be made using minimap2. Hence non-minimap2 indexes will not be
# detected. Alignment is carried out using the map-ont preset of minimap2
# and bam files are automatically created and reads which do not align
# are extracted, saved and converted into fastq files to use downstream.
# minimap2, samtools and bedtools must be in PATH to work
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
    #subprocess.call(runAly, shell=True)
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
    #subprocess.call(runSamSt, shell=True)
    print(runSamEx)
    #subprocess.call(runSamEx, shell=True)
    print(runBamFq)
    #subprocess.call(runBamFq, shell=True)
    return alnOut
##########################################################################
# De-duplication of reads
# Reads will be checked for de-duplication. In many cases, there will be
# duplicated reads or duplicated read names. So this step checks for
# duplicated reads and then renames them. Duplicated reads are not
# deleted as there is no way to know if it is just the read names that 
# are duplicated or the entire reads
##########################################################################
def de_dup(isolate,in_file,file_dir):
    dedupp = os.path.join(FILE_DIRECTORY, "Scripts/remove_duplicated_fastq.py") ################################## come back to see if you need to hard code it like this. No. CHANGE the script to a function and just call the function here!!!!!!
    alnRename = file_dir + "/" + isolate + "VsRef_unmapped_renamed.fastq"
    dedup = (dedupp, "-i", in_file, "-o", alnRename)
    runDeDup = ' '.join(dedup)
    print('Running de-deuplication on ' + isolate)
    print(runDeDup)
    #subprocess.call(runDeDup, shell=True)
    print('De-duplication complete for ' + isolate)
    return alnRename
##########################################################################
# Assemble reads using metaFlye
# Reads which have been extracted from alignments against the reference
# genome will go through de novo genome assembly using the '--meta' flag
# of flye which allows for higher error reads and also allows for
# variance in coverage. Duplicated reads will need to be checked and
# renamed. This is done using a custom python script that checks the name
# of a read and renames it if it is duplicated. Then metaFlye takes place
##########################################################################
def run_flye(in_file,out_dir,threads):
    try:
        fly = ("flye --nano-raw", in_file, "--out-dir", out_dir, "--threads", str(threads), "--meta")
        runFly = ' '.join(fly)
        print(runFly)
        #subprocess.call(runFly, shell=True)
    except FileNotFoundError:
        print('There was an issue with the file. File not found')
    except SyntaxError:
        print('There was a syntax issue. Double check your syntax')
##########################################################################
# Check the stats of the assemblies. Starts with using assembly-stats
# for descriptive statistics. 
##########################################################################
# moved to external script. Left this incase I want to call something in here later
####################################################################################################################################################
# Run the script and functions
####################################################################################################################################################
# Setting up logging
if __name__ == '__main__':
    logger = logging.getLogger('NanoMetaPipe.py: %s' % time.asctime())
    # logger.setLevel(logging.DEBUG)
    logger.setLevel(logging.INFO)
    d_err_handler = logging.StreamHandler(sys.stdout)
    d_err_handler.setLevel(logging.INFO)
    d_err_formatter = logging.Formatter(
        '%(name)s - %(levelname)s - %(message)s')
    d_err_handler.setFormatter(d_err_formatter)
    logger.addHandler(d_err_handler)
    e_err_handler = logging.StreamHandler(sys.stderr)
    e_err_handler.setLevel(logging.ERROR)
    e_err_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    e_err_handler.setFormatter(e_err_formatter)
    logger.addHandler(e_err_handler)
    try:
        logstream = open(OUT_DIR + '/debug.log', 'w')
        d_err_handler_file = logging.StreamHandler(logstream)
        d_err_handler_file.setFormatter(d_err_formatter)
        d_err_handler_file.setLevel(logging.INFO)
        logger.addHandler(d_err_handler_file)
        logstream2 = open(OUT_DIR + '/error.log', 'w')
        e_err_handler_file = logging.StreamHandler(logstream2)
        e_err_handler_file.setFormatter(e_err_formatter)
        e_err_handler_file.setLevel(logging.ERROR)
        logger.addHandler(e_err_handler_file)
    except BaseException:
        outstr = "Debug log file could not be open for logging"
        logger.error(outstr)
        sys.exit(1)
    # Print the arguments to file
    logger.info("Command line: %s", ' '.join(sys.argv))
    logger.info("Starting: %s", time.asctime())
    ##########################################################################
    # Call functions
    ##########################################################################
    # Demultiplexing
    if REDEMULP is True:
        print("You are just re-generating the demultiplexed reads")
        demultip(dem_dir, DEMULP_CHOICE)
        print("Demultiplexed reads regenerated")
        sys.exit(1)
    else:
        demultip(dem_dir, DEMULP_CHOICE)
    # Zip the demultiplexed reads
    print('Zipping the fastq files')
    for i in os.listdir(dem_dir):
        zip_dem = os.path.join(dem_dir, i)
        pgSt = ("pigz --best", zip_dem)
        runPgSt = ' '.join(pgSt)
        print(runPgSt)
        subprocess.call(runPgSt,shell=True)
    print('FastQ files have been zipped')
    # qc the raw demultiplexed reads
    for i in BARCODES:
        ofile = i + ".txt"
        dem_file = dem_dir + "/" + i + ".fastq.gz"
        stats = os.path.join(stats_dir, "Raw_Demultiplexed_Reads", i)
        run_QC(dem_file,i,stats,ofile)
    print(THREADS)
    # filtering: if filter check is true
    if FILTER_PASS is True:
        count = 0
        for i in BARCODES:
            isola = ISOLATES[count]
            ready_path, filtered_stats = filt_qc(dem_dir, i, isola, stats_dir)
            print(
                'The raw demultiplexed reads have been successfully filtered ' + 
                'and saved in ' + ready_path)
            print('The QC stats for the filtered reads are saved in ' +
                  os.path.dirname(filtered_stats))
            print('Please remember that the files are now renamed')
            print(i + ' is now ' + isola)
            print('')
            print('Now running QC of filtered reads')
            ofile = ISOLATES[count] + ".txt"
            dem_file = ready_path + "/" + isola + ".fastq.gz"
            temp = i + "_" + isola
            stats = os.path.join(stats_dir, "Filtered_Demultiplexed_Reads", temp)
            run_QC(dem_file,i,stats,ofile)
            count += 1
    else:
        print('You chose to not filter your reads.')
        print('Your demultiplexed reads will be used to proceed')
        ready_path = rnm_path
        count = 0
        for i in BARCODES:
            isola = ISOLATES[count]
            dem_file = dem_dir + "/" + i + ".fastq.gz"
            rnm_file = ready_path + "/" + isola + ".fastq.gz"
            mv_stuff = ("mv", dem_file, rnm_file)
            rnMv_Stff = ' '.join(mv_stuff)
            print(rnMv_Stff)
            count += 1
    # align against reference genome and de-duplicate
    for i in ISOLATES:
        file = ready_path + "/" + i + ".fastq.gz"
        statss = os.path.join(stats_dir, "Alignment_Vs_Reference", i)
        if os.path.exists(statss):
            pass
        else:
            os.makedirs(statss)
        hs_reads = align(i,file,temp_align_out,aligned_out,THREADS,statss,REFERENCE)
        # de-duplication
        print('Running de-deuplication on ' + i)
        hs_rn_reads = de_dup(i,hs_reads,aligned_out)
        print('De-duplication complete for ' + i)
        print("The reads have been renamed and saved as: " + hs_rn_reads)
        
    # Carry out de novo genome sequencing
    for i in ISOLATES:
        toAlgn = aligned_out + "/" + i + "VsRef_unmapped_renamed.fastq"
        assemOut = os.path.join(OUT_DIR, "DeNoVo_Assembly", i)
        assemSt = os.path.join(stats_dir, "DeNoVo_Assembly", i)
        if os.path.exists(assemOut):
            pass
        else:
            os.makedirs(assemOut)
        if os.path.exists(assemSt):
            pass
        else:
            os.makedirs(assemSt)
        # assemble
        run_flye(toAlgn,assemOut,THREADS)
        print('Assembly completed with metaFlye')
        #### running stats on the raw assembly
        # assembly-stats
        assemFile = os.path.join(assemOut, "assembly.fasta")
        assemO = assemSt + "/RAW_assembly_stats.txt"
        #run_AssemStats(i,assemFile,assemO)
        # quast
        qAssmO = os.path.join(assemSt, "Quast")
        #rawQ = raw_Quast(assemFile,qAssmO,REFERENCE,REF_GFF,str(THREADS))
        #print('Quast assessment done and saved in ' + rawQ)
        # copy the assembly file to a more accessible point
        print('The generated assembly will be copied to a new folder')
        assemS = os.path.join(OUT_DIR, "Assemblies", "Raw")
        if os.path.exists(assemS):
            pass
        else:
            os.makedirs(assemS)
        coAssemFile = os.path.join(assemS, (i + ".fasta"))
        print(coAssemFile)
        #shutil.copy2(assemFile, coAssemFile)
        print('The generated assembly files have been copied to ' + assemS)
        # do some clean up at this point
        # check if cleanup option is given
    if CLEAN is True:
        print('Cleaning up the DeNoVo Assembly')
        print('They are not deleted, just zipped for archiving')
        denClean = os.path.join(OUT_DIR, "DeNoVo_Assembly")
        denCleanr = denClean + ".tar.gz"
        tarr = 'tar --use-compress-program="pigz --best -p ' + str(THREADS) + '"'
        pgSt = (tarr, "-cf", denCleanr, denClean)
        runPgSt = ' '.join(pgSt)
        print(runPgSt)
        subprocess.call(runPgSt,shell=True)
        # check if the zipped folder is deleted. If not then delete it ONLY if the zipped folder exists
        if os.path.exists(denClean):
            if os.path.exists(denCleanr):
                shutil.rmtree(denClean)
            else:
                print('DeNovo Assembly files have not been zipped')
        else:
            if os.path.exists(denCleanr):
                print('DeNovo Assembly has already been zipped')
                pass
            else: 
                print('DeNovo Assembly folder and zipped folder could not be found')
                print('Did the de novo step run properly? You may have issues downstream from this')
        print('The de novo assembly has been zipped')
    else:
        print('No cleanup was chosen so the assembly folders have not been zipped')
        print('Please be aware that they may potentially take up a substantial amount of storage')
    ##### call the racon function for polishing
    # look for the assembly files in the new location
    racAssem = os.path.join(OUT_DIR, "Assemblies", "Racon")
    if os.path.exists(racAssem):
        pass
    else:
        os.makedirs(racAssem)
    for f in Path(assemS).glob('*'):
        f_name = os.path.basename(f).split(".")[0]
        print(f)
        print(f_name)
        # get the reads used to make the assembly
        readsPath = aligned_out + "/" + f_name + "VsRef_unmapped_renamed.fastq"
        # set the output folder for each isolate
        racOUT = os.path.join(OUT_DIR, "Racon_Polishing", f_name)
        print(readsPath)
        print(racOUT)
        # call the racon for 4 iterations
        racStats = runRacon(f_name,str(f),readsPath,racOUT,THREADS,stats_dir)
        print('The stats for the racon alignments can be found at ' + os.path.dirname(racStats))
        print('Copying the final racon iteraction output to the Assemblies folder')
        racFile = os.path.join(racAssem, (f_name + "_iteration4.fasta"))
        racc = racOUT + "/" + f_name + "_iter4.fasta"
        shutil.copy2(racc, racFile)
        print('Racon outputs copied')
        # now call medaka for 1 iteration
        medaOut = os.path.join(OUT_DIR, "Medaka")
        runMedaka(f_name,readsPath,racFile,medaOut,THREADS)




print('Pipeline done')
# print(args)
# print(FILE_DIRECTORY)
# print(args.Isolates[1])
# print(OUT_DIR)
# print(Q_KIT)
# print(G_KIT)
# print(FILTER_PASS)
# print(dem_dir)