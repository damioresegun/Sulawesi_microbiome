#!/usr/bin/env python3
'''
 Name of pipeline: NanoMetaPipe

 Pipeline to take basecalled metagenomic nanopore reads through
 preprocessing to final assembled bins for analysis.

 MacKenzie Institute for Early Diagnosis (2022)
 Author: Damilola Oresegun, Peter Thorpe
'''
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
from Scripts.Preprocessing import *
from Scripts.DNA_processing import *
from Scripts.Deduplication import *
#from Scripts.convert_fa_to_fa import convert_file
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
                               "Does not work with raw reads.",
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
                               "-e isolate1_dna isolate2_cdna " +
                               "Must correspond to the order " +
                               "of barcodes given " +
                               "i.e. barcode01=isolate1_dna. " +
                               "Must give sequence type separated "+
                               "by an underscore",
                               required=True)
    required_args.add_argument("-kr", "--kraken",
                                dest="Kraken_PATH",
                                action="store",
                                type=str,
                                help="Full path to your kraken installation " +
                                "if it is not in your $PATH. If in your $PATH "+
                                "simply write kraken2",
                                required=True)
    required_args.add_argument("-kb", "--kraken_DB",
                                dest="Kraken_DBPATH",
                                action="store",
                                type=str,
                                help="Full path to your kraken database",
                                required=True)
    required_args.add_argument("-r", "--reference",
                               dest="Reference_Genome",
                               action="store",
                               type=str,
                               help="Path to the reference genome to align " +
                               "against. Default is Macaca nemestrina",
                               required=True)
    required_args.add_argument("-s", "--sequence_type",
                                dest="Sequence_Type",
                                action="store",
                                choices=["dna","cdna","both"],
                                type=str,
                                help="The type of sequence you are trying to " +
                                "analyse. Can be dna or cdna. Please indicate " +
                                "this in your isolate name or the pipeline " +
                                "will not initialise! Example of this is: " +
                                "isolate1_dna or isolate2_cdna. No defaults are " +
                                "set. Must signify either or both sequence types",
                                required=True)
    ##########################################################################
    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument("-cr", "--cdna_ref",
                                dest='cdna_reference',
                                action="store",
                                type=str,
                                help="If using both DNA and cDNA sequence types " +
                                "Please provide the path to the transcriptome " +
                                "assembly here. If you are only using cDNA, this " +
                                "parameter is not necessary!")
    optional_args.add_argument("-mcr", "--make_cdna_ref",
                                dest='make_cdna_reference',
                                action="store_true",
                                help="Option to generate a genome-guided " +
                                "transcriptome to remove host sequences from " +
                                "input cDNA sequences. Only works with cDNA option" +
                                ", does not work with -cr option. [Default: off], " +
                                "turn on with just [-mcr/--make_cdna_reference]. " +
                                "To be used in conjunction with -cd option"
                                )
    optional_args.add_argument("-cd", "--ref_cdna_reads",
                                dest='reference_cdna_reads',
                                action="store",
                                nargs='+',
                                help="Path to cDNA reads. If using cDNA reads " +
                                "and wish to make " +
                                "a transcriptome, use this option with the " +
                                "-mcr option."
                                )
    optional_args.add_argument("-ca", "--cdna_adapters",
                                dest="cDNA_Adapters",
                                action="store",
                                type=str,
                                help="Path to a FASTA file containing the " +
                                "adapters used for your cDNA sequencing reads. "+
                                "Must be used with the -mcr option")
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
                               "[Default: qcat]")
    optional_args.add_argument("-fl", "--filter_length",
                               dest='Filter_length',
                               action="store",
                               type=int,
                               default="50",
                               help="Filter length threshold. " +
                               "[Default: 50]")
    optional_args.add_argument("-flw", "--flowcell",
                               dest='Flowcell',
                               action="store",
                               type=str,
                               default="FLO-MIN106",
                               help="The flowcell used for this " +
                               "experiment. [Default: FLO-MIN106]")
    optional_args.add_argument("-fq", "--filter_quality",
                               dest='Filter_Quality',
                               action="store",
                               type=int,
                               default="10",
                               help="Phred quality threshold to filter " +
                               "demultiplexed reads. [Default: 10]")
    optional_args.add_argument("-ft", "--filter",
                               dest='Filter_option',
                               action="store_false",
                               help="Option to carry " +
                               "out filtering. [Default: On]. " +
                               "Turn off with just [-ft]")
    optional_args.add_argument("-g", "--gff",
                               dest='gff',
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
                               ". [Default: LSK109]")
    optional_args.add_argument("-kt", "--kraken_hit_threshold",
                                dest='Kraken_Hit_Threshold',
                                action="store",
                                type=int,
                                default=5,
                                help="A minimum number of groups that must" +
                                "be matched to place a contig into a " +
                                "taxonomic group")
    optional_args.add_argument("-n", "--name",
                               dest='Experiment_Name',
                               action="store",
                               type=str,
                               help="Give the name of " +
                               "the experiment.")
    optional_args.add_argument("-o", "--out_dir",
                               dest='Output_folder', action="store",
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
                               "without 'EXP'. [Default: NBD104]")
    optional_args.add_argument("-rd", "--redo_demulp",
                               dest='Re_demultiplex',
                               action="store_true",
                               help="Option to carry out " +
                               "demultiplexing alone. To be used " +
                               "to re-generate demultiplexed reads " +
                               "if previously deleted. [Default: Off]. " +
                               "Turn on with [-rd]")
    optional_args.add_argument("-t", "--threads",
                               dest='Threads',
                               action="store",
                               type=int,
                               default="24",
                               help="Number of threads. [Default: 24]")
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
                               "deleted. [Default: Off]. Turn on " +
                               "with just [-w]")
    optional_args.add_argument("-mxm", "--max-memory",
                                dest='Maximum_Memory',
                                action="store",
                                default="48G",
                                type=str,
                                help="Maximum memory to use for this " +
                                "pipeline. Bare in mind that you are " +
                                "likely going to be using fairly large " +
                                "read files and these will be read into " +
                                "memory. ZThe classification database " +
                                "will also be read into memory. So " + 
                                "recommended maximum memory should be " +
                                "at least 130G. You can try lower " +
                                "however, this could lead to failure " +
                                "in some steps such as filtering, assembly " +
                                "and classification. [Default: 48G]")
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
REF_GFF = args.gff
SEQ_TYP = args.Sequence_Type
CREF = args.cdna_reference
KRAK = args.Kraken_PATH
KRAKDB = args.Kraken_DBPATH
KRAK_THRESH = args.Kraken_Hit_Threshold
MXMEM = args.Maximum_Memory
MAKCREF = args.make_cdna_reference
CREADS = args.reference_cdna_reads
CADAP = args.cDNA_Adapters
SCPTS = os.path.join(FILE_DIRECTORY, "Scripts") # Scripts folder will be part of the package
#INDEX = os.path.join(FILEDIRETORY, "Index") # Index folder will be part of the package. WILL INCLUDE TruSeq.fa and readme file with links for the SRR sequences and the macaca nemestrina downloads
###########################################################################
# check if the isolate names are properly named
DNA_ISOLATE = []
CDNA_ISOLATE = []
for i in ISOLATES:
    iso = i.lower()
    if (iso.__contains__("_dna")):
        print('You have provided DNA sequences')
        print(i)
        tempd = DNA_ISOLATE
        tempd.append(i)
    elif (iso.__contains__("cdna")):
        print('You have provided cDNA sequences')
        print(i)
        tempc = CDNA_ISOLATE
        tempc.append(i)
    else:
        print('You have not provided the isolates in a satisfactory format')
        print('Do all your isolate names have _dna or _cdna or _dscdna?')
        print('Please look at the help and try again')
        sys.exit(1)
# check if a cDNA reference genome is given
if SEQ_TYP == "dna":
    pass
elif SEQ_TYP == "cdna":
    if MAKCREF is True:
        if not CREADS:
            print('You have chosen to make a transcriptome assembly but provided no reads')
            print('Please provide the reads and run again')
            sys.exit(1)
        else:
            print('You have chosen not to make a transcriptome')
            print(MAKCREF)
            pass
    else:
        CREF = REFERENCE
        pass
elif SEQ_TYP == "both":
    # check if the cdna reference is given
    if not REFERENCE:
        print('Did you provide the two needed references?')
        print('Please remember a DNA reference and a transcriptome reference '+
                'are required')
        print('Please try again')
        sys.exit(1)
    elif not CREF:
        if MAKCREF is True:
            if not CREADS:
                print('You have chosen to make a transcriptome assembly but provided no reads')
                print('Please provide the reads and run again')
                sys.exit(1)
            else:
                print('You are giving both DNA and cDNA reads and have selected to ' +
                'generate a transcriptome assembly.')
            if not CADAP:
                print('You have chosen to make a transcriptome assembly but provided ' +
                        'no adapters for them. Please use the -ca option to add the adapter file')
                sys.exit(1)
            else:
                pass
        else: 
            print('Did you provide the two needed references?')
            print('Please remember a DNA reference and a transcriptome reference '+
                    'are required')
            print('Please try again')
            sys.exit(1)
    else:
        print('All checks are correct. Continuing')
else:
    print('Invalid entry. Please choose between dna, cdna or both')
    sys.exit(1)
print(DNA_ISOLATE)
print(CDNA_ISOLATE)
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
    BRACK_LENGTH = FILT_LENGTH
else:
    pass
    BRACK_LENGTH = args.Filter_length
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
""" if MAKCREF is True:
    if not CREADS:
        print('You have chosen to make a transcriptome assembly but provided no reads')
        print('Please provide the reads and run again')
        sys.exit(1)
    else:
        pass
else:
    pass """
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
''' De-duplication of reads
 Reads will be checked for de-duplication. In many cases, there will be
 duplicated reads or duplicated read names. So this step checks for
 duplicated reads and then renames them. Duplicated reads are not
 deleted as there is no way to know if it is just the read names that 
 are duplicated or the entire reads '''
##########################################################################
def de_dup(isolate,in_file,file_dir):
    dedupp = os.path.join(FILE_DIRECTORY, "Scripts/remove_duplicated_fastq.py") ################################## come back to see if you need to hard code it like this. No. CHANGE the script to a function and just call the function here!!!!!!
    alnRename = file_dir + "/" + isolate + "VsRef_unmapped_renamed.fastq"
    dedup = (dedupp, "-i", in_file, "-o", alnRename)
    #filter_fastq_file(in_file,alnRename)
    runDeDup = ' '.join(dedup)
    print('Running de-deuplication on ' + isolate)
    print(runDeDup)
    #subprocess.call(runDeDup, shell=True)
    print('De-duplication complete for ' + isolate)
    return alnRename
##########################################################################
''' Assemble reads using metaFlye
 Reads which have been extracted from alignments against the reference
 genome will go through de novo genome assembly using the '--meta' flag
 of flye which allows for higher error reads and also allows for
 variance in coverage. Duplicated reads will need to be checked and
 renamed. This is done using a custom python script that checks the name
 of a read and renames it if it is duplicated. Then metaFlye takes place '''
##########################################################################
def run_flye(in_file,out_dir,threads):
    try:
        fly = ("flye --nano-raw", in_file, "--out-dir", out_dir, "--threads", str(threads), "--meta")
        runFly = ' '.join(fly)
        print(runFly)
        subprocess.call(runFly, shell=True)
    except FileNotFoundError:
        print('There was an issue with the file. File not found')
    except SyntaxError:
        print('There was a syntax issue. Double check your syntax')
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
    # Demultiplexing. Function: demultip is in the Preprocessing.py script
    if REDEMULP is True:
        print("You are just re-generating the demultiplexed reads")
        demultip(INP_DIR, dem_dir, DEMULP_CHOICE, THREADS, Q_KIT)
        print("Demultiplexed reads regenerated")
        sys.exit(1)
    else:
        demultip(INP_DIR, dem_dir, DEMULP_CHOICE, THREADS, Q_KIT)
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
        run_QC(dem_file,i,stats,ofile,THREADS)
    # filtering: if filter check is true
    if FILTER_PASS is True:
        count = 0
        for i in BARCODES:
            isola = ISOLATES[count]
            ready_path = filt_qc(dem_dir, i, isola,
                                                OUT_DIR,FILT_LENGTH,FILT_QUAL,)
            print(
                'The raw demultiplexed reads have been successfully filtered ' + 
                'and saved in ' + ready_path)
            print('Please remember that the files are now renamed')
            print(i + ' is now ' + isola)
            print('')
            print('Now running QC of filtered reads')
            ofile = ISOLATES[count] + ".txt"
            dem_file = ready_path + "/" + isola + ".fastq.gz"
            temp = i + "_" + isola
            stats = os.path.join(stats_dir, "Filtered_Demultiplexed_Reads", temp)
            run_QC(dem_file,i,stats,ofile,THREADS)
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
    ''' align against reference genome and de-duplicate
     first check if the sequence type is dna or cdna
     align function is in DNA_processing.py'''
    if SEQ_TYP == "dna":
        for i in DNA_ISOLATE:
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
    elif SEQ_TYP == "cdna":
        for i in CDNA_ISOLATE:
            file = ready_path + "/" + i + ".fastq.gz"
            print(file)
            cdOut = os.path.join(OUT_DIR, "cDNA_Processing")
            if os.path.exists(cdOut):
                pass
            else:
                os.makedirs(cdOut)
            cdProc = os.path.join(SCPTS, "cDNA_Processing.py")
            if MAKCREF is True:
                cdProcy = (cdProc, "-s reads -o", cdOut, "-r", CREADS[0], CREADS[1], 
                            "-a", CADAP, "-hr", REFERENCE, "-hg", REF_GFF, 
                            "-p", str(THREADS), "-m", MXMEM, "-ur", file, "2>&1 | tee bothlogger.txt")
                runCdProcy = ' '.join(cdProcy)
                print(runCdProcy)
                subprocess.call(runCdProcy, shell=True)
                print('cDNA transcriptome generated')
            else:
                cdProcy = (cdProc, "-s genome -o", cdOut, "-p", str(THREADS), "-ur", file, 
                            "-t", CREF)
                runCdProcy = ' '.join(cdProcy)
                print(runCdProcy)
                subprocess.call(runCdProcy, shell=True)
            iaw = i.split("_")[0]
            fileOut = cdOut + "/" + "Alignment" + "/" + iaw + "_cDNA_VsTranscriptome.fastq"
            hfFileOt = aligned_out + "/" + iaw + "_cDNAVsRef_unmapped.fastq" 
            shutil.copy2(fileOut, hfFileOt)
            print('cDNA reads aligned against transcriptome and host reads separated')
            print('Host-free cDNA reads are now in ' + aligned_out)
            #statss = os.path.join(stats_dir, "Alignment_Vs_Reference", i)
            #if os.path.exists(statss):
            #    pass
            #else:
            #    os.makedirs(statss)
            #hs_reads = align(i,file,temp_align_out,aligned_out,THREADS,statss,CREF)
            # de-duplication
            print('Running de-deuplication on ' + i)
            hs_rn_reads = de_dup(i,hfFileOt,aligned_out)
            print('De-duplication complete for ' + i)
            print("The reads have been renamed and saved as: " + hs_rn_reads)
    elif SEQ_TYP == "both":
        for i in DNA_ISOLATE:
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
        for i in CDNA_ISOLATE:
            file = ready_path + "/" + i + ".fastq.gz"
            cdOut = os.path.join(OUT_DIR, "cDNA_Processing")
            if os.path.exists(cdOut):
                pass
            else:
                os.makedirs(cdOut)
            cdProc = os.path.join(SCPTS, "cDNA_Processing.py")
            if MAKCREF is True:
                cdProcy = (cdProc, "-s reads -o", cdOut, "-r", CREADS[0], CREADS[1], 
                            "-a", CADAP, "-hr", REFERENCE, "-hg", REF_GFF, 
                            "-p", str(THREADS), "-m", MXMEM, "-ur", file, "2>&1 | tee bothlogger.txt")
                runCdProcy = ' '.join(cdProcy)
                print(runCdProcy)
                subprocess.call(runCdProcy, shell=True)
                print('cDNA transcriptome generated')
            else:
                cdProcy = (cdProc, "-s genome -o", cdOut, "-p", str(THREADS), "-ur", file, 
                            "-t", CREF)
                runCdProcy = ' '.join(cdProcy)
                print(runCdProcy)
                subprocess.call(runCdProcy, shell=True)
            iaw = i.split("_")[0]
            fileOut = cdOut + "/" + "Alignment" + "/" + iaw + "_cDNA_VsTranscriptome.fastq"
            hfFileOt = aligned_out + "/" + iaw + "_cDNAVsRef_unmapped.fastq" 
            shutil.copy2(fileOut, hfFileOt)
            print('cDNA reads aligned against transcriptome and host reads separated')
            print('Host-free cDNA reads are now in ' + aligned_out)
            #statss = os.path.join(stats_dir, "Alignment_Vs_Reference", i)
            #if os.path.exists(statss):
            #    pass
            #else:
            #    os.makedirs(statss)
            #hs_reads = align(i,file,temp_align_out,aligned_out,THREADS,statss,CREF)
            # de-duplication
            print('Running de-deuplication on ' + i)
            hs_rn_reads = de_dup(i,hfFileOt,aligned_out)
            print('De-duplication complete for ' + i)
            print("The reads have been renamed and saved as: " + hs_rn_reads)
    # Split the pathways
    ############ de novo genome sequencing for the cDNA ##############
    if SEQ_TYP == "cdna":
        # first carry out kraken classification
        for i in CDNA_ISOLATE:
            toAlgn = aligned_out + "/" + i + "VsRef_unmapped_renamed.fastq"
            krakOut = os.path.join(OUT_DIR, "Kraken", i)
            if os.path.exists(krakOut):
                pass
            else:
                os.makedirs(krakOut)
            # kraken
            cKrak = (KRAK, "--db", KRAKDB, toAlgn, "--threads", str(THREADS),
                    "--output", krakOut+"/All_classifications.tsv",
                    "--report", krakOut+"/report.txt", "--use-names",
                    "--unclassified-out", krakOut+"/unclassified.fastq",
                    "--classified-out", krakOut+"/classified.fastq", 
                    "--minimum-hit-groups", str(KRAK_THRESH), "--report-minimizer-data")
            print(runCkrak)
            subprocess.call(runCkrak, shell=True)
            print("Kraken complete")
            ClasassemOut = os.path.join(OUT_DIR, "DeNoVo_Assembly/Classified", i)
            UncClassemOut = os.path.join(OUT_DIR, "DeNoVo_Assembly/Unclassified", i)
            ClasassemSt = os.path.join(stats_dir, "DeNoVo_Assembly/Classified", i)
            UncClassemSt = os.path.join(stats_dir, "DeNoVo_Assembly/Unclassified", i)
            if os.path.exists(ClasassemOut):
                pass
            else:
                os.makedirs(ClasassemOut)
            if os.path.exists(UncClassemOut):
                pass
            else:
                os.makedirs(UncClassemOut)
            if os.path.exists(ClasassemSt):
                pass
            else:
                os.makedirs(ClasassemSt)
            if os.path.exists(UncClassemSt):
                pass
            else:
                os.makedirs(UncClassemSt)
            # assemble classified and unclassified separately
            clasFast = os.path.join(krakOut, "classified.fasta")
            clasFastr = os.path.join(krakOut, "classified_formatted.fasta")
            unclasFast = os.path.join(krakOut, "unclassified.fasta")
            unclasFastr = os.path.join(krakOut, "unclassified_formatted.fasta")
            # format the classified and unclassified
            SeqIO.convert(clasFast, "fasta", clasFastr, "fasta")
            SeqIO.convert(unclasFast, "fasta", unclasFastr, "fasta")
            # assemble classified
            run_flye(clasFastr,ClasassemOut,THREADS)
            # assemble unclassified
            run_flye(unclasFastr,UncClassemOut,THREADS)
            print('Assembly completed with metaFlye')
            #### running stats on the raw assembly
            # assembly-stats
            ClassemFile = os.path.join(ClasassemOut, "assembly.fasta")
            UnclassemFile = os.path.join(UncClassemOut, "assembly.fasta")
            ClassemO = ClasassemSt + "/_classified_assembly_stats.txt"
            UnclassemO = UncClassemSt + "/_unclassified_assembly_stats.txt"
            run_AssemStats(i,ClassemFile,ClassemO)
            run_AssemStats(i,UnclassemFile,UnclassemO)
            # quast
            #ClaqAssmO = os.path.join(ClasassemSt, "Quast")
            #UnclaqAssmO = os.path.join(UncClassemSt, "Quast")
            #ClrawQ = raw_Quast(ClassemFile,ClaqAssmO,REFERENCE,REF_GFF,str(THREADS))
            #print('Quast assessment done and saved in ' + rawQ)
            # copy the assembly file to a more accessible point
            print('The generated assembly will be copied to a new folder')
            assemS = os.path.join(OUT_DIR, "Assemblies", "Raw")
            if os.path.exists(assemS):
                pass
            else:
                os.makedirs(assemS)
            ClacoAssemFile = os.path.join(assemS, (i + "_classified.fasta"))
            UnclacoAssemFile = os.path.join(assemS, (i + "_unclassified.fasta"))
            print(ClacoAssemFile)
            print(UnclacoAssemFile)
            shutil.copy2(ClassemFile, ClacoAssemFile)
            shutil.copy2(UnclassemFile, ClacoAssemFile)
            print('The generated assembly files have been copied to ' + assemS)
            # do some clean up at this point
            # check if cleanup option is given
    ####################### if dna is chosen instead ############################
    elif SEQ_TYP == "dna":
        for i in DNA_ISOLATE:
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
            run_AssemStats(i,assemFile,assemO)
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
            shutil.copy2(assemFile, coAssemFile)
            print('The generated assembly files have been copied to ' + assemS)
            # kraken
            cKrak = (KRAK, "--db", KRAKDB, coAssemFile, "--threads", str(THREADS),
                    "--output", krakOut+"/All_classifications.tsv",
                    "--report", krakOut+"/report.txt", "--use-names",
                    "--unclassified-out", krakOut+"/unclassified.fastq",
                    "--classified-out", krakOut+"/classified.fastq", 
                    "--minimum-hit-groups", str(KRAK_THRESH), "--report-minimizer-data")
            runCkrak = ' '.join(cKrak)
            print(runCkrak)
            subprocess.call(runCkrak, shell=True)
            print("Kraken complete")
    ####################### if both seq types are chosen instead ############################
    elif SEQ_TYP == "both":
        for i in DNA_ISOLATE:
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
            #run_flye(toAlgn,assemOut,THREADS)##############################################################
            print('Assembly completed with metaFlye')
            #### running stats on the raw assembly
            # assembly-stats
            assemFile = os.path.join(assemOut, "assembly.fasta")
            assemO = assemSt + "/RAW_assembly_stats.txt"
            run_AssemStats(i,assemFile,assemO)
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
            #shutil.copy2(assemFile, coAssemFile))##############################################################
            # start the classification
            print('Starting kraken classification of generated contigs from DNA')
            krakOut = os.path.join(OUT_DIR, "Kraken", i)
            if os.path.exists(krakOut):
                pass
            else:
                os.makedirs(krakOut)
            # kraken
            cKrak = (KRAK, "--db", KRAKDB, coAssemFile, "--threads", str(THREADS),
                    "--output", krakOut+"/All_classifications.tsv",
                    "--report", krakOut+"/report.txt", "--use-names",
                    "--unclassified-out", krakOut+"/unclassified.fastq",
                    "--classified-out", krakOut+"/classified.fastq", 
                    "--minimum-hit-groups", str(KRAK_THRESH), "--report-minimizer-data")
            runCkrak = ' '.join(cKrak)
            print(runCkrak)
            subprocess.call(runCkrak, shell=True)
            print("Kraken complete")
        # do the cdna processing pipeline --> kraken classification of reads first and then assembly
        for i in CDNA_ISOLATE:
            print('Starting kraken classification of cDNA reads')
            toAlgn = aligned_out + "/" + i + "VsRef_unmapped_renamed.fastq"
            krakOut = os.path.join(OUT_DIR, "Kraken", i)
            if os.path.exists(krakOut):
                pass
            else:
                os.makedirs(krakOut)
            # kraken
            cKrak = (KRAK, "--db", KRAKDB, toAlgn, "--threads", str(THREADS),
                    "--output", krakOut+"/All_classifications.tsv",
                    "--report", krakOut+"/report.txt", "--use-names",
                    "--unclassified-out", krakOut+"/unclassified.fastq",
                    "--classified-out", krakOut+"/classified.fastq", 
                    "--minimum-hit-groups", str(KRAK_THRESH), "--report-minimizer-data")
            runCkrak = ' '.join(cKrak)
            print(runCkrak)
            subprocess.call(runCkrak, shell=True)
            print("Kraken complete")
            ClasassemOut = os.path.join(OUT_DIR, "DeNoVo_Assembly", i , "Classified")
            UncClassemOut = os.path.join(OUT_DIR, "DeNoVo_Assembly", i, "Unclassified")
            ClasassemSt = os.path.join(stats_dir, "DeNoVo_Assembly", i, "Classified")
            UncClassemSt = os.path.join(stats_dir, "DeNoVo_Assembly", i, "Unclassified")
            if os.path.exists(ClasassemOut):
                pass
            else:
                os.makedirs(ClasassemOut)
            if os.path.exists(UncClassemOut):
                pass
            else:
                os.makedirs(UncClassemOut)
            if os.path.exists(ClasassemSt):
                pass
            else:
                os.makedirs(ClasassemSt)
            if os.path.exists(UncClassemSt):
                pass
            else:
                os.makedirs(UncClassemSt)
            # assemble classified and unclassified separately
            clasFast = os.path.join(krakOut, "classified.fastq")
            clasFastr = os.path.join(krakOut, "classified_formatted.fastq")
            unclasFast = os.path.join(krakOut, "unclassified.fastq")
            unclasFastr = os.path.join(krakOut, "unclassified_formatted.fastq")
            # format the classified and unclassified
            #SeqIO.convert(clasFast, "fasta", clasFastr, "fasta")
            #SeqIO.convert(unclasFast, "fasta", unclasFastr, "fasta")
            # assemble classified
            run_flye(clasFast,ClasassemOut,THREADS)
            # assemble unclassified
            run_flye(unclasFast,UncClassemOut,THREADS)
            print('Assembly completed with metaFlye')
            #### running stats on the raw assembly
            # assembly-stats
            ClassemFile = os.path.join(ClasassemOut, "assembly.fasta")
            UnclassemFile = os.path.join(UncClassemOut, "assembly.fasta")
            ClassemO = ClasassemSt + "/_classified_assembly_stats.txt"
            UnclassemO = UncClassemSt + "/_unclassified_assembly_stats.txt"
            run_AssemStats(i,ClassemFile,ClassemO)
            run_AssemStats(i,UnclassemFile,UnclassemO)
            # quast
            #ClaqAssmO = os.path.join(ClasassemSt, "Quast")
            #UnclaqAssmO = os.path.join(UncClassemSt, "Quast")
            #ClrawQ = raw_Quast(ClassemFile,ClaqAssmO,REFERENCE,REF_GFF,str(THREADS))
            #print('Quast assessment done and saved in ' + rawQ)
            # copy the assembly file to a more accessible point
            print('The generated assembly will be copied to a new folder')
            assemS = os.path.join(OUT_DIR, "Assemblies", "Raw")
            if os.path.exists(assemS):
                pass
            else:
                os.makedirs(assemS)
            ClacoAssemFile = os.path.join(assemS, (i + "_classified.fasta"))
            UnclacoAssemFile = os.path.join(assemS, (i + "_unclassified.fasta"))
            print(ClacoAssemFile)
            print(UnclacoAssemFile)
            shutil.copy2(ClassemFile, ClacoAssemFile)
            shutil.copy2(UnclassemFile, ClacoAssemFile)
        print('The generated assembly files have been copied to ' + assemS)
        print('Assembly stages completed')
        print('Remember: cDNA reads were classified and then assembled')
        print('Remember: Both classfied and unclassified cDNA reads were assembled with Flye')
        print('Remember: The DNA genomes have not been classified yet')
    ##### call the racon function for polishing
    # look for the assembly files in the new location
    racAssem = os.path.join(OUT_DIR, "Assemblies", "Racon")
    if os.path.exists(racAssem):
        pass
    else:
        os.makedirs(racAssem)   
    if SEQ_TYP == "dna":
        print('DNA chosen and assemblies will be polished with Racon and Medaka')
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
    if SEQ_TYP == "cdna":
        print('You chose cDNA only, the assemblies will not be polished')
    if SEQ_TYP == "both":
        for i in DNA_ISOLATE:
            readsPath = aligned_out + "/" + i + "VsRef_unmapped_renamed.fastq"
            # set the output folder for each isolate
            racOUT = os.path.join(OUT_DIR, "Racon_Polishing", i)
            print(readsPath)
            print(racOUT)
            # call the racon for 4 iterations
            racStats = runRacon(i,str(f),readsPath,racOUT,THREADS,stats_dir)
            print('The stats for the racon alignments can be found at ' + os.path.dirname(racStats))
            print('Copying the final racon iteraction output to the Assemblies folder')
            racFile = os.path.join(racAssem, (i + "_iteration4.fasta"))
            racc = racOUT + "/" + i + "_iter4.fasta"
            shutil.copy2(racc, racFile)
            print('Racon outputs copied')
            # now call medaka for 1 iteration
            medaOut = os.path.join(OUT_DIR, "Medaka")
            runMedaka(i,readsPath,racFile,medaOut,THREADS)         
    if CLEAN is True:
            print('Cleaning up the DeNoVo Assembly')
            print('They are not deleted, just zipped for archiving')
            denClean = os.path.join(OUT_DIR, "DeNoVo_Assembly")
            denCleanr = denClean + ".tar.gz"
            tarr = 'tar --use-compress-program="pigz --best -p ' + str(THREADS) + '"'
            pgSt = (tarr, "-cf", denCleanr, denClean)
            runPgSt = ' '.join(pgSt)
            print(runPgSt)
            #subprocess.call(runPgSt,shell=True)
            # check if the zipped folder is deleted. If not then delete it ONLY if the zipped folder exists
            if os.path.exists(denClean):
                if os.path.exists(denCleanr):
                    #shutil.rmtree(denClean)
                    print('')
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




print('Pipeline done')
# print(args)
# print(FILE_DIRECTORY)
# print(args.Isolates[1])
# print(OUT_DIR)
# print(Q_KIT)
# print(G_KIT)
# print(FILTER_PASS)
# print(dem_dir)