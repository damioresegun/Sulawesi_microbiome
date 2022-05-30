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
from Scripts.Racon_Medaka import runRacon, runMedaka
from Scripts.Preprocessing import demultip, filt_qc, run_QC
from Scripts.DNA_processing import align
from Scripts.Deduplication import filter_fastq_file
##########################################################################
def get_args():
    parser = argparse.ArgumentParser(description="Pipeline for metagenomics " +
                                     "data pre-processing and analysis for " +
                                     "nanopore sequenced data." +
                                     "Ensure that the tools required are in " +
                                     "in the PATH!", add_help=False)
    file_directory = os.path.realpath(__file__).split("NewNanoMetaPipe")[0]
    if not os.path.isfile(os.path.join(file_directory, "NewNanoMetaPipe.py")):
        file_directory = os.path.join(file_directory, "NewNanoMetaPipe")
    if not os.path.isfile(os.path.join(file_directory, "NewNanoMetaPipe.py")):
        print("Can't locate the correct path to the NewNanoMetaPipe pipeline script")
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
    optional_args.add_argument("-dfl", "--dna_filter_length",
                               dest='dna_Filter_length',
                               action="store",
                               type=int,
                               default="500",
                               help="Filter length threshold. " +
                               "[Default: 500]")
    optional_args.add_argument("-cfl", "--cdna_filter_length",
                               dest='cdna_Filter_length',
                               action="store",
                               type=int,
                               default="100",
                               help="Filter length threshold. " +
                               "[Default: 100]")
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
'''check if the isolate names are properly named'''
# make a list for the DNA and CDNA isolates
DNA_ISOLATE = []
CDNA_ISOLATE = []
for i in ISOLATES:
    iso = i.lower()
    # if the isolate has dna in its name
    if (iso.__contains__("_dna")):
        print('You have provided DNA sequences')
        print(i)
        # add the isolate to the DNA list
        tempd = DNA_ISOLATE
        tempd.append(i)
    # if the isolate has cdna in its name
    elif (iso.__contains__("cdna")):
        print('You have provided cDNA sequences')
        print(i)
        # add the isolate to the cDNA list
        tempc = CDNA_ISOLATE
        tempc.append(i)
    else:
        print('You have not provided the isolates in a satisfactory format')
        print('Do all your isolate names have _dna or _cdna or _dscdna?')
        print('Please look at the help and try again')
        sys.exit(1)
''' check if a cDNA reference genome is given '''
# set the sequence type based on the user's options
if SEQ_TYP == "dna":
    pass
elif SEQ_TYP == "cdna":
    # if the sequence type is cdna, check if the user wants to make a transcriptome
    if MAKCREF is True:
        # if the user wants to make a transcriptome, check if they provide reads
        if not CREADS:
            print('You have chosen to make a transcriptome assembly but provided no reads')
            print('Please provide the reads and run again')
            sys.exit(1)
        else:
            print('You have chosen not to make a transcriptome')
            print(MAKCREF)
            pass
    else:
        # if the user is not making a transcriptome, make the given reference the 
        # cdna transcriptome
        CREF = REFERENCE
        pass
# if the user choose both dna and cda, check inputs
elif SEQ_TYP == "both":
    # check if the dna reference is given
    if not REFERENCE:
        print('Did you provide the two needed references?')
        print('Please remember a DNA reference and a transcriptome reference '+
                'are required')
        print('Please try again')
        sys.exit(1)
    # check if the cdna reference is given
    elif not CREF:
        # check if the transcriptome assembly is to be made
        if MAKCREF is True:
            if not CREADS:
                print('You have chosen to make a transcriptome assembly but provided no reads')
                print('Please provide the reads and run again')
                sys.exit(1)
            else:
                print('You are giving both DNA and cDNA reads and have selected to ' +
                'generate a transcriptome assembly.')
            # check if the adaptors are given for the cdna
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
'''set some variables'''
# check the demultiplexer choice
# if its qcat
if DEMULP_CHOICE == "qcat":
    G_KIT = args.Sequencing_kit
    # set the correct 
    Q_KIT = "NBD103/" + args.Expansion_kit
# if its guppy
else:
    G_KIT = "SQK-" + args.Sequencing_kit
    Q_KIT = "EXP-" + args.Expansion_kit
# check if the user chose to carry out filtering
if FILTER_PASS is True:
    # set the dna filter length
    DNA_FILT_LENGTH = args.dna_Filter_length
    CDNA_FILT_LENGTH = args.cdna_Filter_length
    # set the quality
    FILT_QUAL = args.Filter_Quality
    # set bracken length to the filter length
    DBRACK_LENGTH = DNA_FILT_LENGTH
    CBRACK_LENGTH = DNA_FILT_LENGTH
else:
    pass
# set bracken length to the default filter length
    DBRACK_LENGTH = DNA_FILT_LENGTH
    CBRACK_LENGTH = DNA_FILT_LENGTH
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
# make/check folder to hold demultiplexed reads
dem_dir = os.path.join(OUT_DIR, "Demultiplexed")
if os.path.exists(dem_dir):
    pass
else:
    os.makedirs(dem_dir)
# make/check folder to hold stats outputs
stats_dir = os.path.join(OUT_DIR, "Stats")
if os.path.exists(stats_dir):
    pass
else:
    os.mkdir(stats_dir)
# make/check folder to hold temporary alignment files
# folder will be deleted if the user chooses to clean up
temp_align_out = os.path.join(OUT_DIR, "Isolate_Reads_Aligned_Vs_Reference")
if os.path.exists(temp_align_out):
    pass
else:
    os.makedirs(temp_align_out)
# make/check folder to hold cleaned reads
aligned_out = os.path.join(OUT_DIR, "Host_Free_Reads")
if os.path.exists(aligned_out):
    pass
else:
    os.makedirs(aligned_out)
# if the user does not want to filter
# use the reads as is
if FILTER_PASS is False:
    rnm_path = os.path.join(OUT_DIR, "Isolate_Demultiplexed_Reads")
    if os.path.exists(rnm_path):
        pass
    else:
        os.makedirs(rnm_path)
##########################################################################
''' Filtering reads
 Reads will be filtered.  '''
##########################################################################
# define the function
def dna_filter():
    # set isolate to the indexed DNA_ISOLATE
    count = 0
    isola = DNA_ISOLATE[count]
    ''' run the filt_qc function that is in the 'Preprocessing.py' script
    and set it to a variable '''
    ready_path = filt_qc(dem_dir, BARCODES[count], isola,
                                        OUT_DIR,DNA_FILT_LENGTH,FILT_QUAL,)
    print('The raw demultiplexed reads have been successfully filtered ' + 
        'and saved in ' + ready_path)
    print('Please remember that the files are now renamed')
    print(BARCODES[count] + ' is now ' + isola)
    print('')
    print('Now running QC of filtered reads')
    # make a txt file with the barcode for filename
    ofile = DNA_ISOLATE[count] + ".txt"
    dem_file = ready_path + "/" + isola + ".fastq.gz"
    # set a temporary filename
    temp = BARCODES[count] + "_" + isola
    stats = os.path.join(stats_dir, "Filtered_Demultiplexed_Reads", temp)
    # run the run_QC function on the filtered reads
    run_QC(dem_file,BARCODES[count],stats,ofile,THREADS)
    count += 1
    return ready_path
def cdna_filter():
    # set isolate to the indexed CDNA_ISOLATE
    count = 0
    cisola = CDNA_ISOLATE[count]
    ''' run the filt_qc function that is in the 'Preprocessing.py' script
    and set it to a variable '''
    ready_path = filt_qc(dem_dir, BARCODES[count], cisola,
                                        OUT_DIR,CDNA_FILT_LENGTH,FILT_QUAL,)
    print('The raw demultiplexed reads have been successfully filtered ' + 
        'and saved in ' + ready_path)
    print('Please remember that the files are now renamed')
    print(BARCODES[count] + ' is now ' + cisola)
    print('')
    print('Now running QC of filtered reads')
    # make a txt file with the barcode for filename
    ofile = CDNA_ISOLATE[count] + ".txt"
    dem_file = ready_path + "/" + cisola + ".fastq.gz"
    # set a temporary filename
    temp = BARCODES[count] + "_" + cisola
    stats = os.path.join(stats_dir, "Filtered_Demultiplexed_Reads", temp)
    # run the run_QC function on the filtered reads
    run_QC(dem_file,BARCODES[count],stats,ofile,THREADS)
    count += 1
    return ready_path
##########################################################################
''' Assemble reads using metaFlye
 Reads which have been extracted from alignments against the reference
 genome will go through de novo genome assembly using the '--meta' flag
 of flye which allows for higher error reads and also allows for
 variance in coverage. Duplicated reads will need to be checked and
 renamed. This is done using a custom python script that checks the name
 of a read and renames it if it is duplicated. Then metaFlye takes place '''
##########################################################################
# define the flye function
def run_flye(in_file,out_dir,threads):
    try:
        # set the command
        fly = ("flye --nano-raw", in_file, "--out-dir", out_dir, "--threads", str(threads), "--meta")
        # join the string
        runFly = ' '.join(fly)
        print(runFly)
        # run the command
        subprocess.call(runFly, shell=True)
        # catch errors
    except FileNotFoundError:
        print('There was an issue with the file. File not found')
    except SyntaxError:
        print('There was a syntax issue. Double check your syntax')
##########################################################################
''' Alignment functions. Two functions for the DNA and cDNA read option 
DNA_align() carries out alignment using the list of DNA isolates. It makes
    a stats directory to carry out some descriptive stats. It also carries
    out deduplication using a function from an external script.
cDNA_align() acts in two forms, depending on the user's wishes. If the
    user chose to not create a transcriptome and has provided one already
    then the cDNA reads are aligned against the transcriptome. If the user
    chose to make a transcriptome, then a de novo transcriptome assembly
    process is started. For this the user needs to have provided reads to
    generated the transcriptome and other options like max memory. Once
    finished, deduplication occurs and then moved to a central folder.'''
##########################################################################
def DNA_align(ready_path):
    for i in DNA_ISOLATE:
            # set the right filtered/demultiplexed file
            file = ready_path + "/" + i + ".fastq.gz"
            # set/make the stats directory
            statss = os.path.join(stats_dir, "Alignment_Vs_Reference", i)
            if os.path.exists(statss):
                pass
            else:
                os.makedirs(statss)
            # run the alignment of DNA against the reference
            hs_reads = align(i,file,temp_align_out,aligned_out,THREADS,statss,REFERENCE)
            # Run de-duplication
            print('Running de-deuplication on ' + i)
            alnRename = aligned_out + "/" + i + "_DNAVsRef_unmapped_renamed.fastq"
            ''' filter_fastq_file function carries out deduplication and is in the
            'Deduplication.py' script '''
            filter_fastq_file(hs_reads,alnRename)
            print('De-duplication complete for ' + i)
            print("The reads have been renamed and saved as: " + alnRename)
def cDNA_align(ready_path):
    for i in CDNA_ISOLATE:
            # set the right filtered/demultiplexed file
            file = ready_path + "/" + i + ".fastq.gz"
            # make/set the cDNA_processing folder
            cdOut = os.path.join(OUT_DIR, "cDNA_Processing")
            if os.path.exists(cdOut):
                pass
            else:
                os.makedirs(cdOut)
            # set the path to the cdna processing script
            cdProc = os.path.join(SCPTS, "cDNA_Processing.py")
            # check if the user chose to make a transcriptome assembly
            if MAKCREF is True:
                # make the command 
                cdProcy = (cdProc, "-s reads -o", cdOut, "-r", CREADS[0], CREADS[1], 
                            "-a", CADAP, "-hr", REFERENCE, "-hg", REF_GFF, 
                            "-p", str(THREADS), "-m", MXMEM, "-ur", file, "2>&1 | tee bothlogger.txt")
                # join them together
                runCdProcy = ' '.join(cdProcy)
                print(runCdProcy)
                # run the command
                subprocess.call(runCdProcy, shell=True)
                print('cDNA transcriptome generated')
            else:
                ''' if the user does not want to make a transcriptome, then align their reads
                against the given transcriptome using the cdna_processing script '''
                cdProcy = (cdProc, "-s genome -o", cdOut, "-p", str(THREADS), "-ur", file, 
                            "-t", CREF)
                runCdProcy = ' '.join(cdProcy)
                print(runCdProcy)
                # run the command
                subprocess.call(runCdProcy, shell=True)
            # get just the isolate name without the 'DNA or cDNA'
            iaw = i.split("_")[0]
            # copy the fastq of unmapped sequences to the host_free_reads folders
            fileOut = cdOut + "/" + "Alignment" + "/" + iaw + "_cDNA_VsTranscriptome.fastq"
            hfFileOt = aligned_out + "/" + iaw + "_cDNAVsRef_unmapped.fastq"
            # copy and rename in the new location
            shutil.copy2(fileOut, hfFileOt)
            print('cDNA reads aligned against transcriptome and host reads separated')
            print('Host-free cDNA reads are now in ' + aligned_out)
            ''' filter_fastq_file function carries out deduplication and is in the
            'Deduplication.py' script '''
            print('Running de-deuplication on ' + i)
            # set the path for the deduplicated fastq
            hfRnm = aligned_out + "/" + iaw + "_cDNAVsRef_unmapped_renamed.fastq"
            filter_fastq_file(hfFileOt,hfRnm)
            print('De-duplication complete for ' + i)
            print("The reads have been renamed and saved as: " + hfRnm)
####################################################################################################################################################
'''Kraken classification'''
#def dna_Kraken():
    
def cdna_Kraken(isolate,seq):
    # set the path to the deduplicated fastq- of the isolate given
    toAlgn = aligned_out + "/" + isolate + "_" + seq +"VsRef_unmapped_renamed.fastq"
    # outputs of this function will be saved in a folder called 'Kraken'
    krakOut = os.path.join(OUT_DIR, "Kraken", isolate)
    # check if the kraken folder already exists
    if os.path.exists(krakOut):
        pass
    else:
        os.makedirs(krakOut)
    # build the kraken command
    cKrak = (KRAK, "--db", KRAKDB, toAlgn, "--threads", str(THREADS),
            "--output", krakOut+"/All_classifications.tsv",
            "--report", krakOut+"/report.txt", "--use-names",
            "--unclassified-out", krakOut+"/unclassified.fastq",
            "--classified-out", krakOut+"/classified.fastq", 
            "--minimum-hit-groups", str(KRAK_THRESH), "--report-minimizer-data")
    runCkrak = ' '.join(cKrak)
    print(runCkrak)
    # run the kraken command
    subprocess.call(runCkrak, shell=True)
    print("Kraken complete")
    # return the path to the isolate kraken folder
    return krakOut
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
    ##########################################################################
    # Call functions
    ##########################################################################
    '''Demultiplexing. Function: demultip is in the Preprocessing.py script'''
    # check if the user just wants to re-try demultiplexing
    if REDEMULP is True:
        print("You are just re-generating the demultiplexed reads")
        # run the demultiplexing function from 'Preprocessing.py' script
        demultip(INP_DIR, dem_dir, DEMULP_CHOICE, THREADS, Q_KIT)
        # finish demultiplexing and exit
        print("Demultiplexed reads regenerated")
        sys.exit(1)
    else:
        # run demultiplexing function
        demultip(INP_DIR, dem_dir, DEMULP_CHOICE, THREADS, Q_KIT)
    ''' Zip the demultiplexed reads '''
    print('Zipping the fastq files')
    # zip each file in the demultiplexed file
    for i in os.listdir(dem_dir):
        zip_dem = os.path.join(dem_dir, i)
        pgSt = ("pigz --best", zip_dem)
        runPgSt = ' '.join(pgSt)
        print(runPgSt)
        # run the command
        subprocess.call(runPgSt,shell=True)
    print('FastQ files have been zipped')
    #############################################################################
    ''' qc the raw demultiplexed reads '''
    #############################################################################
    # for each barcode...
    for i in BARCODES:
        # make a txt file with the barcode for filename
        ofile = i + ".txt"
        # select the zipped demultiplexed file
        dem_file = dem_dir + "/" + i + ".fastq.gz"
        # set the stats path
        stats = os.path.join(stats_dir, "Raw_Demultiplexed_Reads", i)
        # run the runQC function that is in the 'Preprocessing.py' script
        run_QC(dem_file,i,stats,ofile,THREADS)
    #############################################################################
    ''' filtering: if filter check is true '''
    #############################################################################
    if FILTER_PASS is True:
        #count = 0
        # if the user sequence type is just DNA
        if SEQ_TYP == "dna":
            ready_path = dna_filter()
        # if the user sequence type is just cdna 
        elif SEQ_TYP == "cdna":
            ready_path = cdna_filter()
        # if both sequence types are to be used
        elif SEQ_TYP == "both":
            ''' dna filtering first '''
            ready_path = dna_filter()
            """ # set isolate to the indexed DNA_ISOLATE
            isola = DNA_ISOLATE[count]
            ''' run the filt_qc function that is in the 'Preprocessing.py' script
            and set it to a variable '''
            ready_path = filt_qc(dem_dir, BARCODES[count], isola,
                                                OUT_DIR,DNA_FILT_LENGTH,FILT_QUAL,)
            print('The raw demultiplexed reads have been successfully filtered ' + 
                'and saved in ' + ready_path)
            print('Please remember that the files are now renamed')
            print(BARCODES[count] + ' is now ' + isola)
            print('')
            print('Now running QC of filtered reads')
            ofile = DNA_ISOLATE[count] + ".txt"
            dem_file = ready_path + "/" + isola + ".fastq.gz"
            temp = BARCODES[count] + "_" + isola
            stats = os.path.join(stats_dir, "Filtered_Demultiplexed_Reads", temp)
            run_QC(dem_file,BARCODES[count],stats,ofile,THREADS) """
            # cdna filtering
            ready_path = cdna_filter()
            """ cisola = CDNA_ISOLATE[count]
            ready_path = filt_qc(dem_dir, BARCODES[count], cisola,
                                                OUT_DIR,CDNA_FILT_LENGTH,FILT_QUAL,)
            print('The raw demultiplexed reads have been successfully filtered ' + 
                'and saved in ' + ready_path)
            print('Please remember that the files are now renamed')
            print(BARCODES[count] + ' is now ' + cisola)
            print('')
            print('Now running QC of filtered reads')
            # make a txt file with the barcode for filename
            ofile = CDNA_ISOLATE[count] + ".txt"
            dem_file = ready_path + "/" + cisola + ".fastq.gz"
            # set a temporary filename
            temp = BARCODES[count] + "_" + cisola
            stats = os.path.join(stats_dir, "Filtered_Demultiplexed_Reads", temp)
            # run the run_QC function on the filtered reads
            run_QC(dem_file,BARCODES[count],stats,ofile,THREADS)
            count += 1 """
        '''if reads are not being filtered then move the demultiplexed file'''
    else:
        print('You chosen to not filter your reads.')
        print('Your demultiplexed reads will be used to proceed')
        ready_path = rnm_path
        count = 0
        for i in BARCODES:
            isola = ISOLATES[count]
            dem_file = dem_dir + "/" + i + ".fastq.gz"
            rnm_file = ready_path + "/" + isola + ".fastq.gz"
            # set the moving command
            mv_stuff = ("mv", dem_file, rnm_file)
            rnMv_Stff = ' '.join(mv_stuff)
            print(rnMv_Stff)
            # run the command
            subprocess.call(rnMv_Stff, shell=True)
            count += 1
    ###########################################################################################
    ''' align against reference genome and de-duplicate
     first check if the sequence type is dna or cdna
     align function is in DNA_processing.py'''
    ###########################################################################################
     # if dna is the sequence type
    if SEQ_TYP == "dna":
        # use the DNA isolate list
        print('Aligning reads against the reference')
        DNA_align(ready_path)
        print('Alignment done')
    elif SEQ_TYP == "cdna":
        print('Aligning reads against the reference')
        cDNA_align(ready_path)
        print('Alignment done')
    elif SEQ_TYP == "both":
        # do dna alignment
        print('Aligning reads against the DNA reference')
        DNA_align(ready_path)
        print('Alignment done')
        # do cdna alignment
        print('Aligning reads against the transcriptome reference')
        cDNA_align(ready_path)
        print('Alignment done')
    ###########################################################################################
    '''Kraken and assembly classification'''
    ###########################################################################################
    # if the user chose to use cdna reads
    if SEQ_TYP == "cdna":
        for i in CDNA_ISOLATE:
            # run the kraken classification function
            sequt = "cDNA"
            iaw = i.split("_")[0]
            ckrakOut = cdna_Kraken(iaw, sequt)
            print('Kraken completed. Formatting outputs')
            # set the paths to the outputs of kraken and the newly formatted versions
            clasFast = os.path.join(ckrakOut, "classified.fastq")
            unclasFast = os.path.join(ckrakOut, "unclassified.fastq")
            ''' de novo assembly with metaflye'''
            print('Starting de novo assembly')
            # set the output locations for metaflye
            ClasassemOut = os.path.join(OUT_DIR, "DeNoVo_Assembly", i, "Classified")
            UncClassemOut = os.path.join(OUT_DIR, "DeNoVo_Assembly", i, "Unclassified")
            if os.path.exists(ClasassemOut):
                pass
            else:
                os.makedirs(ClasassemOut)
            if os.path.exists(UncClassemOut):
                pass
            else:
                os.makedirs(UncClassemOut)
            # assemble classified
            run_flye(clasFast,ClasassemOut,THREADS)
            # assemble unclassified
            run_flye(unclasFast,UncClassemOut,THREADS)
            print('Assembly completed with metaFlye')
             # running stats on the raw assemblies
            print('Now running general stats on the generated assemblies')
            # make the folders to hold stats for the assembled genomes
            ClasassemSt = os.path.join(stats_dir, "DeNoVo_Assembly", i, "Classified")
            UncClassemSt = os.path.join(stats_dir, "DeNoVo_Assembly", i, "Unclassified")
            if os.path.exists(ClasassemSt):
                pass
            else:
                os.makedirs(ClasassemSt)
            if os.path.exists(UncClassemSt):
                pass
            else:
                os.makedirs(UncClassemSt)
            # assembly-stats
            ClassemFile = os.path.join(ClasassemOut, "assembly.fasta")
            UnclassemFile = os.path.join(UncClassemOut, "assembly.fasta")
            ClassemO = ClasassemSt + "/_classified_assembly_stats.txt"
            UnclassemO = UncClassemSt + "/_unclassified_assembly_stats.txt"
            # run the assembly stats function that is in the AssemblyQC script
            run_AssemStats(i,ClassemFile,ClassemO)
            run_AssemStats(i,UnclassemFile,UnclassemO)
            print('The generated assembly will be copied to a new folder')
            # set/make the folder to transfer the raw assemblies for easier access
            assemS = os.path.join(OUT_DIR, "Assemblies", "Raw")
            if os.path.exists(assemS):
                pass
            else:
                os.makedirs(assemS)
            # set/make the 'Assemblies' folder to only hold the assembled genomes
            ClacoAssemFile = os.path.join(assemS, (i + "_classified.fasta"))
            UnclacoAssemFile = os.path.join(assemS, (i + "_unclassified.fasta"))
            # copy the assemblies to the new location
            shutil.copy2(ClassemFile, ClacoAssemFile)
            shutil.copy2(UnclassemFile, ClacoAssemFile)
            print('The generated assembly files have been copied to ' + assemS)
    elif SEQ_TYP == "dna":
        # run the classification function
        for i in DNA_ISOLATE:
            sequt = "DNA"
            iaw = i.split("_")[0]
            ckrakOut = cdna_Kraken(iaw, sequt)
            print('Kraken completed. Formatting outputs')
            # set the paths to the outputs of kraken and the newly formatted versions
            clasFast = os.path.join(ckrakOut, "classified.fastq")
            unclasFast = os.path.join(ckrakOut, "unclassified.fastq")
            ''' de novo assembly with metaflye'''
            print('Starting de novo assembly')
            # set the output locations for metaflye
            ClasassemOut = os.path.join(OUT_DIR, "DeNoVo_Assembly", i, "Classified")
            UncClassemOut = os.path.join(OUT_DIR, "DeNoVo_Assembly", i, "Unclassified")
            if os.path.exists(ClasassemOut):
                pass
            else:
                os.makedirs(ClasassemOut)
            if os.path.exists(UncClassemOut):
                pass
            else:
                os.makedirs(UncClassemOut)
            # assemble classified
            run_flye(clasFast,ClasassemOut,THREADS)
            # assemble unclassified
            run_flye(unclasFast,UncClassemOut,THREADS)
            print('Assembly completed with metaFlye')
             # running stats on the raw assemblies
            print('Now running general stats on the generated assemblies')
            # make the folders to hold stats for the assembled genomes
            ClasassemSt = os.path.join(stats_dir, "DeNoVo_Assembly", i, "Classified")
            UncClassemSt = os.path.join(stats_dir, "DeNoVo_Assembly", i, "Unclassified")
            if os.path.exists(ClasassemSt):
                pass
            else:
                os.makedirs(ClasassemSt)
            if os.path.exists(UncClassemSt):
                pass
            else:
                os.makedirs(UncClassemSt)
            # assembly-stats
            ClassemFile = os.path.join(ClasassemOut, "assembly.fasta")
            UnclassemFile = os.path.join(UncClassemOut, "assembly.fasta")
            ClassemO = ClasassemSt + "/_classified_assembly_stats.txt"
            UnclassemO = UncClassemSt + "/_unclassified_assembly_stats.txt"
            # run the assembly stats function that is in the AssemblyQC script
            run_AssemStats(i,ClassemFile,ClassemO)
            run_AssemStats(i,UnclassemFile,UnclassemO)
            print('The generated assembly will be copied to a new folder')
            # set/make the folder to transfer the raw assemblies for easier access
            assemS = os.path.join(OUT_DIR, "Assemblies", "Raw")
            if os.path.exists(assemS):
                pass
            else:
                os.makedirs(assemS)
            # set/make the 'Assemblies' folder to only hold the assembled genomes
            ClacoAssemFile = os.path.join(assemS, (i + "_classified.fasta"))
            UnclacoAssemFile = os.path.join(assemS, (i + "_unclassified.fasta"))
            # copy the assemblies to the new location
            shutil.copy2(ClassemFile, ClacoAssemFile)
            shutil.copy2(UnclassemFile, ClacoAssemFile)
            print('The generated assembly files have been copied to ' + assemS)
    elif SEQ_TYP == "both":
        for i in CDNA_ISOLATE:
            # run the kraken classification function
            sequt = "cDNA"
            iaw = i.split("_")[0]
            ckrakOut = cdna_Kraken(iaw, sequt)
            print('Kraken completed. Formatting outputs')
            # set the paths to the outputs of kraken and the newly formatted versions
            clasFast = os.path.join(ckrakOut, "classified.fastq")
            unclasFast = os.path.join(ckrakOut, "unclassified.fastq")
            ''' de novo assembly with metaflye'''
            # sometimes cdna does not classify so check if this occurred
            if os.path.exists(clasFast):
                print('Starting de novo assembly')
                # set the output locations for metaflye
                ClasassemOut = os.path.join(OUT_DIR, "DeNoVo_Assembly", i, "Classified")
                UncClassemOut = os.path.join(OUT_DIR, "DeNoVo_Assembly", i, "Unclassified")
                if os.path.exists(ClasassemOut):
                    pass
                else:
                    os.makedirs(ClasassemOut)
                if os.path.exists(UncClassemOut):
                    pass
                else:
                    os.makedirs(UncClassemOut)
                # assemble classified
                run_flye(clasFast,ClasassemOut,THREADS)
                # assemble unclassified
                run_flye(unclasFast,UncClassemOut,THREADS)
                print('Assembly completed with metaFlye')
                # running stats on the raw assemblies
                print('Now running general stats on the generated assemblies')
                # make the folders to hold stats for the assembled genomes
                ClasassemSt = os.path.join(stats_dir, "DeNoVo_Assembly", i, "Classified")
                UncClassemSt = os.path.join(stats_dir, "DeNoVo_Assembly", i, "Unclassified")
                if os.path.exists(ClasassemSt):
                    pass
                else:
                    os.makedirs(ClasassemSt)
                if os.path.exists(UncClassemSt):
                    pass
                else:
                    os.makedirs(UncClassemSt)
                # assembly-stats
                ClassemFile = os.path.join(ClasassemOut, "assembly.fasta")
                UnclassemFile = os.path.join(UncClassemOut, "assembly.fasta")
                ClassemO = ClasassemSt + "/_classified_assembly_stats.txt"
                UnclassemO = UncClassemSt + "/_unclassified_assembly_stats.txt"
                # run the assembly stats function that is in the AssemblyQC script
                run_AssemStats(i,ClassemFile,ClassemO)
                run_AssemStats(i,UnclassemFile,UnclassemO)
                print('The generated assembly will be copied to a new folder')
                # set/make the folder to transfer the raw assemblies for easier access
                assemS = os.path.join(OUT_DIR, "Assemblies", "Raw")
                if os.path.exists(assemS):
                    pass
                else:
                    os.makedirs(assemS)
                # set/make the 'Assemblies' folder to only hold the assembled genomes
                ClacoAssemFile = os.path.join(assemS, (i + "_classified.fasta"))
                UnclacoAssemFile = os.path.join(assemS, (i + "_unclassified.fasta"))
                # copy the assemblies to the new location
                shutil.copy2(ClassemFile, ClacoAssemFile)
                shutil.copy2(UnclassemFile, ClacoAssemFile)
                print('The generated assembly files have been copied to ' + assemS)
            else:
                print('cDNA classification did not work properly.')
                print('This will likely break things somewhere. Watch out')
        for i in DNA_ISOLATE:
            sequt = "DNA"
            iaw = i.split("_")[0]
            ckrakOut = cdna_Kraken(iaw, sequt)
            print('Kraken completed. Formatting outputs')
            # set the paths to the outputs of kraken and the newly formatted versions
            clasFast = os.path.join(ckrakOut, "classified.fastq")
            unclasFast = os.path.join(ckrakOut, "unclassified.fastq")
            ''' de novo assembly with metaflye'''
            print('Starting de novo assembly')
            # set the output locations for metaflye
            ClasassemOut = os.path.join(OUT_DIR, "DeNoVo_Assembly", i, "Classified")
            UncClassemOut = os.path.join(OUT_DIR, "DeNoVo_Assembly", i, "Unclassified")
            if os.path.exists(ClasassemOut):
                pass
            else:
                os.makedirs(ClasassemOut)
            if os.path.exists(UncClassemOut):
                pass
            else:
                os.makedirs(UncClassemOut)
            # assemble classified
            run_flye(clasFast,ClasassemOut,THREADS)
            # assemble unclassified
            run_flye(unclasFast,UncClassemOut,THREADS)
            print('Assembly completed with metaFlye')
             # running stats on the raw assemblies
            print('Now running general stats on the generated assemblies')
            # make the folders to hold stats for the assembled genomes
            ClasassemSt = os.path.join(stats_dir, "DeNoVo_Assembly", i, "Classified")
            UncClassemSt = os.path.join(stats_dir, "DeNoVo_Assembly", i, "Unclassified")
            if os.path.exists(ClasassemSt):
                pass
            else:
                os.makedirs(ClasassemSt)
            if os.path.exists(UncClassemSt):
                pass
            else:
                os.makedirs(UncClassemSt)
            # assembly-stats
            ClassemFile = os.path.join(ClasassemOut, "assembly.fasta")
            UnclassemFile = os.path.join(UncClassemOut, "assembly.fasta")
            ClassemO = ClasassemSt + "/_classified_assembly_stats.txt"
            UnclassemO = UncClassemSt + "/_unclassified_assembly_stats.txt"
            # run the assembly stats function that is in the AssemblyQC script
            run_AssemStats(i,ClassemFile,ClassemO)
            run_AssemStats(i,UnclassemFile,UnclassemO)
            print('The generated assembly will be copied to a new folder')
            # set/make the folder to transfer the raw assemblies for easier access
            assemS = os.path.join(OUT_DIR, "Assemblies", "Raw")
            if os.path.exists(assemS):
                pass
            else:
                os.makedirs(assemS)
            # set/make the 'Assemblies' folder to only hold the assembled genomes
            ClacoAssemFile = os.path.join(assemS, (i + "_classified.fasta"))
            UnclacoAssemFile = os.path.join(assemS, (i + "_unclassified.fasta"))
            # copy the assemblies to the new location
            shutil.copy2(ClassemFile, ClacoAssemFile)
            shutil.copy2(UnclassemFile, ClacoAssemFile)
            print('The generated assembly files have been copied to ' + assemS)
