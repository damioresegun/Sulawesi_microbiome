#!/usr/bin/env python3
'''
NanoMetaPipe
Author: Damilola R Oresegun
Contributors: Peter Thorpe

A pipeline to take basecalled metagenomic nanopore reads through preprocessing
to final assembled bins for microbiome analysis. 
The user will provide either raw reads or basecalled reads with the path to 
the output folder to place generated outputs. Other options include different
thresholds necessary to trimming, filtering and classification steps.
'''
#  import packages
import os
import sys
import configparser
import subprocess
import argparse
import logging
import logging.handlers
import time
import shutil
from pathlib import Path
#from PipelineDevScripts.NanoMetaPipe_assemblyApproach import CDNA_ISOLATE, DNA_ISOLATE
from Scripts.PreChecks import filterOptions, isolateList, seqCheck
from Scripts.Racon_Medaka import runRacon, runMedaka
from Scripts.Preprocessing import cdna_filter, demultip, dna_filter, filt_qc, run_QC
from Scripts.DNA_processing import DNA_align, align, cDNA_align
from Scripts.Tools import krakBrak, run_flye, zipFiles, makeDirectory, run_AssemStats, raw_Quast
#############################################################################################
def get_args():
    '''
    Function to take in the different arguments necessary to run the pipeline script. 
    No inputs are needed but its returns a list of argument variables and the file directory
    where the script is saved to use as the working directory.
    '''
    parser = argparse.ArgumentParser(description = "A pipeline to take basecalled " + 
                "metagenomic nanopore reads through preprocessing to final " + 
                "assembled bins for microbiome analysis. The user will provide " + 
                "either raw reads or basecalled reads with the path to the output " +
                "folder to place generated outputs. Other options include different " +
                "thresholds necessary to trimming, filtering and classification steps.",
                add_help = False)
    file_directory = os.path.realpath(__file__).split("NanoMetaPipe")[0]
    if not os.path.isfile(os.path.join(file_directory, "NanoMetaPipe.py")):
        file_directory = os.path.join(file_directory, "NanoMetaPipe")
    if not os.path.isfile(os.path.join(file_directory, "NanoMetaPipe.py")):
        print("Can't locate the correct path to the NanoMetaPipe pipeline script")
    #########################################################################################
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument("-b", "--basecalled", 
                                dest = "Basecalled_Reads", 
                                action = "store", type = str,
                                help = "Folder of basecalled reads in FASTQ format. " +
                                "Reads must be basecalled previously!. " +
                                "Example: path/to/exp_folder/pass. This option does not " +
                                "work with raw reads. ", required = True)
    required_args.add_argument("-c", "--barcodes", 
                                dest = "Barcodes",
                                action = "store", nargs = "+",
                                help = "List of barcodes used e.g. -c barcode01 barcode02",
                                required = True)
    required_args.add_argument("-e", "--isolates", 
                                dest = "Isolates",
                                action = "store", nargs = "+",
                                help = "List of isolates used. Must give sequence type " +
                                "separated by underscore e.g. -e isolate1_dna isolate2_cdna. " +
                                "Must correspond to the order of barcodes given i.e. " +
                                "barcode01=isolate1_dna.", required = True)
    required_args.add_argument("-kr", "--kraken", 
                                dest = "Kraken_PATH",
                                action = "store", type = str, 
                                help = "Full path to the install kraken package. " +
                                "If kraken is in the $PATH, simply write kraken2",
                                required = True)
    required_args.add_argument("-kb", "--kraken_DB", 
                                dest = "Kraken_DB_Path",
                                action = "store", type = str,
                                help = "Full path to the installed kraken database",
                                required = True)
    required_args.add_argument("-r", "--reference", 
                                dest = "Reference_Genome",
                                action = "store", type = str,
                                help = "Full path to the reference genome to align against.",
                                required = True)
    required_args.add_argument("-s", "--sequence_type", 
                                dest = "Sequence_Type",
                                action = "store", 
                                choices = ["dna", "cdna", "both"],
                                type = str, 
                                help = "The type of sequence you want to be " +
                                "analysed. DNA and cDNA are supported individually and " +
                                "together in a multiplexed dataset. Ensure to identify " +
                                "the sequence types in the isolate names given i.e. " +
                                "Sequence files should be given as isolate1_dna or  "+
                                "isolate2_cdna or isolate2_dscdna.", required = True)
    required_args.add_argument("-br", "--bracken", 
                                dest = "Bracken_PATH",
                                action = "store", type = str, 
                                help = "Full path to the installed bracken. " +
                                "If in the $PATH, simply write bracken", required = True)
    required_args.add_argument("-nd", "--ncbi_db",
                                dest = "NCBI_DB",
                                action = "store", type = str,
                                help = "The full path to the local NCBI nt database", 
                                required = True)
    #########################################################################################
    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument("-bt", "--bracken-hit-threshold", 
                                dest = "Bracken_Hit_Threshold", 
                                action = "store", type = str, default = 20,
                                help = "A minimum number of kmers that must be matched to " +
                                "place a contig into a taxonomic group by bracken " +
                                "re-estimation. Default is [20]")
    optional_args.add_argument("-cr", "--cdna_ref", 
                                dest = "cDNA_Reference",
                                action = "store", type = str, 
                                help = "If using both DNA and CDNA sequence types, please " +
                                "provide the path to the transcriptome assembly here. " +
                                "If only using cDNA sequence files, this option is not " +
                                "necessary!")
    optional_args.add_argument("-mcr", "--make_cdna_ref", 
                                dest = "make_cDNA_Reference",
                                action = "store_true", 
                                help = "Generate a genome-guided transcriptome to remove " +
                                "host sequences from input cDNA sequences. Only works with " +
                                "with the detection of indicated cDNA input isolates. " +
                                "This option does not work with the '-cr' option. " +
                                "[Default:Off]. Activate this option with just '-mcr' or " +
                                "'--make_cdna_ref'. To be used in conjunction with the " +
                                "'-cd' option")
    optional_args.add_argument("-cd", "--ref_cdna_reads", 
                                dest = "reference_cDNA_Reads",
                                action = "store", nargs = "+", 
                                help = "Path to cDNA reads to generate a host " +
                                "transcriptome. Must be used with '-mcr' option")
    optional_args.add_argument("-ca", "--cdna_adapter", 
                                dest = "cDNA_Adapters", 
                                action = "store", type = str,
                                help = "Full path to a FASTA file containing the " +
                                "adapters used for your cDNA sequencing reads. Must be " +
                                "used with the '-mcr' option")
    optional_args.add_argument("-d", "--demultiplexer", 
                                dest = "Demultiplexer",
                                action = "store", type = str, 
                                default = "qcat", choices = ["qcat", "guppy"],
                                help = "Choice of demultiplexer to use. Guppy or Qcat. " +
                                "If you wish to use guppy, you must have access to a GPU " +
                                "set up with tensorflow and other requirements set by " +
                                "Oxford Nanopore. [Default: qcat]")
    optional_args.add_argument("-dfl", "--dna_filter_length", 
                                dest = "dna_Filter_Length",
                                action = "store", type = int, default = "500", 
                                help = "Filter length threshold for DNA input sequences " +
                                "[Default: 500]")
    optional_args.add_argument("-cfl", "--cdna_filter_length", 
                                dest = "cDNA_Filter_Length",
                                action = "store", type = int, default = "100",
                                help = "Filter length threshold for cDNA input sequences " +
                                "[Default: 100]")
    optional_args.add_argument("-flw", "--flowcell", 
                                dest = "Flowcell",
                                action = "store", type = str, default = "FLO-MIN106",
                                help = "The flowcell used for the sequencing experiment. " +
                                "[Default: FLO-MIN106]")
    optional_args.add_argument("-fq", "--filter_quality", 
                                dest = "Filter_Quality",
                                action = "store", type = int, default = "10",
                                help = "Phred quality threshold to filter demultiplexed " +
                                "reads. [Default: 10]")
    optional_args.add_argument("-ft", "--filter", 
                                dest = "Filter_Option",
                                action = "store_false",
                                help = "Carries out filtering of reads after " +
                                "demultiplexing. [Default: On]. Turn off with just '-ft'")
    optional_args.add_argument("-g", "--gff", 
                                dest = "GFF",
                                action = "store", type = str,
                                help = "The GFF file associated with the reference " +
                                "genome provided with the '-r' option")
    optional_args.add_argument("-k", "--kit", 
                                dest = "Sequencing_Kit",
                                action = "store", type = str,
                                default = "LSK109",
                                help = "The sequencing kit used for the experiment. " +
                                "The SQK prefix is not necessary i.e. for SQK-LSK109, " +
                                "only 'LSK109' should be stated")
    optional_args.add_argument("-kt", "--kraken_Hit_Threshold", 
                                dest = "Kraken_Hit_Threshold",
                                action = "store", type = int, default = 5,
                                help = "A minimum number of groups that must be matched " +
                                "to place a contig into a taxonomic group")
    optional_args.add_argument("-km", "--kraken_mode",
                                dest = "Kraken_Mode",
                                action = "store",
                                choices = ["reads", "assembly", "both"],
                                type = str, default = "assembly",
                                help = "Choose the method to use with kraken taxonomic " +
                                "classification. Classification can be done with just " +
                                "the host-free reads, the metagenome assembly or both. " +
                                "Default is [assembly].")
    optional_args.add_argument("-n", "--name", 
                                dest = "Experiment_Name",
                                action = "store", type = str,
                                help = "Give a name for the experiment. This will be used " +
                                "to name files and folders where appropriate")
    optional_args.add_argument("-o", "--out_dir", 
                                dest = "Output_Folder",
                                action = "store", 
                                default = os.path.join(file_directory, "Pipeline_Output"),
                                type = str,
                                help = "Full path to the output directory. Default will " +
                                "create the output file in the current working directory")
    optional_args.add_argument("-p", "--expansion", 
                                dest = "Expansion_Kit",
                                action = "store", type = str, default = "NBD104",
                                help = "The expansion kit used for barcoding isolates " +
                                "in this experiment. The 'EXP' prefix is not necessary. " +
                                "i.e. EXP-NBD104 should be given as NBD104. " +
                                "[Default: NBD104]")
    optional_args.add_argument("-rd", "--redo_demulp", 
                                dest = "Redo_Demultiplex",
                                action = "store_true", 
                                help = "Carry out demultiplexing alone. To be used to " +
                                "re-generate demultiplexed reads if previously deleted. " +
                                "[Default: Off], turn on with '-rd' ")
    optional_args.add_argument("-t", "--threads", 
                                dest = "Threads",
                                action = "store", type = int,
                                default = "24", 
                                help = "Number of threads to use. [Default: 24]")
    optional_args.add_argument("-w", "--cleanup", 
                                dest = "Cleanup",
                                action = "store_true",
                                help = "To clean up temporary files and other generated " +
                                "files as the pipeline progresses. Large generated files " +
                                "and folders would be deleted e.g. Basecalled outputs, " +
                                "demultiplexed reads and others will be deleted after " +
                                "being used. Log files will detail which files and " +
                                "folders were deleted. [Default: Off], turn on with '-w'")
    optional_args.add_argument("-mxm", "--max_memory", 
                                dest = "Maximum_Memory",
                                action = "store", type = str,
                                default = "150G",
                                help = "Maximum memory to usefor the pipeline. Large " +
                                "input files and large classification databases may " +
                                "require a large amount of memory. Default is based " +
                                "on size of complete Refseq Kraken database which is " +
                                "~132G. [Default: 150G]")
    optional_args.add_argument("-h", "--help", action = "help",
                                default = argparse.SUPPRESS,
                                help = "Displays this help message")
    #########################################################################################
    args = parser.parse_args()
    return args, file_directory


#############################################################################################
# set global variables
args, FILE_DIRECTORY = get_args()
INP_DIR = args.Basecalled_Reads
OUT_DIR = args.Output_Folder
EXP_NAME = args.Experiment_Name
REFERENCE = args.Reference_Genome
BARCODES = args.Barcodes
ISOLATES = args.Isolates
DEMULP_CHOICE = args.Demultiplexer
FLOWCELL = args.Flowcell
FILTER_PASS = args.Filter_Option
THREADS = args.Threads
CLEAN = args.Cleanup
REDEMULP = args.Redo_Demultiplex
REF_GFF = args.GFF
SEQ_TYP = args.Sequence_Type
CREF = args.cDNA_Reference
KRAK = args.Kraken_PATH
KRAKDB = args.Kraken_DB_Path
KRAK_THRESH = args.Kraken_Hit_Threshold
MXMEM = args.Maximum_Memory
MAKCREF = args.make_cDNA_Reference
CREADS = args.reference_cDNA_Reads
CADAP = args.cDNA_Adapters
BRAK = args.Bracken_PATH
BRAKTHRESH = args.Bracken_Hit_Threshold 
NCBI_DB = args.NCBI_DB      
KMODE = args.Kraken_Mode                  
SCPTS = os.path.join(FILE_DIRECTORY, "Scripts") # Scripts folder will be part of the package
# INDEX = os.path.join(FILEDIRETORY, "Index") # Index folder will be part of the package. WILL INCLUDE TruSeq.fa and readme file with links for the SRR sequences and the macaca nemestrina downloads
####################################################################################################################################################
''' Run the script and functions '''
####################################################################################################################################################
# Setting up logging
if __name__ == '__main__':
    logger = logging.getLogger("NanoMetaPipe.py: %s" % time.asctime())
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
        logstream = open(OUT_DIR + '/logINFO.log', 'w')
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
#############################################################################################
''' Check inputs are as expected '''
#############################################################################################
# make lists for DNA and cDNA isolates
print = logger.info
DNA_ISOLATE = []
CDNA_ISOLATE = []
# check the isolate names given fit the needed format
for isolate in ISOLATES:
    iso = isolate.lower()
    # if the isolate has DNA in its name
    if (iso.__contains__("_dna")):
        logger.info('You have provided DNA sequences')
        logger.info(isolate)
        # add the isolate to the DNA list
        DNA_ISOLATE.append(isolate)
    # if the isolate has cDNA in its name
    elif (iso.__contains__("cdna")):
        logger.info('You have provided cDNA sequences')
        logger.info(isolate)
        # add the isolate to the cDNA list
        CDNA_ISOLATE.append(isolate)
    else:
        logger.info('You have not provided the isolates in a satisfactory format')
        logger.info('Do all your isolate names have _dna or _cdna or _dscdna?')
        print('Please look at the help message and try again')
        sys.exit(1)
#############################################################################################
''''Not working'''
""" DNA_ISOLATE = []
CDNA_ISOLATE = []
DNA_ISOLATE.append, CDNA_ISOLATE.append, DNAPres, CDNAPres = isolateList(ISOLATES)
if DNAPres is True:
    print('You have provided DNA sequences')
    print(DNA_ISOLATE)
elif CDNAPres is True:
    print('You have provided cDNA sequences')
    print(CDNA_ISOLATE)
else:
    print('You have not provided the isolates in a satisfactory format')
    print('Do all your isolate names have _dna or _cdna or _dscdna?')
    print('Please look at the help and try again')
    sys.exit(1) """
#############################################################################################
''' check if the a cDNA reference genome is given '''
# set the sequence type based on the user's options
noCref, Cref, CrefMake, noRef, noAdap, Adap, makCref, AllCheck, NoSeqType = seqCheck(SEQ_TYP, 
                        MAKCREF, CREADS, REFERENCE, CREF, CADAP)
if noCref is True:
    logger.info('You have chosen to make a transcriptome assembly but provided no reads.' +
    'Please provide the reads and run again')
    sys.exit(1)
elif Cref is True:
    logger.info('You are giving both DNA and cDNA reads and have selected to ' +
                'generate a transcriptome assembly.')
    logger.info(MAKCREF)
elif CrefMake is True:
    CREF = REFERENCE
elif noRef is True:
    logger.info('Did you provide the two needed references?')
    logger.info('Please remember a DNA reference and a transcriptome reference ' +
                'are required')
    sys.exit(1)
elif noCref is True:
    logger.info('You have chosen to make a transcriptome assembly but provided no reads')
    logger.info('Please provide the reads and run again')
    sys.exit(1)
elif noAdap is True:
    logger.info('You have chosen to make a transcriptome assembly but provided ' +
                'no adapters for them. Please use the -ca option to add the adapter file')
    sys.exit(1)
elif makCref is False:
    logger.info('Did you provide the two needed references?')
    logger.info('Please remember a DNA reference and a transcriptome reference ' +
                    'are required')
    logger.info('Please try again')
    sys.exit(1)
elif NoSeqType is True:
    logger.info('Invalid sequence type entry. Please choose between dna, cdna or both')
    sys.exit(1)
elif AllCheck is True:
    logger.info('All checks are correct. Continuing')
    pass
else: 
    logger.info("Some checks have failed. Please check your options again")
'''set some variables'''
# check the demultiplexer choice
G_KIT, Q_KIT, DBRACK_LENGTH, CBRACK_LENGTH, FILT_QUAL, DNA_FILT_LENGTH, CDNA_FILT_LENGTH  = filterOptions(
    DEMULP_CHOICE, args, FILTER_PASS)
if CLEAN is True:
    logger.info("You have chosen to do clean up. Large files and directories" +
          "will be deleted as the pipelines progress. The most " +
          "consequential of these include the deletion of the " +
          "basecalled reads as these can simply be remade with " +
          "the basecalling script and demultiplexing can be " +
          "with the -rd parameter")
#############################################################################################
''' Making global directories '''
#############################################################################################
# folder to hold demultiplexed reads
dem_dir = os.path.join(OUT_DIR, "Demultiplexed")
makeDirectory(dem_dir)
# folder to hold stats outputs
stats_dir = os.path.join(OUT_DIR, "Stats")
makeDirectory(stats_dir)
# folder to hold temporary alignment files. To be deleted with cleanup option
temp_align_out = os.path.join(OUT_DIR, "Isolate_Reads_Aligned_Vs_Reference")
makeDirectory(temp_align_out)
# folder to hold cleaned reads
aligned_out = os.path.join(OUT_DIR, "Host_Free_Reads")
makeDirectory(aligned_out)
# if filter option is turned off, use raw reads
if FILTER_PASS is False:
    rnm_path = os.path.join(OUT_DIR, "Isolate_Demultiplexed_Reads")
    makeDirectory(rnm_path)
#############################################################################################
''' Call functions '''
#############################################################################################
'''Demultiplexing. Function: demultip is in the Preprocessing.py script'''
    # check if the user just wants to re-try demultiplexing
if REDEMULP is True:
    logger.info('You are just re-generating the demultiplexed reads')
    # run the demultiplexing function from 'Preprocessing.py' script
    demultip(INP_DIR, dem_dir, DEMULP_CHOICE, THREADS, Q_KIT)
    # finish demultiplexing and exit
    logger.info('Demultiplexed reads regenerated')
    sys.exit(1)
else:
    # run demultiplexing function
    demultip(INP_DIR, dem_dir, DEMULP_CHOICE, THREADS, Q_KIT)
    pass
#############################################################################################
''' Zip the demultiplexed reads '''
#############################################################################################
logger.info("Zipping the fastq files")
zipFiles(dem_dir, THREADS)
logger.info("Files zipped")
#############################################################################################
'''QC raw reads'''
#############################################################################################
for barcode in BARCODES:
    # make a txt file with the barcode for filename
        ofile = barcode + ".txt"
        # select the zipped demultiplexed file
        dem_file = dem_dir + "/" + barcode + ".fastq.gz"
        # set the stats path
        stats = os.path.join(stats_dir, "Raw_Demultiplexed_Reads", barcode)
        # run the runQC function that is in the 'Preprocessing.py' script
        run_QC(dem_file, barcode, stats, ofile, THREADS)
#############################################################################################
''' filtering: if filter check is true '''
#############################################################################################
if FILTER_PASS is True:
    # if the user sequence type is just DNA
    if SEQ_TYP == "dna":
        ready_path = dna_filter(DNA_ISOLATE, dem_dir, BARCODES, 
                                OUT_DIR, DNA_FILT_LENGTH, FILT_QUAL, 
                                stats_dir, THREADS, SCPTS)
        logger.info('The raw DNA demultiplexed reads have been successfully filtered ' + 
                    'and saved in ' + ready_path)
        logger.info('Please remember that the files are now renamed')
        logger.info('Stats have been done and saved in ' + stats_dir)
    # if the user sequence type is just cdna 
    elif SEQ_TYP == "cdna":
        ready_path = cdna_filter(CDNA_ISOLATE, dem_dir, BARCODES, OUT_DIR,
                                CDNA_FILT_LENGTH, FILT_QUAL, stats_dir, THREADS, SCPTS)
        logger.info('The raw cDNA demultiplexed reads have been successfully filtered ' + 
                    'and saved in ' + ready_path)
        logger.info('Please remember that the files are now renamed')
        logger.info('Stats have been done and saved in ' + stats_dir)
    # if both sequence types are to be used
    elif SEQ_TYP == "both":
        ''' dna filtering first '''
        ready_path = os.path.join(OUT_DIR, "Filtered_Demultiplexed_Reads")
        ready_path = dna_filter(DNA_ISOLATE, dem_dir, BARCODES, 
                                OUT_DIR, DNA_FILT_LENGTH, FILT_QUAL, 
                                stats_dir, THREADS, SCPTS)
        logger.info('The raw DNA demultiplexed reads have been successfully filtered ' + 
                    'and saved in ' + ready_path)
        # cdna filtering
        ready_path = cdna_filter(CDNA_ISOLATE, dem_dir, BARCODES, OUT_DIR,
                                CDNA_FILT_LENGTH, FILT_QUAL, stats_dir, THREADS, SCPTS)
        logger.info('The raw cDNA demultiplexed reads have been successfully filtered ' + 
                    'and saved in ' + ready_path)
        logger.info('Please remember that the files are now renamed')
        logger.info('Stats have been done and saved in ' + stats_dir)
    '''if reads are not being filtered then move the demultiplexed file'''
else:
    logger.info('You chosen not to filter your reads.')
    logger.info('Your demultiplexed reads will be used to proceed')
    ready_path = rnm_path
    count = 0
    for barcode in BARCODES:
        isola = ISOLATES[count]
        dem_file = dem_dir + "/" + barcode + ".fastq.gz"
        rnm_file = ready_path + "/" + isola + ".fastq.gz"
        # set the moving command
        rnMv_Stff = ' '.join(["mv", dem_file, rnm_file])
        print(rnMv_Stff)
        # run the command
        subprocess.call(rnMv_Stff, shell=True)
        count += 1
#############################################################################################
''' align against reference genome and de-duplicate first check if the sequence type 
    is dna or cdna align function is in DNA_processing.py'''
#############################################################################################
if SEQ_TYP == "dna":
    # use the DNA isolate list
    print('Aligning reads against the reference')
    DNA_align(ready_path, DNA_ISOLATE, stats_dir, temp_align_out, 
                aligned_out, THREADS, REFERENCE)
    logger.info('Alignment complete and saved in ' + aligned_out)
    logger.info('De-duplication complete')
elif SEQ_TYP == "cdna":
    print('Aligning reads against the reference')
    cDNA_align(ready_path, CDNA_ISOLATE, CREADS, OUT_DIR, SCPTS,
                MAKCREF, CADAP, REFERENCE, CREF, REF_GFF,
                THREADS, MXMEM, aligned_out)
    logger.info('Alignment complete and saved in ' + aligned_out)
    logger.info('De-duplication complete')
elif SEQ_TYP == "both":
    # do dna alignment
    print('Aligning reads against the DNA reference')
    DNA_align(ready_path, DNA_ISOLATE, stats_dir, temp_align_out, 
                aligned_out, THREADS, REFERENCE)
    # do cdna alignment
    print('Aligning reads against the transcriptome reference')
    cDNA_align(ready_path, CDNA_ISOLATE, CREADS, OUT_DIR, SCPTS,
                MAKCREF, CADAP, REFERENCE, CREF, REF_GFF,
                THREADS, MXMEM, aligned_out)
    logger.info('DNA and cDNA Alignment complete and saved in ' + aligned_out)
    logger.info('De-duplication complete')
#############################################################################################
'''Assemble the reads'''
#############################################################################################
for isolate in ISOLATES:
    toAlgn = os.path.join(aligned_out, isolate + "VsRef_unmapped_renamed.fastq")
    assemOut = os.path.join(OUT_DIR, "DeNoVo_Assembly", isolate)
    assemSt = os.path.join(stats_dir, "DeNoVo_Assembly", isolate)
    makeDirectory(assemOut)
    makeDirectory(assemSt)
    # assemble
    run_flye(toAlgn, assemOut, THREADS)
    logger.info('Assembly completed with metaFlye')
    # assembly-stats
    assemFile = os.path.join(assemOut, "assembly.fasta")
    assemO = os.path.join(assemSt, "RAW_assembly_stats.txt")
    run_AssemStats(assemFile, assemO)
    # quast
    qAssmO = os.path.join(assemSt, "Quast")
    rawQ = raw_Quast(assemFile, qAssmO, str(THREADS), NCBI_DB)
    logger.info('Quast assessment done and saved in ' + rawQ)
     # copy the assembly file to a more accessible point
    logging.info('The generated assembly will be copied to a new folder')
    assemS = os.path.join(OUT_DIR, "Assemblies", "Raw")
    makeDirectory(assemS)
    coAssemFile = os.path.join(assemS, isolate + ".fasta")
    shutil.copy2(assemFile, coAssemFile)
    logging.info('The generated assembly files have been copied to ' + assemS)
#############################################################################################
'''Kraken classification'''
#############################################################################################
# if the user chose to use cdna reads
if SEQ_TYP == "cdna":
    logger.info('Carrying out Kraken and Bracken classification for cDNA isolates')
    for isolate in CDNA_ISOLATE:
        # run the kraken classification function
        if KMODE == "assembly":
            assemFile = os.path.join(OUT_DIR, "DeNoVo_Assembly", isolate, "assembly.fasta")
            readFile = os.path.join(OUT_DIR, "DeNoVo_Assembly", isolate, "assembly.fasta")## needed but not used
            ckrakOut = krakBrak(KRAK, KRAKDB, BRAK, isolate, assemFile, readFile, OUT_DIR,
                            CBRACK_LENGTH, KRAK_THRESH, BRAKTHRESH, THREADS, KMODE)
        elif KMODE == "reads":
            readFile = os.path.join(aligned_out, isolate + "VsRef_unmapped_renamed.fastq")
            assemFile = os.path.join(aligned_out, isolate + "VsRef_unmapped_renamed.fastq") ## needed but not used
            ckrakOut = krakBrak(KRAK, KRAKDB, BRAK, isolate, assemFile, readFile, OUT_DIR,
                            CBRACK_LENGTH, KRAK_THRESH, BRAKTHRESH, THREADS, KMODE)
        elif KMODE == "both":
            assemFile = os.path.join(OUT_DIR, "DeNoVo_Assembly", isolate, "assembly.fasta")
            readFile = os.path.join(aligned_out, isolate + "VsRef_unmapped_renamed.fastq")
            ckrakOut1, ckrakOut2 = krakBrak(KRAK, KRAKDB, BRAK, isolate, assemFile, readFile, OUT_DIR,
                            CBRACK_LENGTH, KRAK_THRESH, BRAKTHRESH, THREADS, KMODE)
        logger.info('Kraken and Bracken completed for ' + isolate)
        logger.info('Kraken and Bracken outputs saved in ' + ckrakOut)
elif SEQ_TYP == "dna":
    logger.info('Carrying out Kraken and Bracken classification for DNA isolates')
    for isolate in DNA_ISOLATE:
        if KMODE == "assembly":
            assemFile = os.path.join(OUT_DIR, "DeNoVo_Assembly", isolate, "assembly.fasta")
            readFile = os.path.join(OUT_DIR, "DeNoVo_Assembly", isolate, "assembly.fasta")## needed but not used
            krakOut = krakBrak(KRAK, KRAKDB, BRAK, isolate, assemFile, readFile, OUT_DIR,
                            DBRACK_LENGTH, KRAK_THRESH, BRAKTHRESH, THREADS, KMODE)
        elif KMODE == "reads":
            readFile = os.path.join(aligned_out, isolate + "VsRef_unmapped_renamed.fastq")
            assemFile = os.path.join(aligned_out, isolate + "VsRef_unmapped_renamed.fastq") ## needed but not used
            krakOut = krakBrak(KRAK, KRAKDB, BRAK, isolate, assemFile, readFile, OUT_DIR,
                            DBRACK_LENGTH, KRAK_THRESH, BRAKTHRESH, THREADS, KMODE)
        elif KMODE == "both":
            assemFile = os.path.join(OUT_DIR, "DeNoVo_Assembly", isolate, "assembly.fasta")
            readFile = os.path.join(aligned_out, isolate + "VsRef_unmapped_renamed.fastq")
            krakOut1, krakOut2 = krakBrak(KRAK, KRAKDB, BRAK, isolate, assemFile, readFile, OUT_DIR,
                                DBRACK_LENGTH, KRAK_THRESH, BRAKTHRESH, THREADS, KMODE)
        logger.info('Kraken and Bracken completed for ' + isolate)
        logger.info('Kraken and Bracken outputs saved in ' + krakOut)
elif SEQ_TYP == "both":
    logger.info('Carrying out Kraken and Bracken classification for cDNA isolates')
    for isolate in CDNA_ISOLATE:
        # run the kraken classification function
        if KMODE == "assembly":
            assemFile = os.path.join(OUT_DIR, "DeNoVo_Assembly", isolate, "assembly.fasta")
            readFile = os.path.join(OUT_DIR, "DeNoVo_Assembly", isolate, "assembly.fasta")## needed but not used
            ckrakOut = krakBrak(KRAK, KRAKDB, BRAK, isolate, assemFile, readFile, OUT_DIR,
                            CBRACK_LENGTH, KRAK_THRESH, BRAKTHRESH, THREADS, KMODE)
        elif KMODE == "reads":
            readFile = os.path.join(aligned_out, isolate + "VsRef_unmapped_renamed.fastq")
            assemFile = os.path.join(aligned_out, isolate + "VsRef_unmapped_renamed.fastq") ## needed but not used
            ckrakOut = krakBrak(KRAK, KRAKDB, BRAK, isolate, assemFile, readFile, OUT_DIR,
                            CBRACK_LENGTH, KRAK_THRESH, BRAKTHRESH, THREADS, KMODE)
        elif KMODE == "both":
            assemFile = os.path.join(OUT_DIR, "DeNoVo_Assembly", isolate, "assembly.fasta")
            readFile = os.path.join(aligned_out, isolate + "VsRef_unmapped_renamed.fastq")
            ckrakOut1, ckrakOut2 = krakBrak(KRAK, KRAKDB, BRAK, isolate, assemFile, readFile, OUT_DIR,
                                CBRACK_LENGTH, KRAK_THRESH, BRAKTHRESH, THREADS, KMODE)
        logger.info('Kraken and Bracken completed for ' + isolate)
        logger.info('Kraken and Bracken outputs saved in ' + ckrakOut)
    logger.info('Carrying out Kraken and Bracken classification for DNA isolates')
    for isolate in DNA_ISOLATE:
        if KMODE == "assembly":
            assemFile = os.path.join(OUT_DIR, "DeNoVo_Assembly", isolate, "assembly.fasta")
            readFile = os.path.join(OUT_DIR, "DeNoVo_Assembly", isolate, "assembly.fasta")## needed but not used
            krakOut = krakBrak(KRAK, KRAKDB, BRAK, isolate, assemFile, readFile, OUT_DIR,
                            DBRACK_LENGTH, KRAK_THRESH, BRAKTHRESH, THREADS, KMODE)
        elif KMODE == "reads":
            readFile = os.path.join(aligned_out, isolate + "VsRef_unmapped_renamed.fastq")
            assemFile = os.path.join(aligned_out, isolate + "VsRef_unmapped_renamed.fastq") ## needed but not used
            krakOut = krakBrak(KRAK, KRAKDB, BRAK, isolate, assemFile, readFile, OUT_DIR,
                            DBRACK_LENGTH, KRAK_THRESH, BRAKTHRESH, THREADS, KMODE)
        elif KMODE == "both":
            assemFile = os.path.join(OUT_DIR, "DeNoVo_Assembly", isolate, "assembly.fasta")
            readFile = os.path.join(aligned_out, isolate + "VsRef_unmapped_renamed.fastq")
            krakOut1, krakOut2 = krakBrak(KRAK, KRAKDB, BRAK, isolate, assemFile, readFile, OUT_DIR,
                            DBRACK_LENGTH, KRAK_THRESH, BRAKTHRESH, THREADS, KMODE)
        logger.info('Kraken and Bracken completed for ' + isolate)
        logger.info('Kraken and Bracken outputs saved in ' + krakOut)
logger.info("Taxonomic classification completed")
if SEQ_TYP == "cdna":
    if KMODE == "assembly" or KMODE == "reads":
        logger.info('Kraken and Bracken outputs are saved in ' + os.path.dirname(ckrakOut))
        logger.info
    elif KMODE == "both":
        logger.info('Kraken and Bracken outputs are saved in ' + os.path.dirname(ckrakOut1) + ' and ' + os.path.dirname(ckrakOut2))
elif SEQ_TYP == "dna":
    if KMODE == "assembly" or KMODE == "reads":
        logger.info('Kraken and Bracken outputs are saved in ' + os.path.dirname(krakOut))
    elif KMODE == "both":
        logger.info('Kraken and Bracken outputs are saved in ' + os.path.dirname(krakOut1) + ' and ' + os.path.dirname(krakOut2))
elif SEQ_TYP == "both":
    if KMODE == "assembly" or KMODE == "reads":
        logger.info('Kraken and Bracken outputs for the cdna samples are saved in ' + os.path.dirname(ckrakOut))
        logger.info('Kraken and Bracken outputs for the dna samples are saved in ' + os.path.dirname(krakOut))
    elif KMODE == "both":
        logger.info('Kraken and Bracken outputs for the cdna assemblies are saved in ' + os.path.dirname(ckrakOut1))
        logger.info('Kraken and Bracken outputs for the cdna reads are saved in ' + os.path.dirname(ckrakOut2))
        logger.info('Kraken and Bracken outputs for the dna assemblies are saved in ' + os.path.dirname(krakOut1))
        logger.info('Kraken and Bracken outputs for the dna reads are saved in ' + os.path.dirname(krakOut2))
    
