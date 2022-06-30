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
from webbrowser import get
from PipelineDevScripts.NanoMetaPipe_assemblyApproach import BARCODES, DEMULP_CHOICE, EXP_NAME, FILE_DIRECTORY, FILTER_PASS, FLOWCELL, INP_DIR, ISOLATES, MAKCREF, MXMEM, OUT_DIR, REF_GFF, REFERENCE
from Scripts.AssemblyQC import run_AssemStats, raw_Quast
from Scripts.Racon_Medaka import runRacon, runMedaka
from Scripts.Preprocessing import demultip, filt_qc, run_QC
from Scripts.DNA_processing import align
from Scripts.Deduplication import filter_fastq_file
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
                                action = "store", type = str, default = "NDB104",
                                help = "The expansion kit used for barcoding isolates " +
                                "in this experiment. The 'EXP' prefix is not necessary. " +
                                "i.e. EXP-NDB104 should be given as NBD104. " +
                                "[Default: NDB104]")
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
SCPTS = os.path.join(FILE_DIRECTORY, "Scripts") # Scripts folder will be part of the package
# INDEX = os.path.join(FILEDIRETORY, "Index") # Index folder will be part of the package. WILL INCLUDE TruSeq.fa and readme file with links for the SRR sequences and the macaca nemestrina downloads
#############################################################################################

####################################################################################################################################################
# Run the script and functions
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

    logger.info("The arguments are: %s", args)