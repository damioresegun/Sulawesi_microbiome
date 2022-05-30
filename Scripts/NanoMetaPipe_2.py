#!/usr/bin/env python3
#
# Name of pipeline: NanoMetaPipe
#
# Pipeline to take basecalled metagenomic nanopore reads through
# preprocessing to final assembled bins for analysis.
#
# MacKenzie Institute for Early Diagnosis (2022)
# Author: Damilola Oresegun, Peter Thorpe

import sys
import os
import configparser
import subprocess
import argparse
import mappy as mp
import logging.handlers
import time
from Bio import SeqIO
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
                               default=os.path.join(file_directory,
                                                    "JapaneseMacaqueData",
                                                    "Basecalled",
                                                    "pass"),
                               type=str,
                               help="Folder of basecalled reads " +
                               "in FASTQ format. Example: path/to/" +
                               "exp_folder/pass. " +
                               "Does not work with raw reads [-r]. " +
                               "Default is for tests only")
    required_args.add_argument("-o", "--out_dir",
                               dest="Output_folder", action="store",
                               default=os.path.join(file_directory,
                                                    "JapaneseMacaqueData"),
                               type=str,
                               help="Path to the output directory. " +
                               "Default will create the output file " +
                               "in your current working directory.")
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
    optional_args.add_argument("-n", "--name",
                               dest='Experiment_Name',
                               action="store",
                               type=str,
                               help="Give the name of " +
                               "the experiment.")
    optional_args.add_argument("-c", "--barcodes",
                               dest='Barcodes',
                               action="store",
                               type=list,
                               default=['barcode01', 'barcode02'],
                               help="List of barcodes used")
    optional_args.add_argument("-e", "--isolates",
                               dest='Isolates',
                               action="store",
                               type=list,
                               default=['isolate1_dna', 'isolate2_cDNA'],
                               help="List of isolates used corresponding " +
                               "to the order of barcodes given " +
                               "i.e. barcode01=isolate1.")
    optional_args.add_argument("-k", "--kit",
                               dest='Sequencing_kit',
                               action="store",
                               type=str,
                               default='LSK109',
                               help="The sequencing kit used, without SQK" +
                               ". The default is LSK109")
    optional_args.add_argument("-p", "--expansion",
                               dest='Expansion_kit',
                               action="store",
                               type=str,
                               default="NBD104",
                               help="The expansion kit used for " +
                               "barcoding the isolates sequenced " +
                               "without 'EXP'. Default is NBD104")
    optional_args.add_argument("-flw", "--flowcell",
                               dest='Flowcell',
                               action="store",
                               type=str,
                               default="FLO-MIN106",
                               help="The flowcell used for this " +
                               "experiment. Default is FLO-MIN106")
    optional_args.add_argument("-ft", "--filter",
                               dest='Filter_option',
                               action="store_true",
                               default=True,
                               help="Option to carry " +
                               "out filtering. Default is True. " +
                               "Turn off with --ft False")
    optional_args.add_argument("-fl", "--filter_length",
                               dest='Filter_length',
                               action="store",
                               type=int,
                               default="50",
                               help="Filter length threshold. " +
                               "Default is 50")
    optional_args.add_argument("-fq", "--filter_quality",
                               dest='Filter_Quality',
                               action="store",
                               type=int,
                               default="10",
                               help="Phred quality threshold to filter " +
                               "demultiplexed reads. Default is 10")
    optional_args.add_argument("-t", "--threads",
                               dest='Threads',
                               action="store",
                               type=int,
                               default="24",
                               help="Number of threads. Default is 24")
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
WK_DIR = os.getcwd()
INP_DIR = args.Basecalled_reads
OUT_DIR = args.Output_folder
EXP_NAME = args.Experiment_Name
BARCODES = args.Barcodes
ISOLATES = args.Isolates
DEMULP_CHOICE = args.Demultiplexer
FLOWCELL = args.Flowcell
FILTER_PASS = args.Filter_option
THREADS = args.Threads
# set conditional variables
if EXP_NAME is None:
    pass
    PREFIX = os.path.split(OUT_DIR)[-1].split()[0]
else:
    OUT_DIR = os.path.join(OUT_DIR, EXP_NAME)
    PREFIX = EXP_NAME
#
if DEMULP_CHOICE == "Qcat":
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
# log_dir = os.path.join(OUT_DIR, "LogFiles")
# os.mkdir(log_dir, mode = 0o666)
##########################################################################
# Demultiplexing
##########################################################################


def demultip(dem_dir):
    try:
        if DEMULP_CHOICE == "qcat":
            dem_INP = INP_DIR + "/*"
            qdem = ("cat", dem_INP, "| qcat -b", dem_dir, "--detect-middle -t",
                    str(THREADS), "--trim -k", Q_KIT)
            runQdem = ' '.join(qdem)
            print(runQdem)
            # subprocess.call(runQdem, shell=True)
            print('Demultiplexing complete')
        elif DEMULP_CHOICE == "guppy":
            gdem = ("guppy_barcoder -i", INP_DIR, "-s", dem_dir, "--barcode_kits", Q_KIT,
                    "-r -q 0 -t", str(THREADS), "--compress_fastq -x auto",
                    "--detect_mid_strand_barcodes --trim_barcodes --trim_adapters")
            runGdem = ' '.join(gdem)
            print(runGdem)
            # subprocess.call(runGdem, shell=True)
            print('Demultiplexing complete')
    except OSError as error:
        print(error)
##########################################################################
# Filtering
##########################################################################


def filt_qc(dem_dir, barcode, isolate, stats):
    try:
        temp = barcode + "_" + isolate
        stats_dir = os.path.join(stats, "Filtered_Demultiplexed_Reads", temp)
        if os.path.exists(stats_dir):
            pass
        else:
            os.mkdir(stats_dir)
        print('Starting NanoFilt')
        temp = barcode + ".fastq.gz"
        file_in = os.path.join(dem_dir, temp)
        filt_out = os.path.join(OUT_DIR, "Filtered_RawReads")
        if os.path.exists(filt_out):
            pass
        else:
            os.mkdir(filt_out)
        temp = isolate + ".fastq.gz"
        filt_file_out = os.path.join(filt_out, temp)
        filtSt = ("gunzip -c", file_in, "|NanoFilt -l", FILT_LENGTH, "-q", FILT_QUAL,
                  "| gzip >", filt_file_out)
        runFiltSt = ' '.join(filtSt)
        print(runFiltSt)
        return filt_out, stats_dir
##########################################################################
# Pre-Processing QC
##########################################################################


##########################################################################
# Run the script and functions
##########################################################################
# Setting up logging
if __name__ == '__main__':
    logger = logging.getLogger('NanoMetaPipe.py: %s' % time.asctime())
    # logger.setLevel(logging.DEBUG)
    logger.setLevel(logging.INFO)
    d_err_handler = logging.StreamHandler(sys.stdout)
    d_err_handler.setLevel(logging.INFO)
    d_err_formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
    d_err_handler.setFormatter(d_err_formatter)
    logger.addHandler(d_err_handler)
    e_err_handler = logging.StreamHandler(sys.stderr)
    e_err_handler.setLevel(logging.ERROR)
    e_err_formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
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
    except:
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
    demultip(dem_dir)
    print(THREADS)
    # filtering
    # if filter check is true
    count = 0
    for i in BARCODES:
        isola = ISOLATES[count]
        filtered_path, filtered_statsfilt_qc(dem_dir, i, isola, stats)
        print('The raw demultiplexed reads have been successfully filtered and saved in ' + filtered_path)
        print('The QC stats for the filtered reads are saved in ' + os.path.split(filtered_stats)[-1].split()[0])
        count += 1
        print('Please remember that the files are now renamed')
        print(i + ' is now ' + isola)
        print('')

# print(args)
# print(FILE_DIRECTORY)
# print(args.Isolates[1])
# print(OUT_DIR)
# print(Q_KIT)
# print(G_KIT)
# print(FILTER_PASS)
# print(dem_dir)
