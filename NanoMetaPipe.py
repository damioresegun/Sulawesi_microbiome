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
from Scripts.AssemblyQC import run_AssemStats, raw_Quast
from Scripts.Racon_Medaka import runRacon, runMedaka
from Scripts.Preprocessing import demultip, filt_qc, run_QC
from Scripts.DNA_processing import align
from Scripts.Deduplication import filter_fastq_file
#############################################################################################
def get_args():
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
        print("Can't locate the correct path to the NewNanoMetaPipe pipeline script")
    #########################################################################################
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument("-b", "--basecalled", dest = "Basecalled reads", 
                                action = "store", type = str,
                                help = "Folder of basecalled reads in FASTQ format. " +
                                "Reads must be basecalled previously!. " +
                                "Example: path/to/exp_folder/pass. This option does not " +
                                "work with raw reads. ", required = True)
    required_args.add_argument("-c", "--barcodes", dest = "Barcodes",
                                action = "store", nargs = "+",
                                help = "List of barcodes used e.g. -c barcode01 barcode02",
                                required = True)
    required_args.add_argument("-e", "--isolates", dest = "Isolates",
                                action = "store", nargs = "+",
                                help = "List of isolates used. Must give sequence type " +
                                "separated by underscore e.g. -e isolate1_dna isolate2_cdna. " +
                                "Must correspond to the order of barcodes given i.e. " +
                                "barcode01=isolate1_dna.", required = True)
    required_args.add_argument("-kr", "kraken", dest = "Kraken_PATH",
                                ation = "store", type = str, 
                                help = "Full path to the install kraken package. " +
                                "If kraken is in the $PATH, simply write kraken2",
                                required = True)
    required_args.add_argument("-kb", "--kraken_DB", dest = "Kraken_DB_path",
                                action = "store", type = str,
                                help = "Full path to the installed kraken database",
                                require = True)
    required_args.add_argument("-r", "--reference", dest = "Reference_Genome",
                                action = "store", type = str,
                                help = "Full path to the reference genome to align against.",
                                required = True)
    required_args.add_argument("-s", "--sequence_type", dest = "Sequence_Type",
                                action = "store", choices = ["dna", "cdna", "both"],
                                type = str, 
                                help = "The type of sequence you want to be " +
                                "analysed. DNA and cDNA are supported individually and " +
                                "together in a multiplexed dataset. Ensure to identify " +
                                "the sequence types in the isolate names given i.e. " +
                                "Sequence files should be given as isolate1_dna or  "+
                                "isolate2_cdna or isolate2_dscdna.", required = True)
    required_args.add_argument("-br", "--bracken", dest = "Bracken_PATH",
                                action = "store", type = str, 
                                help = "Full path to the installed bracken. " +
                                "If in the $PATH, simply write bracken", required = True)
    #########################################################################################
    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument("-bt", "--bracken-hit-threshold", 
                                dest = "Bracken_Hit_Threshold", action = "store",
                                type = str, default = 20,
                                help = "A minimum number of kmers that must be matched to " +
                                "place a contig into a taxonomic group by bracken " +
                                "re-estimation. Default is [20]")
    optional_args.add_argument("-cr", "--cdna_ref", dest = "cDNA_Reference",
                                action = "store", type = str, 
                                help = "If using both DNA and CDNA sequence types, please " +
                                "provide the path to the transcriptome assembly here. " +
                                "If only using cDNA sequence files, this option is not " +
                                "necessary!")
    optional_args.add_argument("-mcr", "--make_cdna_ref", dest = "make_cDNA_Reference",
                                action = "store_true", 
                                help = "Generate a genome-guided transcriptome to remove " +
                                "host sequences from input cDNA sequences. Only works with " +
                                "with the detection of indicated cDNA input isolates. " +
                                "This option does not work with the '-cr' option. " +
                                "[Default:Off]. Activate this option with just '-mcr' or " +
                                "'--make_cdna_ref'. To be used in conjunction with the " +
                                "'-cd' option")
    optional_args.add_argument("-cd", "--ref_cdna_reads", dest = "reference_cDNA_Reads",
                                action = "store", nargs = "+", 
                                help = "Path to cDNA reads to generate a host " +
                                "transcriptome. Must be used with '-mcr' option")
                        