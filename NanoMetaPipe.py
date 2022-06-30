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
from Scripts.Preprocessing import demultip, filt_qc run_QC
from Scripts.DNA_Processing import align
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
    required_args.add_argument("-b", "--basecalled",
                                dest = "Basecalled reads", 
                                action = "store",
                                type = str,
                                help = "Folder of basecalled reads in FASTQ format. " +
                                "Reads must be basecalled previously!. " +
                                "Example: path/to/exp_folder/pass. This option does not " +
                                "work with raw reads. ", required = True)
    required_
   