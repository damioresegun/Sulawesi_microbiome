#!/usr/bin/env python3
#
# Name of pipeline: NanoMetaPipe
#
# Pipeline to take raw metagenomic nanopore reads through 
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
#from pycits.tools import convert_fq_to_fa
#from pycits.metapy_tools import make_folder
########################################################################
def get_args():
    parser = argparse.ArgumentParser(description="Pipeline for metagenomics " +
                                     "data pre-processing and analysis for " +
                                     "nanopore sequenced data", add_help=False)
    file_directory = os.path.realpath(__file__).split("NanoMetaPipe")[0]
    if not os.path.isfile(os.path.join(file_directory, "NanoMetaPipe.py")):
        file_directory = os.path.join(file_directory, "NanoMetaPipe")
    if not os.path.isfile(os.path.join(file_directory, "NanoMetaPipe.py")):
        print("Can't locate the correct path to the pipeline script
    args = parser.parse_args()
    return args, file_directory

args, FILE_DIRECTORY = get_args()
print(args)
print(FILE_DIRECTORY)