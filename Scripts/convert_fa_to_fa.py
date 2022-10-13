#!/usr/bin/env python
# coding: utf-8
# Title: script to reformat fastq to fasta format
# Why: after pear/ flash have assembled the reads.
# Need to convert from fastq to fa
# author: Peter Thorpe (PT) and Leighton Pritchard
# November 2016. The James Hutton Insitute, Dundee, UK.
# Shared by PT to DRO July 2019
# Adapted: Damilola Oresegun (DRO)
# Why: adapted this script to reformat fast* files and convert fastq to fasta
# Adapted: July 2019. School of Medicine, The University of St Andrews.

import gzip
from Bio import SeqIO

#os imports
import os
from sys import stdin,argv
import sys
from optparse import OptionParser

def convert_file(in_file, out_file, mode):
    if mode == "fq":
        if in_file.endswith('.gz'):
            print("Yes")
            file = gzip.open(in_file, "rt")
            SeqIO.convert(file, "fastq", out_file, "fastq")
        else:
            # Convert fastq to reformatted fastq
            SeqIO.convert(in_file, "fastq", out_file, "fastq")
    elif mode == "fa":
        # Uncomment the one below to convert fastq to fasta
        SeqIO.convert(in_file, "fasta", out_file, "fasta")
    elif mode == "fqfa":
        SeqIO.convert(in_file, "fastq", out_file, "fasta")
    
if "-v" in sys.argv or "--version" in sys.argv:
    print ("v0.0.2")
    sys.exit(0)


usage = """Use as follows:

$ convert_fq_to_fa -i in.fastq -o outfile.fasta -m mode

requires biopython
"""

parser = OptionParser(usage=usage)

parser.add_option("-i", dest="in_file", default=None,
                  help="fastq file")
parser.add_option("-o", "--output", dest="out_file", default=None,
                  help="Output filename, fastq",
                  metavar="FILE")
parser.add_option("-m", "--mode", dest="mode", default=None,
                 help="Which mode to run in. 'fq' to format fastq to fastq. "+
                 "'fa' to formart fasta to fasta. 'fqfa' to format fastq to fasta")


(options, args) = parser.parse_args()

in_file = options.in_file
out_file = options.out_file
mode = options.mode

(options, args) = parser.parse_args()
# run the program
convert_file(in_file, out_file, mode)
