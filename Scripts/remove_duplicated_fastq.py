#!/usr/bin/env python3
################################################################################
# Title: Remove duplicated reads names from a fastq file
#
# Author: Dami Oresegun -- adapted from Peter Thorpe's
# script
################################################################################
#
"""This script uses Biopython to iterate through a fastq
file and remove duplicate names before carrying out
de novo genome assembly on Flye"""
# import modules
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import os
from sys import stdin,argv
import sys
import argparse
from optparse import OptionParser

# set options
def get_args():
    parser = argparse.ArgumentParser(description="Quick script to " +
                                    "rename duplicated reads in a " +
                                    "FASTQ file.")
    ############################################################################
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument("-i", "--input", dest="Input",default=None,
                               type=str, action="store",
                               help="Input FASTQ file",required=True)
    required_args.add_argument("-o", "--output", dest="Output", default=None,
                               type=str, action="store",
                               help="Output filename.fastq. Can be full path",
                               required=True)
    args = parser.parse_args()
    return args

def filter_fastq_file(in_fastq, out_fastq):
    # parse the fastq file, collect names in a set
    # if the name is not in the set, write out the entry
    inp_file = open(in_fastq)
    out_file = open(out_fastq, "w")
    count = 0
    name_set = set([])
    # use enumerate to count iterations
    for i, (title, seq, qual) in enumerate(FastqGeneralIterator(inp_file)):
        # print the title of the read
        if title.rstrip("\n") not in name_set:
            # add the title to the name set
            name_set.add(title.rstrip("\n"))
            # write to file
            out_file.write("@%s\n%s\n+\n%s\n" % (title, seq.upper(), qual))
        else:
            # change the title name by adding counter
            title = title.rstrip("\n") + "_" + str(count)
            out_file.write("@%s\n%s\n+\n%s\n" % (title, seq.upper(), qual))
            #print(title + " is a duplicate. It has been changed")
            count +=1
    print("The total number of duplicated reads is: " + str(count))
    out_file.close()
    inp_file.close()
##############################################################################################
# run the function
args = get_args()
in_file = args.Input
out_file = args.Output

if __name__ == '__main__':
    filter_fastq_file(in_file,out_file)
    print("Fastq renaming complete")

            
