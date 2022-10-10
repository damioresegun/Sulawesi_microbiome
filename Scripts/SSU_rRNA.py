#!/bin/env python3
'''
Script to carry out a hidden markov model search for the Plasmodium SSU rRNA.
The script will either take in or generate a HMM profile to map input reads against
Requirements:
- list of candidate accession IDs OR a previously made HMM profile
- hmmer v3+
- pigz
- Download_NCBI_Accession_to_Fasta.py script
'''
import os
import sys
import subprocess
import argparse

#######################################################################################################
def get_args():
    '''
    Take in the input arguments and give a help message when prompted
    '''
    parser = argparse.ArgumentParser(description = "Short script to carry out a hidden markov " +
                    "model search and mapping for Plasmodium SSU within a metagenomic community " +
                    "The user will need to provide as input, a list of candidates to build "+
                    "the model with OR a previously generated model. Input reads will also need " +
                    "to be given.",
                    add_help = False)
    file_directory = os.path.realpath(__file__).split("SSU_rRNA")[0]
    if not os.path.isfile(os.path.join(file_directory, "SSU_rRNA.py")):
        file_directory = os.path.join(file_directory, "SSU_rRNA")
    if not os.path.isfile(os.path.join(file_directory, "SSU_rRNA.py")):
        print("Can't locate the correct path to the script")
    #######################################################################################################
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument("-f", "--fasta",
                                dest = "FASTAfile", type = str,
                                action = "store",
                                required = True,
                                help = "Full path to the reads to map against the profile")
    required_args.add_argument("-m", "--mode",
                                dest = "Mode", type = str,
                                choices = ["create", "run"],
                                required = True, 
                                help = "Choose which mode to run the script. In 'create' " +
                                "the script will require a candidates file (-c) to make " +
                                "a HMM profile. In 'run' mode, the script will require a " +
                                "previouly made HMM profile (-p) to run the rest of the script")
    required_args.add_argument("-o", "--output",
                                dest = "OutputDIR", type = str,
                                action = "store",
                                required = True,
                                help = "Full path to the directory to generate outputs into")
    #######################################################################################################
    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument("-c", "--candidates",
                                dest = "Candidates", type = str,
                                action = "store",
                                help = "Full path to a txt file of a list accession IDs for " +
                                        "candidate sequences to generate the HMM profile. Only " +
                                        "use this if you want to generate a HMM profile. Use " +
                                        "'-p' to enter previously generated HMM profiles.")
    optional_args.add_argument("-p", "--profile", 
                                dest = "Profile", type = str,
                                action = "store",
                                help = "Full path to a previously made HMM profile. This script " +
                                "was designed for profile made from DNA sequences! Not proteins!!!")
    optional_args.add_argument("-s", "--isolate",
                                dest = "Isolate", type = str,
                                action = "store",
                                default = "Isolate",
                                help = "Give a name for the isolate. If not chosen, the " +
                                "script will attempt to find the isolate name from either " +
                                "filename or the directory name. Default is [Isolate]")
    optional_args.add_argument("-t", "--threads",
                                dest = "Threads", type = int,
                                action = "store",
                                default = 2,
                                help = "Number of threads to use. Default is 2")
    optional_args.add_argument("-h", "--help",
                                action = "help",
                                default = argparse.SUPPRESS,
                                help = "Displays this help message")
    #######################################################################################################
    args = parser.parse_args()
    return args, file_directory
#######################################################################################################
# set global variables
args, FILE_DIRECTORY = get_args()
FFILE = args.FASTAfile
OUTDIR = args.OutputDIR
SMODE = args.Mode
CFILE = args.Candidates
PFILE = args.Profile
ISOLATE = args.Isolate
THREADS = args.Threads
####################################################################################################################################################
''' Run the script and functions '''
####################################################################################################################################################
def checkInputs(smode,ffile):
    if smode == "create":
        print("Script will run in create mode. Checking for candidate file")
        if not CFILE:
            print("Could not find the candidates file. Please use the -c option")
            print("The candidates files needs to be a list of accession IDs to download from NCBI")
        else:
            print("Candidate file found. Proceeding with script")
    if smode == "run":
        print("Script will in 'run' mode. Checking for profile file")
        if not PFILE:
            print("HMM Profile file could not be found. Please use the -p option")
        else:
            print("HMM profile not empty. The script will proceed")
            print("However, be aware that if the file used with -p is not a HMM profile, " +
            "the script will fail downstream.")
    if ffile.endswith('.gz'):
        print("Input reads file is Gzipped. Unzipping now")
        deZip = ' '.join(["pigz -d", ffile, "-p"+str(THREADS)])
        print(deZip)
        subprocess.call(deZip, shell = True)
        ffile = os.path.abspath(ffile).split(".gz")[0]
        print("Input reads file unzipped")
    if ffile.endswith('.fastq'):
        print("Input reads is a FASTQ. Will convert to FASTA")
        FAFILE = os.path.abspath(ffile).split(".fastq")[0]
        makQ = ' '.join(["sed -n '1~4s/^@/>/p; 2~4p' ", ffile, ">", 
                    os.path.abspath(ffile).split(".fastq")[0]+".fasta"])
        print(makQ)
        subprocess.call(makQ, shell = True)
        ffile = os.path.abspath(ffile).split(".fastq")[0]+".fasta"

checkInputs(SMODE,FFILE)
    


#check if the input reads end with fastq. if so then convert to fasta and then use