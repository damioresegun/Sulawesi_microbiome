#!/bin/env python3
'''
Script to carry out a hidden markov model search for the Plasmodium SSU rRNA.
The script will either take in or generate a HMM profile to map input reads against
Requirements:
- list of candidate accession IDs OR a previously made HMM profile
- hmmer v3+
- pigz
- muscle
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
SCPTS = os.path.join(FILE_DIRECTORY, "Scripts")
####################################################################################################################################################
''' Run the script and functions '''
####################################################################################################################################################
def checkInputs(smode,ffile):
    if smode == "create":
        print("Script will run in create mode. Checking for candidate file")
        if not CFILE:
            print("Could not find the candidates file. Please use the -c option")
            print("The candidates files needs to be a list of accession IDs to download from NCBI")
            sys.exit()
        else:
            print("Candidate file found. Proceeding with script")
    if smode == "run":
        print("Script will in 'run' mode. Checking for profile file")
        if not PFILE:
            print("HMM Profile file could not be found. Please use the -p option")
            sys.exit()
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
        makQ = ' '.join(["sed -n '1~4s/^@/>/p; 2~4p' ", ffile, ">", 
                    os.path.abspath(ffile).split(".fastq")[0]+".fasta"])
        print(makQ)
        subprocess.call(makQ, shell = True)
        ffile = os.path.abspath(ffile).split(".fastq")[0]+".fasta"
    print("All checks done. Proceeding...")


def downloadCandidates(can_list, downloadScript, outdir):
    '''
    Function to take in the candidate file list of accessions and then carry out downloading the sequences
    Sequences will be downloaded into a single fasta file to be used for multiple sequence alignment
    Requires:
    - candidate list file
    - path to the download script
    - path to output directory to place fasta
    - returns path to the generated fasta file
    '''
    try:
        # make the output file
        outfile = os.path.join(outdir + "DownloadedCandidates.fasta")
        print("Downloading candidate accession list")
        runDown = ' '.join(["python", downloadScript, "-t", can_list, "-o", outfile])
        print(runDown)
        subprocess.call(runDown, shell = True)
        print("Download complete. Checking the output file")
    except (FileExistsError, FileExistsError) as err:
        print(type(err), err)
    try:
        # checking the file size
        if (os.path.getsize(outfile)) == 0:
            print("Fasta file has been downloaded but the file is empty")
            print("Please check your inputs once again")
            sys.exit()
        elif (os.path.getsize(outfile)) > 0:
            print("Fasta file successfully downloaded and populated")
    except (FileExistsError, FileExistsError):
        print(type(err), err)    
    return outfile



checkInputs(SMODE,FFILE)
downSC = os.path.join(SCPTS, "Download_NCBI_Accession_to_Fasta.py")
downloadCandidates(CFILE, downSC, OUTDIR)
