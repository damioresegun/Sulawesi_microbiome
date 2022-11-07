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
from pathlib import Path
from Bio import SeqIO

from Scripts.Tools import makeDirectory

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
    required_args.add_argument("-db", "--database",
                                dest = "Database", type = str,
                                action = "store",
                                required = True,
                                help = "Path to the downloaded SILVA database")
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
    optional_args.add_argument("-pp", "--profile_prefix",
                                dest = "Profile_prefix", type = str,
                                action = "store",
                                default = "MyHMMProfile",
                                help = "Only to be used when creating a profile. Use this option " +
                                "provide a specific prefix name to save the create HMM profile. " +
                                "Default is [MyHMMProfile]")
    optional_args.add_argument("-s", "--isolate",
                                dest = "Isolate", type = str,
                                action = "store",
                                default = "Isolate",
                                help = "Give a name for the isolate. If not chosen, the " +
                                "script will attempt to find the isolate name from either " +
                                "filename or the directory name. Default is [Isolate]")
    optional_args.add_argument("-incE", "--inclusion_evalue",
                                dest = "InclusionExpect",
                                action = "store", type = float,
                                default = 0.01,
                                help = "The inclusion expectation value threshold. E-values " +
                                "above the value given are excluded from the outputs. " +
                                "Default is 0.01 which is an expectation of 1 false positive in " +
                                "100 query searches.")
    optional_args.add_argument("-E", "--e_value",
                                dest = "ExpectationVal",
                                action = "store", type = float,
                                default = 5,
                                help = "Expectation value for false positives. This option "+
                                "works to allow some noise through and for the user to "+
                                "determine the most appropriate hit based on their data needs." +
                                " Default is 5 i.e. No more than 5 false positives may be " +
                                "returned per query sequence")
    optional_args.add_argument("-t", "--threads",
                                dest = "Threads", type = int,
                                action = "store",
                                default = 2,
                                help = "Number of threads to use. Default is 2")
    optional_args.add_argument("-vt", "--vsearch_threshold", 
                                dest = "Vsearch_Threshold",
                                action = "store", type = float,
                                default = 0.65,
                                help = "Enter a number between 0 and 1 to act as a similarity " +
                                "threshold for vsearch global search. Note: the higher the " + 
                                "threshold value, the more stringent and thus the less the " +
                                "number of reads returned in the BLAST output")
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
PFPREF = args.Profile_prefix
INCE = args.InclusionExpect
EVAL = args.ExpectationVal
DBASE = args.Database
VTHRES = args.Vsearch_Threshold
SCPTS = os.path.join(FILE_DIRECTORY, "Scripts")
####################################################################################################################################################
''' Run the script and functions '''
####################################################################################################################################################
def checkInputs(smode,ffile):
    '''
    Function to check the inputs into the script. Checks if the user wants to create a HMM profile 
    or if the user provides a HMM profile. It also checks if the provided profile has already been indexed.
    It also checks for zipped input files and for fastq files and converts these to fasta. 
    '''
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
            chekHmInd1 = PFILE+".h3f"
            chekHmInd2 = PFILE+".h3i"
            if not chekHmInd1 or not chekHmInd2:
              print("HMM profile not empty.")
              print("HMM profile is not indexed. This is good.")
              print("The script will proceed")
              print("However, be aware that if the file used with -p is not a HMM profile, " +
            "the script will fail downstream.")
            else:
              print("The HMM profile has been previously indexed")
              print("The indexes will be deleted")
              rmV = ' '.join(["rm", PFILE+".*"])
              subprocess.call(rmV, shell = True)
    if ffile.endswith('.gz'):
        print("Input reads file is Gzipped. Unzipping now")
        deZip = ' '.join(["pigz -d", ffile, "-p"+str(THREADS)])
        print(deZip)
        subprocess.call(deZip, shell = True)
        ffile = os.path.abspath(ffile).split(".gz")[0]
        print("Input reads file unzipped")
        makQ = ' '.join(["sed -n '1~4s/^@/>/p; 2~4p' ", ffile, ">", 
                    os.path.abspath(ffile).split(".fastq")[0]+".fasta"])
        print(makQ)
        subprocess.call(makQ, shell = True)
        ffile = os.path.abspath(ffile).split(".fastq")[0]+".fasta"
    if ffile.endswith('.fastq'):
        print("Input reads is a FASTQ. Will convert to FASTA")
        makQ = ' '.join(["sed -n '1~4s/^@/>/p; 2~4p' ", ffile, ">", 
                    os.path.abspath(ffile).split(".fastq")[0]+".fasta"])
        print(makQ)
        subprocess.call(makQ, shell = True)
        ffile = os.path.abspath(ffile).split(".fastq")[0]+".fasta"
    print("All checks done. Proceeding...")
    return ffile


def fastaExtract(hitlist, allReads, output):
    '''
    Function to carry out extraction of fasta sequences based on a list of readIDs
    Takes in a txt file containing a list of readIDs, the full fasta file containing all reads
    and requires the path to the output fasta file
    '''
    readsList = open(hitlist, 'r')
    outputfile = open(output, 'w')
    wanted = set()
    with readsList as f:
        for line in f:
            line = line.strip()
            if line != "":
                wanted.add(line)
    fasta_sequences = SeqIO.parse(open(allReads),'fasta')
    with outputfile as i:
        for seq in fasta_sequences:
            if seq.id in wanted:
                SeqIO.write([seq], i, "fasta")
    readsList.close()
    outputfile.close()

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


def alnMuscle(candidates, prefix, outdir, threads):
    '''
    Function to carry out multiple sequence alignments of candidate SSU sequences
    that were downloaded to create a HMM profile
    Takes in the candidates fasta file, a prefix to save the output files with,
    the output directory and the number of threads to run the commands.

    Function uses muscle for multiple sequence alignment and then trimal to trim
    the alignment file 
    '''
    # check if the output directory exists
    makeDirectory(outdir)
    alnOut = os.path.join(outdir, prefix+"_candidates_aligned.aln")
    runMusc = ' '.join(["muscle -align", candidates, "-output", alnOut, "--threads", str(threads)])
    print(runMusc)
    #subprocess.call(runMusc, shell = True)
    print("Muscle multiple sequence alignment complete")
    print("Starting trimal")
    trmOut = os.path.join(outdir, prefix+"_candidates_aligned_clean.fasta")
    runTrml = ' '.join(["trimal -in", alnOut, "-out", trmOut, "-gappyout"])
    print(runTrml)
    #subprocess.call(runTrml, shell = True)
    print("Multiple sequence alignment completed")
    return trmOut


def hmmprf(alignedFile, outDIR, prefix, threads):
    '''
    Function to build a HMM profile using an inputted alignment file.
    The function uses hmmbuild to make the profile based on the input
    alignment file. 
    Returns the path to the generated hmm
    '''
    # set the output file
    outfile = os.path.join(outDIR, prefix+".hmm")
    runBuild = ' '.join(["hmmbuild --dna --cpu", str(threads), "-n", prefix, outfile, alignedFile])
    print(runBuild)
    subprocess.call(runBuild, shell = True)
    print("HMM Profile built.")
    return outfile


def hmmAln(profile, outdir, isolate, fasta):
    '''
    Function to carry out alignment of the isolate reads against the generated or provided HMM profile.
    First the profile is indexed and then nhmmscan is carried out to scan the reads against the profile.
    Outputs multiple files: RawMatches, PerSequenceMatches and Short_PerHit text files. 
    '''
    # index the profile
    print("Indexing the HMM profile")
    indBuild = ' '.join(["hmmpress", profile])
    print(indBuild)
    #subprocess.call(indBuild, shell = True)
    print("HMM Profile constructed and indexed")
    # then nhmmscan
    outIso = os.path.join(outdir, isolate)
    makeDirectory(outIso)
    runHmScan = ' '.join(["nhmmscan --incE", str(INCE), "-E", str(EVAL), "--cpu", str(THREADS),
                            "-o", outIso+"/RawMatches.txt", "--tblout", outIso+"/PerSequenceMatches.txt",
                            "--dfamtblout", outIso+"/Short_PerHit.txt", profile, fasta])
    print(runHmScan)
    #subprocess.call(runHmScan, shell = True)
    outRaw = os.path.join(outIso, "RawMatches.txt")
    outPerSeq = os.path.join(outIso, "PerSequenceMatches.txt")
    outPerHit = os.path.join(outIso, "Short_PerHit.txt")
    print("Hmmer scan complete")
    return outRaw, outPerSeq, outPerHit


def getReads(input, output, isolate, fasta):
    '''
    Function to extract the reads that have been identified to align against the HMM profile
    Takes in the path to the PerSequenceMatches file, the output folder, the isolate name
    and the host-free reads in FASTA format. 
    Function first extracts the read IDs of reads that align to the HMM profile and then calls
    a secondary script (fastaExtractor.py) to then use that list of readIDs to make a new FASTA
    file of just reads that are identified to correspond to the SSU. 
    Outputs the list of read IDs and makes a fasta file containing the reads which align to
    the HMM profile.
    Returns list of readIDs and Reads fasta
    '''
    if Path(input).is_file() and os.path.getsize(input) > 0:
        print("PerSequenceMatches file is not empty, script continues")
    else:
        print("PerSequenceMatches file may not exist or is empty. Please try again")
        sys.exit()
    print("Extracting the read IDs that align to the HMM profile")
    makeDirectory(output)
    extID = os.path.join(output, isolate+"_List_of_Hit_Reads.txt")
    runExt = ' '.join(["cat", input, "| cut -f19 -d ' ' | tail -n +3 | uniq | head -n -2 >", extID])
    print(runExt)
    #subprocess.call(runExt, shell=True) 
    print("Identified read IDs have been extracted. Corresponding reads will now be extracted")
    extReads = os.path.join(output, isolate+"_SSU_reads.fasta")
    print(fastaExtract(extID, fasta, extReads))
    #runReExt = ' '.join(["python3", fastaExtract, extID, fasta, extReads])
    #print(runReExt)
    #subprocess.call(runReExt, shell = True)
    print("SSU reads have been extracted into "+ extReads)
    return extID, extReads


def vsrch(reads, DBASE, outdir, threads):
    '''
    Function to carry out vsearch global search of isolate SSU reads against the SILVA database
    SILVA database must be downloaded and provided by the user. Function first carries out 
    dereplication of the reads and the database before carrying out the search
    Requires: SSU reads, SILVA database, Output directory, threads
    '''
    # check if the reads and the database has been previously dereplicated
    checkReds = os.path.abspath(reads).split(".fasta")[0]+"_derep.fasta"
    checkDB = os.path.abspath(DBASE).split(".fasta")[0]+"_derep.fasta"
    if Path(checkReds).is_file(): 
        if os.path.getsize(checkReds) > 0:
            print("The reads have already been dereplicated and the file is not empty")
            deRepRead = checkReds
        else:
            print("The input reads appear to have a dereplicated file but it is empty")
            print("This input file will be deleted and reconstituted")
            # set the output file
            deRepRead = checkReds
            # dereplicate the reads
            runDeRepReads = ' '.join(["vsearch --derep_fulllength", reads, "--output", deRepRead, "--sizeout"])
            print(runDeRepReads)
            #subprocess.call(runDeRepReads, shell = True)
            print("Reads have been dereplicated")
    else:
        print("The reads file has not been dereplicated. Dereplication will now take place")
        # set the output file
        deRepRead = checkReds
        # dereplicate the reads
        runDeRepReads = ' '.join(["vsearch --derep_fulllength", reads, "--output", deRepRead, "--sizeout"])
        print(runDeRepReads)
        #subprocess.call(runDeRepReads, shell = True)
        print("Reads have been dereplicated")
    if Path(checkDB).is_file(): 
        if os.path.getsize(checkDB) > 0:
            print("The database has already been dereplicated and the file is not empty")
            deRepDB = checkDB
        else:
            print("The database appears to have a dereplicated file but it is empty")
            print("This empty dereplicated file will be deleted and reconstituted")
            # set the output file
            deRepDB = checkDB
            # dereplicate the reads
            runDeRepDB = ' '.join(["vsearch --derep_fulllength", DBASE, "--output", deRepDB, "--sizeout"])
            print(runDeRepDB)
            #subprocess.call(runDeRepDB, shell = True)
            print("Database has been dereplicated")
    else:
        print("The database has not been dereplicated. Dereplication will now take place")
        # set the output file
        deRepDB = checkDB
        # dereplicate the reads
        runDeRepDB = ' '.join(["vsearch --derep_fulllength", DBASE, "--output", deRepDB, "--sizeout"])
        print(runDeRepDB)
        #subprocess.call(runDeRepDB, shell = True)
        print("Database has been dereplicated")
    print("Proceeding to carry out global search")
    try:
        runGlobSrch = ' '.join(["vsearch --usearch_global", deRepRead, "--uc", outdir+"/VSEARCH_clusters_UC.tab",
                                "--blast6out", outdir+"/VSEARCH_BLAST_alignedReads.txt", "--db", deRepDB, "--id", str(VTHRES),
                                "--threads", str(threads)])
        print(runGlobSrch)
        #subprocess.call(runGlobSrch, shell = True)
        print("Global search completed and results have been saved in the output directory")
        blastOut = outdir+"/VSEARCH_clusters_UC.tab"
        clustOut = outdir+"/VSEARCH_BLAST_alignedReads.txt"
        return blastOut, clustOut
    except (FileExistsError, FileExistsError) as err:
        print(type(err), err)


if __name__ == '__main__':
    FFILE = checkInputs(SMODE,FFILE)
    # if the user chooses to create profile:
    if SMODE == "create":
        print("You have chosen generate a profile. Proceeding")
        # download the candidates
        downSC = os.path.join(SCPTS, "Download_NCBI_Accession_to_Fasta.py")
        candidateFile = downloadCandidates(CFILE, downSC, OUTDIR)
        # carry out multiple sequence alignment
        trmmed = alnMuscle(candidateFile, PFPREF, OUTDIR, THREADS)
        # build the profile
        hmmprofile = hmmprf(trmmed, OUTDIR, PFPREF, THREADS)
        print("HMM Profile has been constructed")
        PFILE = hmmprofile
    else: 
        print("You have provided a profile. Proceeding...")
    print(FFILE)
    # carry out the alignment of the profile against the input fasta file
    outRaw, outPerSeq, outPerHit = hmmAln(PFILE, OUTDIR, ISOLATE, FFILE)
    # extract reads and determine the species associated
    extIDs, extReads = getReads(outPerSeq, OUTDIR, ISOLATE, FFILE)
    # carry out vsearch of the identified reads against the SILVA database
    blastTab, clustTab = vsrch(extReads, DBASE, OUTDIR, THREADS)
    print("Reads have been put through vsearch global search and results saved in " + OUTDIR)
    print("Have a look at " + blastTab)