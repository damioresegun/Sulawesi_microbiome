#!/usr/bin/env python3
'''Script to hold various small tools and recurring commands'''

import subprocess
import os
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def makeDirectory(directory):
    '''Function checks if a directory exists. If so, nothing is done,
    if the folder does not exist, the folder and any parent folder,
    is also made
    Input: Path to directory to be checked/made
    Output: None
    Usage: makeDirectory(path/to/directory)
    '''
    if os.path.exists(directory):
        pass
    else:
        os.makedirs(directory)


def zipFiles(directory, threads):
    '''Function to zip files using multi-threads. Function is set up
    to take in only files, not folder. Input is the path to the 
    directory holding the files to be zipped. For each file in that
    folder, zipping will be done using pigz
    Input: path to the directory holding the files to be zipped
    Output: None
    Usage: zipFiles(path/to/directory)
    '''
    for file in os.listdir(directory):
        zipDir = os.path.join(directory, file)
        runZip = ' '.join(["pigz --best", zipDir, "-p", str(threads)])
        print(runZip)
        subprocess.call(runZip, shell=True)
        print('FastQ files have been zipped')


def filter_fastq_file(in_fastq, out_fastq):
    ''' The function uses Biopython to iterate through a fastq
        file and remove duplicate names before carrying out
        de novo genome assembly on Flye. Function is adapted from
        Peter Thorpe's script. The function parses the fastq file, 
        collect names in a set, if the name is not in the set, 
        write out the entry 
        Input: FastQ file to filter
        Output: Filtered FastQ file
    '''
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


def run_flye(in_file, out_dir, threads):
    '''
    Function to assemble reads using metaFlye.
    Reads which have been extracted from alignments against the reference
    genome will go through de novo genome assembly using the '--meta' flag
    of flye which allows for higher error reads and also allows for
    variance in coverage. Duplicated reads will need to be checked and
    renamed. This is done using a custom python script that checks the name
    of a read and renames it if it is duplicated. Then metaFlye takes place

    Input: Input reads file, path to output directory, number of threads
    Output: Output directory holding isolate assembly folders
    '''
    try:
        # set the command
        runFly = ' '.join(["flye --nano-raw", in_file, "--out-dir", out_dir, "--threads", str(threads), "--meta"])
        print(runFly)
        # run the command
        subprocess.call(runFly, shell=True)
        # catch errors
    except FileNotFoundError:
        print('There was an issue with the file. File not found')
    except SyntaxError:
        print('There was a syntax issue. Double check your syntax')


def run_AssemStats(in_file,out_file):
    '''
    Function to carry out assembly-stats of a fasta/fastq file. 
    File must be uncompressed. 
    Input: Path to fasta/fastq file
    Output: Test file containing the descriptive stats calculated
    '''
    try:
        asSt = ("assembly-stats", in_file, ">", out_file)
        runAsSt = ' '.join(asSt)
        print(runAsSt)
        subprocess.call(runAsSt, shell=True)
    except FileNotFoundError:
        print('File not found error. Check your input file again')


def raw_Quast(assembly,output_dir,threads):
    # make output folder
    try:
        out_dir = os.path.join(output_dir, "Raw")
        if os.path.exists(out_dir):
            print('Output folder for this already exists')
            print('Checking it contains any quast outputs')
            check_out = os.path.join(out_dir, "quast.log")
            if os.path.exists(check_out):
                print('Quast files exist. Nothing will done')
                pass
            else:
                print('No Quast files found in the existing folder. Proceeding with Quast')
                runRawQ = ' '.join(["quast -t", threads, assembly, "-o", out_dir, "-f --circos --meta"])
                print(runRawQ)
                subprocess.call(runRawQ, shell=True)
        else:
            os.makedirs(out_dir)
            runRawQ = ' '.join(["quast -t", threads, assembly, "-o", out_dir, "-f --circos --meta"])
            print(runRawQ)
            subprocess.call(runRawQ, shell=True)
    except FileNotFoundError:
        print('Could not find the file')
    return out_dir


def krakBrak(krak, krakdb, brak, isolate, assembly, 
                outDir, blength, krakThres, brakThres, threads):
    '''
    Function to carry out kraken classsification and bracken re-estimation
    '''
    # outputs of this function will be saved in a folder called 'Kraken'
    krakOut = os.path.join(outDir, "Kraken", isolate)
    makeDirectory(krakOut)
    # build the kraken command
    runCkrak = ' '.join([krak, "--db", krakdb, assembly, "--threads", str(threads),
            "--output", krakOut + "/All_classifications.tsv",
            "--report", krakOut + "/FULLreport.txt", "--use-names",
            "--unclassified-out", krakOut + "/unclassified.fastq",
            "--classified-out", krakOut + "/classified.fastq", 
            "--minimum-hit-groups", str(krakThres), 
            "--report-minimizer-data --confidence 0.1"])
    print(runCkrak)
    subprocess.call(runCkrak, shell = True)
    print("Kraken complete")
    # make a classic kraken command output
    runCutKrak = ' '.join(["cut -f1-3,6-8", krakOut + "/FULLreport.txt", 
                            ">", krakOut + "/ClassicFullReport.txt"])
    subprocess.call(runCutKrak, shell=True)
    # carry out bracken
    runBrak = ' '.join([brak, "-d", krakdb, "-i", krakOut + "/ClassicFullReport.txt", "-o", 
            krakOut + "/bracken_report.txt", "-t", str(brakThres), "-w",
            krakOut + "/bracken_KrakenReport.txt", "-r", str(blength)])
    print(runBrak)
    #
    subprocess.call(runBrak, shell=True)
    print("Bracken complete")
    # convert the bracken report to Krona
    kroIso = os.path.join(krakOut, isolate)
    runKrona = ' '.join(["kreport2krona.py -r ", krakOut + "/bracken_KrakenReport.txt",
                         "-o", kroIso + ".krona"])
    runKtImp = ' '.join(["ktImportText", kroIso + ".krona", "-o", kroIso + ".html"])
    print(runKrona)
    print(runKtImp)
    subprocess.call(runKrona, shell = True)
    subprocess.call(runKtImp, shell = True)
    print("Krona charts generated and saved in " + krakOut)
    # return the path to the isolate kraken folder
    return krakOut

