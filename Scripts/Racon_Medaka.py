#!/usr/bin/env python3
# Script to carry out racon and medaka polishing
import os
import subprocess
import shutil

# racon first
def runRacon(isolate,assemPath,readPath,outPath,THREADS,stats):
    try:
        print('Starting racon process')
        # make the auxiliary folders to save alignments
        outAlign = os.path.join(outPath, "Alignments")
        if os.path.exists(outAlign):
            pass
        else:
            os.makedirs(outAlign)
        # make the stats folder
        statsOut = os.path.join(stats, "RaconPolishing", isolate)
        if os.path.exists(statsOut):
            pass
        else:
            os.makedirs(statsOut)
        #################### iteration 1 ####################
        print('Starting racon iteration 1')
        iter = outAlign + "/" + isolate + "_ReadsVsAssem"
        # align the reads against the raw assembly
        mini = ("minimap2 -ax map-ont", assemPath, readPath, ">", iter+"_iter1.sam", "-t", str(THREADS))
        runMini = ' '.join(mini)
        print(runMini)
        subprocess.call(runMini, shell=True)
        # call the racon aligner for racon iteration 1
        iterOut = outPath + "/" + isolate + "_iter1.fasta"
        racR = ("racon", readPath, iter+"_iter1.sam", assemPath, "-t", str(THREADS), ">", iterOut)
        runRacR = ' '.join(racR)
        print(runRacR)
        subprocess.call(runRacR, shell=True)
        # check the alignment stats
        #statsFile = statsOut + "/" +
        samSt = ("samtools flagstat --threads", str(THREADS), iter+"_iter1.sam", ">", statsOut+"/"+isolate+"_iter1_stats.txt")
        runSamSt = ' '.join(samSt)
        print(runSamSt)
        subprocess.call(runSamSt, shell=True)
        print('Iteration 1 complete')
        print('Starting racon iteration 2')
        #################### iteration 2 ####################
        mini = ("minimap2 -ax map-ont", iterOut, readPath, ">", iter+"_iter2.sam", "-t", str(THREADS))
        runMini = ' '.join(mini)
        print(runMini)
        subprocess.call(runMini, shell=True)
        # call the racon aligner for racon iteration 1
        iterOut = outPath + "/" + isolate + "_iter2.fasta"
        racR = ("racon", readPath, iter+"_iter2.sam", assemPath, "-t", str(THREADS), ">", iterOut)
        runRacR = ' '.join(racR)
        print(runRacR)
        subprocess.call(runRacR, shell=True)
        # check the alignment stats
        #statsFile = statsOut + "/" +
        samSt = ("samtools flagstat --threads", str(THREADS), iter+"_iter2.sam", ">", statsOut+"/"+isolate+"_iter2_stats.txt")
        runSamSt = ' '.join(samSt)
        print(runSamSt)
        subprocess.call(runSamSt, shell=True)
        print('Iteration 2 complete')
        print('Starting racon iteration 3')
        #################### iteration 3 ####################
        mini = ("minimap2 -ax map-ont", iterOut, readPath, ">", iter+"_iter3.sam", "-t", str(THREADS))
        runMini = ' '.join(mini)
        print(runMini)
        subprocess.call(runMini, shell=True)
        # call the racon aligner for racon iteration 1
        iterOut = outPath + "/" + isolate + "_iter3.fasta"
        racR = ("racon", readPath, iter+"_iter3.sam", assemPath, "-t", str(THREADS), ">", iterOut)
        runRacR = ' '.join(racR)
        print(runRacR)
        subprocess.call(runRacR, shell=True)
        # check the alignment stats
        #statsFile = statsOut + "/" +
        samSt = ("samtools flagstat --threads", str(THREADS), iter+"_iter3.sam", ">", statsOut+"/"+isolate+"_iter3_stats.txt")
        runSamSt = ' '.join(samSt)
        print(runSamSt)
        subprocess.call(runSamSt, shell=True)
        print('Iteration 3 complete')
        print('Starting racon iteration 4')
        #################### iteration 4 ####################
        mini = ("minimap2 -ax map-ont", iterOut, readPath, ">", iter+"_iter4.sam", "-t", str(THREADS))
        runMini = ' '.join(mini)
        print(runMini)
        subprocess.call(runMini, shell=True)
        # call the racon aligner for racon iteration 1
        iterOut = outPath + "/" + isolate + "_iter4.fasta"
        racR = ("racon", readPath, iter+"_iter4.sam", assemPath, "-t", str(THREADS), ">", iterOut)
        runRacR = ' '.join(racR)
        print(runRacR)
        subprocess.call(runRacR, shell=True)
        # check the alignment stats
        #statsFile = statsOut + "/" +
        samSt = ("samtools flagstat --threads", str(THREADS), iter+"_iter4.sam", ">", statsOut+"/"+isolate+"_iter4_stats.txt")
        runSamSt = ' '.join(samSt)
        print(runSamSt)
        subprocess.call(runSamSt, shell=True)
        print('Iteration 4 complete')
        print('Racon successfully completed')
    except FileNotFoundError:
        print('Could not find the file in the path provided. Please check again')
    return statsOut

### medaka
def runMedaka(isolate,readPath,raconOut,out_path,THREADS):
    '''Function to call medaka for just one iteration after racon. The function takes in
    the isolate name, path to the reads used to make the assembly, the path to the output 
    of racon, the path to save the outputs of medaka and the number of threads to use. 
    Note: the function is optimised for the output of 4 iterations of racon'''
    # check if the output folder exists
    try:
        medOut = os.path.join(out_path, isolate)
        if os.path.exists(medOut):
            pass
        else:
            os.makedirs(medOut)
        print('Starting medaka')
        medy = ("medaka_consensus -i", readPath, "-d", raconOut, "-o", medOut, "-t", str(THREADS))
        runMedy = ' '.join(medy)
        print(runMedy)
        subprocess.call(runMedy, shell=True)
    except FileNotFoundError:
        print('Could not find the file')
    return medOut