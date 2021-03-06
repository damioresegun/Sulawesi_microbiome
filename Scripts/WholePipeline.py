import sys
import os
import configparser
import subprocess
import argparse
import mappy as mp
# usage: python3 WholePipeline.py path/to/config_file
conf = sys.argv[1]
# read in the config file to get information on the run
config_file = configparser.ConfigParser()
config_file.read_file(open(conf))
####### PATHS #################################################
# get path to the reference genome
REFPATH = config_file.get('PATHS', 'REFPATH')
# get the path to the raw reads
RAWPATH = config_file.get('PATHS', 'RPATH')
# get the path to the output folder
SAVEPATH = config_file.get('PATHS', 'SPATH')
# get the path to folder holding sub-scripts
SCPTS = config_file.get('PATHS', 'SCPTS')
# get the path to the basecalled data if available. If not, will be
# defined later in the script
BASEPATH = config_file.get('PATHS', 'DPATH')
# get path to the reference GFF file
REFGFF = config_file.get('PATHS', 'REFGFF')
# get the path to the ONT guppy
ONT = config_file.get('PATHS', 'ONT')
# get the path to the BUSCO container
busDOCK = config_file.get('PATHS', 'busDock')
###### get information on the experiment #################################
# what is the name of the experiment
EXPMT = config_file.get('EXPERIMENT INFO', 'EXPMT')
# what isolates were used?
ISOLATES = config_file.get('EXPERIMENT INFO', 'ISOLATES').split(',')
FILT_LENGTH = config_file.get('EXPERIMENT INFO', 'FILT_LENGTH')
FILT_QUAL = config_file.get('EXPERIMENT INFO', 'FILT_QUAL')
# what sequencing preparation was used
KIT = config_file.get('EXPERIMENT INFO', 'KIT')
KIT2 = config_file.get('EXPERIMENT INFO', 'KIT2')
# what barcoding library was used
BCDLIB = config_file.get('EXPERIMENT INFO', 'BCDLIB')
BCDLIB2 = config_file.get('EXPERIMENT INFO', 'BCDLIB2')
# what barcodes were used
BRCDES=config_file.get('EXPERIMENT INFO', 'BCODES').split(',')
# what flowcell was used
FLCLL = config_file.get('EXPERIMENT INFO', 'FLOWCELL')
######### environments #################################
tConv = config_file.get('ENVIRONMENTS', 'tConEnv')
THRDS = config_file.get('PROCESSING', 'THREADS')
#################################################
# confirm information
print("Your raw reads are saved in: " + RAWPATH)
print("Your basecalled reads are saved in: " + BASEPATH)
os.system("mkdir -p " + SAVEPATH)
print("The results will be saved in: " + SAVEPATH)
print(
    "Your reference genome is in: " +
    REFPATH +
    " and its GFF is in: " +
    REFGFF)
print("Your experiment name is: " + EXPMT)
print("Your sequencing kit name is: " + KIT)
print("Your flowcell is: " + FLCLL)
#################################################
print('The basecalling will need to be done independently on the GPU node. Make sure you have done this if you are starting from fast5 files')
###### demultiplexing #################################################
# carry out demultiplexing
demCal = SCPTS + "/Demultiplexing.sh"
demultipScrip = ['bash', demCal, tConv, BASEPATH, SAVEPATH, THRDS, BCDLIB2]
#subprocess.run(demultipScrip)
print("Demultiplexing done")
print("The demultiplexing outputs are saved in " + SAVEPATH + "/Demultiplexed")
######################################################
# STATS
######################################################
# def run_QC(SAVEPATH):
    # demultiplexed = SAVEPATH + "/Demultiplexed"
    # # run nanoStat
    # #!mkdir -p "$SAVEPATH/Stats"
    # for i in BRCDES:
        # statss = SAVEPATH + "/Stats/Demultiplexed/" + i
        # os.system('mkdir -p ' + statss)
        # file = demultiplexed + "/" + i + ".fastq.gz"
        # ofile = i + ".txt"
        # nanSt = ("NanoStat", "--fastq", file, "--outdir", statss, "-n", ofile)
        # runNanSt = ' '.join(nanSt)
        # print(runNanSt)
        # subprocess.call(runNanSt, shell=True)
        # print('nanoStat complete for ' + i)
        # print('Proceeding to nanoQC')
        # nanQ = ("nanoQC", "-o", statss, file, "-l", "50")
        # runNanQ = ' '.join(nanQ)
        # print(runNanQ)
        # subprocess.call(runNanQ, shell=True)
        # print('nanoQC complete for ' + i)
        # print('Proceeding to FastQC')
        # fatq = ("fastqc", "-t", THRDS, "-o", statss, file)
        # runFatq = ' '.join(fatq)
        # print(runFatq)
        # subprocess.call(runFatq, shell=True)
        # print('FastQC complete')
        
# run the QC scripts
#run_QC(SAVEPATH)
def run_QC(SAVEPATH,BRCDE):
    demultiplexed = SAVEPATH + "/Demultiplexed"
    # run nanoStat
    #!mkdir -p "$SAVEPATH/Stats"
    statss = SAVEPATH + "/Stats/Raw_Demultiplexed_Reads/" + BRCDE
    os.system('mkdir -p ' + statss)
    file = demultiplexed + "/" + BRCDE + ".fastq.gz"
    ofile = BRCDE + ".txt"
    nanSt = ("NanoStat", "--fastq", file, "--outdir", statss, "-n", ofile)
    runNanSt = ' '.join(nanSt)
    print(runNanSt)
    subprocess.call(runNanSt, shell=True)
    print('nanoStat complete for ' + BRCDE)
    print('Proceeding to nanoQC')
    nanQ = ("nanoQC", "-o", statss, file)
    runNanQ = ' '.join(nanQ)
    print(runNanQ)
    subprocess.call(runNanQ, shell=True)
    print('nanoQC complete for ' + BRCDE)
    print('Proceeding to FastQC')
    fatq = ("fastqc", "-t", THRDS, "-o", statss, file)
    runFatq = ' '.join(fatq)
    print(runFatq)
    subprocess.call(runFatq, shell=True)
    print('FastQC complete')
    return statss
# for i in BRCDES:
    # filtered_stats=run_QC(SAVEPATH,i)
    # print('The QC stats for the raw reads are saved in ' + filtered_stats)

def filt_qc(SAVEPATH,BRCDE,ISOLATE):
    demultiplexed = SAVEPATH + "/Demultiplexed"
    # filtering
    statss = SAVEPATH + "/Stats/Filtered_Demultiplexed_Reads/" + BRCDE + "_" + ISOLATE
    os.system('mkdir -p ' + statss)
    print('Filtering by the defined length you provided')
    print('Starting NanoFilt')
    file = demultiplexed + "/" + BRCDE + ".fastq.gz"
    filt_out = SAVEPATH + "/Filtered_RawReads"
    os.system('mkdir -p ' + filt_out)
    filt_file = filt_out + "/" + ISOLATE + ".fastq.gz"
    filtSt = ("gunzip -c", file, "| NanoFilt -l", FILT_LENGTH, "-q", FILT_QUAL, "|", "gzip >", filt_file)
    runFiltSt = ' '.join(filtSt)
    print(runFiltSt)
    #subprocess.call(runFiltSt, shell=True)
    print('Filtering complete. QC will be done for your new filtered reads')
    print('The files have been renamed according to the names provided in the config file')
    print('This means that ' + BRCDE + ' has been renamed to ' + ISOLATE)
    # qc filtered
    ofile = ISOLATES[count] + ".txt"
    nanSt = ("NanoStat", "--fastq", filt_file, "--outdir", statss, "-n", ofile)
    runNanSt = ' '.join(nanSt)
    print(runNanSt)
    subprocess.call(runNanSt, shell=True)
    print('nanoStat complete for ' + BRCDE)
    print('Proceeding to nanoQC')
    nanQ = ("nanoQC", "-o", statss, filt_file)
    runNanQ = ' '.join(nanQ)
    print(runNanQ)
    subprocess.call(runNanQ, shell=True)
    print('nanoQC complete for ' + BRCDE)
    print('Proceeding to FastQC')
    fatq = ("fastqc", "-t", THRDS, "-o", statss, filt_file)
    runFatq = ' '.join(fatq)
    print(runFatq)
    subprocess.call(runFatq, shell=True)
    print('FastQC complete')
    return filt_out, statss
# count = 0
# for i in BRCDES:
    # isola = ISOLATES[count]
    # filtered_path,filtered_stats=filt_qc(SAVEPATH,i,isola)
    # print('The raw demultiplexed reads have been successfully filtered and saved in ' + filtered_path)
    # print('The QC stats for the filtered reads are saved in ' + filtered_stats)
    # count+=1
    # print('Please remember that the files are now renamed')
    # print(i + ' is now ' + isola)
    # print('')
########################################################
#   ALIGNMENTS
########################################################
# def align(SAVEPATH):
    # global REFPATH
    # global THRDS
    # for i in BRCDES:
        # file = SAVEPATH + "/Demultiplexed/" + i + ".fastq.gz"
        # aligned_out = SAVEPATH + "/Reads_Aligned_Vs_Reference"
        # aligned_files = SAVEPATH + "/Host_Free_Reads"
        # statss = SAVEPATH + "/Stats/Alignment_Vs_Reference/" + i
        # alignCal = SCPTS + "/Alignment.sh"
        # alignScrp = ['bash', alignCal, file, aligned_out, aligned_files, THRDS, statss, REFPATH]
        # print(alignScrp)
        # subprocess.run(alignScrp)
        # #indx_path = SAVEPATH + "/Reference_index.mmi"
        # #indx = mp.Aligner(REFPATH, preset="map-ont", n_threads=int(THRDS), fn_idx_out=indx_path)
        # #if not indx: 
        # #    raise Exception("ERROR: failed to build/load index")
        # #for name, seq in mp.fastx_read(file):
        # #    print("Aligning sequence " + name + ":")
        # #    indx.map(file, n_threads=THRDS)
        # print('Alignment done for ' + i)
# align(SAVEPATH) #### not working yet
def align(REFPATH,filterd,SAVEPATH,THRDS,ISOLATE):
    file = filterd + "/" + ISOLATE + ".fastq.gz"
    aligned_out = SAVEPATH + "/Filtered_Reads_Aligned_Vs_Reference"
    aligned_files = SAVEPATH + "/Host_Free_Reads"
    statss = SAVEPATH + "/Stats/Alignment_Vs_Reference/" + ISOLATE
    alignCal = SCPTS + "/Alignment.sh"
    alignScrp = ['bash', alignCal, file, aligned_out, aligned_files, THRDS, statss, REFPATH]
    print(alignScrp)
    subprocess.run(alignScrp)
    #indx_path = SAVEPATH + "/Reference_index.mmi"
    #indx = mp.Aligner(REFPATH, preset="map-ont", n_threads=int(THRDS), fn_idx_out=indx_path)
    #if not indx: 
    #    raise Exception("ERROR: failed to build/load index")
    #for name, seq in mp.fastx_read(file):
    #    print("Aligning sequence " + name + ":")
    #    indx.map(file, n_threads=THRDS)
    print('Alignment done for ' + ISOLATE)
filtered_path="/mnt/shared/scratch/doresegu/private/JCS_MetaGenome_Project/JapaneseMacaqueData/Filtered_RawReads"
for i in ISOLATES:
    align(REFPATH,filtered_path,SAVEPATH,THRDS,i)