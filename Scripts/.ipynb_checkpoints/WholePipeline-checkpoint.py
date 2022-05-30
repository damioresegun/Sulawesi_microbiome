import sys
import os
import configparser
import subprocess
import argparse
### usage: python3 WholePipeline.py path/to/config_file
conf=sys.argv[1]
# read in the config file to get information on the run
config_file = configparser.ConfigParser()
config_file.read_file(open(conf))
####### PATHS #################################################
# get path to the reference genome
REFPATH=config_file.get('PATHS', 'REFPATH')
# get the path to the raw reads
RAWPATH=config_file.get('PATHS', 'RPATH')
# get the path to the output folder
SAVEPATH=config_file.get('PATHS','SPATH')
# get the path to folder holding sub-scripts
SCPTS=config_file.get('PATHS','SCPTS')
# get the path to the basecalled data if available. If not, will be defined later in the script
BASEPATH=config_file.get('PATHS','DPATH')
# get path to the reference GFF file
REFGFF=config_file.get('PATHS','REFGFF')
# get the path to the ONT guppy
ONT=config_file.get('PATHS','ONT')
# get the path to the BUSCO container
busDOCK=config_file.get('PATHS','busDock')
###### get information on the experiment #################################################
# what is the name of the experiment
EXPMT=config_file.get('EXPERIMENT INFO', 'EXPMT')
# what sequencing preparation was used
KIT=config_file.get('EXPERIMENT INFO', 'KIT')
KIT2=config_file.get('EXPERIMENT INFO', 'KIT2')
# what barcoding library was used
BCDLIB=config_file.get('EXPERIMENT INFO', 'BCDLIB')
BCDLIB2=config_file.get('EXPERIMENT INFO', 'BCDLIB2')
# what barcodes were used
BRCDE=config_file.get('EXPERIMENT INFO', 'BCODES')
# what flowcell was used
FLCLL=config_file.get('EXPERIMENT INFO', 'FLOWCELL')
#environments
tConv=config_file.get('ENVIRONMENTS', 'tConEnv')
THRDS=config_file.get('PROCESSING', 'THREADS')
#################################################
# confirm information
print("Your raw reads are saved in: "+ RAWPATH)
print("Your basecalled reads are saved in: " + BASEPATH)
os.system("mkdir -p " + SAVEPATH)
print("The results will be saved in: " + SAVEPATH)
print("Your reference genome is in: "+ REFPATH + " and its GFF is in: " + REFGFF)
print("Your experiment name is: " + EXPMT)
print("Your sequencing kit name is: " + KIT)
print("Your flowcell is: " + FLCLL)
#################################################
### basecalling #################################################
# set the path to the basecalling script
bascal=SCPTS+"/Basecalling.sh"
# call the basecalling bash script
basecallingScript=['bash',bascal,tConv,RAWPATH,SAVEPATH,FLCLL,KIT,ONT]
#print(basecallingScript)
subprocess.run(basecallingScript)
if BASEPATH=='':
    BASEPATH=SAVEPATH+"/Basecalled"
print("Basecalling is done")
###### demultiplexing #################################################
# carry out demultiplexing
demCal=SCPTS+"/Demultiplexing.sh"
demultipScrip=['bash',demCal,tConv,BASEPATH,SAVEPATH,THRDS,BCDLIB2]
subprocess.run(demultipScrip)
print("Demultiplexing done")
print("The demultiplexing outputs are saved in "+SAVEPATH+"/Demultiplexed")