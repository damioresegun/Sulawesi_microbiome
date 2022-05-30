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
# ####### PATHS #################################################
# # get path to the reference genome
# REFPATH=config_file.get('PATHS', 'REFPATH')
# # get the path to the raw reads
# RAWPATH=config_file.get('PATHS', 'RPATH')
# # get the path to the output folder
# SAVEPATH=config_file.get('PATHS','SPATH')
# # get the path to folder holding sub-scripts
# SCPTS=config_file.get('PATHS','SCPTS')
# # get the path to the basecalled data if available. If not, will be defined later in the script
# BASEPATH=config_file.get('PATHS','DPATH')
# # get path to the reference GFF file
# REFGFF=config_file.get('PATHS','REFGFF')
# # get the path to the ONT guppy
# ONT=config_file.get('PATHS','ONT')
# # get the path to the BUSCO container
# busDOCK=config_file.get('PATHS','busDock')
# ###### get information on the experiment #################################################
# # what is the name of the experiment
# EXPMT=config_file.get('EXPERIMENT INFO', 'EXPMT')
# # what sequencing preparation was used
# KIT=config_file.get('EXPERIMENT INFO', 'KIT')
# KIT2=config_file.get('EXPERIMENT INFO', 'KIT2')
# # what barcoding library was used
# BCDLIB=config_file.get('EXPERIMENT INFO', 'BCDLIB')
# BCDLIB2=config_file.get('EXPERIMENT INFO', 'BCDLIB2')
# # what barcodes were used
# BRCDE=config_file.get('EXPERIMENT INFO', 'BCODES')
# # what flowcell was used
# FLCLL=config_file.get('EXPERIMENT INFO', 'FLOWCELL')
# #environments
# tConv=config_file.get('ENVIRONMENTS', 'tConEnv')
# THRDS=config_file.get('PROCESSING', 'THREADS')
#################################################
#################################################
### basecalling #################################################
# set the path to the basecalling script
def basecalling_function(RAWPATH,SAVEPATH,FLCLL,KIT,ONT):
    #global SCPTS
    #global bascal
    #global tConEnv
    # global RAWPATH
    # global SAVEPATH
    # global FLCLL
    # global KIT
    # global ONT
    try:
        Base_OUT = SAVEPATH+"/Basecalled"
        Log_Ouut = SAVEPATH+"/LogFiles"
        guppy = ONT+"/guppy_basecaller"
        #statss=SAVEPATH+"/Stats"
        #os.system('mkdir -p '+statss)
        os.system('mkdir -p '+Base_OUT)
        os.system('mkdir -p '+Log_Ouut)
        bascal = SCPTS+"/Basecalling.sh"
        # call the basecalling bash script
        #basecallingScript=['bash',bascal,tConv,RAWPATH,SAVEPATH,FLCLL,KIT,ONT]
        #subprocess.run(basecallingScript)
        sys.stdout = open(Log_Ouut+'/Basecalling_output.txt', 'w')
        command=('bsub "', guppy, "-i", RAWPATH, "--save_path", Base_OUT, "--flowcell", FLCLL, "--kit", KIT, "--disable_pings", "-r -v -q 0 -x auto --trim_adapters", '"')
        #myCommand = ('bsub "', program, "-i", guppy_input, "-s", out_dir.path, "--flowcell", args.flowcell, "--kit", args.kit, "--cpu_threads_per_caller", str(8), "--num_callers", str(8), '"')
        run=' '.join(command)
        print(run)
        #fin=call(run, shell=True, stdout=fo)
        #fin
        sys.stdout.close()
        print('Basecalling is done')
        print("Your raw reads are saved in: "+ RAWPATH)
        print("The results will be saved in: " + Base_OUT)
    except FileNotFoundError:
        print('There was an issue with the file. File not found')
    except SyntaxError:
        print('There was a syntax issue. Double check your syntax')
    except IndentationError:
        print('Indent Error. Check your tabs')
    except:
        print('An unknown issue occurred. Please check and try again')
    finally:
        print("Basecalling complete. Basecalling results are saved in "+SAVEPATH+"/Basecalled.")
# call the function
basecalling_function(config_file.get('PATHS', 'RPATH'), config_file.get('PATHS','SPATH'), config_file.get('EXPERIMENT INFO', 'FLOWCELL'), config_file.get('EXPERIMENT INFO', 'KIT'), config_file.get('PATHS','ONT'))