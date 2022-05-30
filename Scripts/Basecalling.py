import sys
import os
import configparser
import subprocess
import argparse
# usage: python3 Basecalling.py path/to/config_file
# e.g. python3 ~/scratch/private/JCS_MetaGenome_Project/Scripts/Basecalling.py ~/scratch/private/JCS_MetaGenome_Project/Scripts/metagen_config.txt 

### Note: Not currently working due to an unknown error. Get help from Pete. Bash version works though
conf = sys.argv[1]
# read in the config file to get information on the run
config_file = configparser.ConfigParser()
config_file.read_file(open(conf))
RAWPATH = config_file.get('PATHS', 'RPATH')
SAVEPATH = config_file.get('PATHS', 'SPATH')
ONT = config_file.get('PATHS', 'ONT')
KIT = config_file.get('EXPERIMENT INFO', 'KIT')
FLCLL = config_file.get('EXPERIMENT INFO', 'FLOWCELL')
SCPTS = config_file.get('PATHS', 'SCPTS')

def basecalling_function(RAWPATH, SAVEPATH, FLCLL, KIT, ONT):
    try:
        Base_OUT = SAVEPATH + "/Basecalled"
        Log_Ouut = SAVEPATH + "/LogFiles"
        guppy = ONT + "/guppy_basecaller"
        os.system('mkdir -p ' + Base_OUT)
        os.system('mkdir -p ' + Log_Ouut)
        bascal=SCPTS+"/Basecalling.sh"
        basecallingScript=['bash',bascal,RAWPATH,SAVEPATH,FLCLL,KIT,ONT]
        print(basecallingScript)
        subprocess.run(basecallingScript)
        print('Basecalling is done')
        print("Your raw reads are saved in: " + RAWPATH)
        print("The results will be saved in: " + Base_OUT)
    except FileNotFoundError:
        print('There was an issue with the file. File not found')
    except SyntaxError:
        print('There was a syntax issue. Double check your syntax')
    except IndentationError:
        print('Indent Error. Check your tabs')
    except BaseException:
        print('An unknown issue occurred. Please check and try again')
    finally:
        print(
            "Basecalling complete. Basecalling results are saved in " +
            SAVEPATH +
            "/Basecalled.")
      
basecalling_function(RAWPATH, SAVEPATH, FLCLL, KIT, ONT)

# def basecalling_function(RAWPATH, SAVEPATH, FLCLL, KIT, ONT):
    # try:
        # Base_OUT = SAVEPATH + "/Basecalled"
        # Log_Ouut = SAVEPATH + "/LogFiles"
        # guppy = ONT + "/guppy_basecaller"
        # os.system('mkdir -p ' + Base_OUT)
        # os.system('mkdir -p ' + Log_Ouut)
        # #sys.stdout = open(Log_Ouut+'/Basecalling_output.txt', 'w')
        # command = (guppy, "-i", RAWPATH, "--save_path", Base_OUT, "--flowcell", FLCLL,
                   # "--kit", KIT, "--disable_pings", "-r -v -q 0 -x auto --trim_adapters")
        # run = ' '.join(command)
        # print(run)
        # #subprocess.run('bash -c', run, shell=True)
        # subprocess.run(run)
        # #sys.stdout.close()
        # print('Basecalling is done')
        # print("Your raw reads are saved in: " + RAWPATH)
        # print("The results will be saved in: " + Base_OUT)
    # except FileNotFoundError:
        # print('There was an issue with the file. File not found')
    # except SyntaxError:
        # print('There was a syntax issue. Double check your syntax')
    # except IndentationError:
        # print('Indent Error. Check your tabs')
    # except BaseException:
        # print('An unknown issue occurred. Please check and try again')
    # finally:
        # print(
            # "Basecalling complete. Basecalling results are saved in " +
            # SAVEPATH +
            # "/Basecalled.")

