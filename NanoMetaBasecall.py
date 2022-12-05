#!/usr/bin/env python3
#
# NanoMetaBasecall v1
#
# Short script to carry out basecalling on a GPU node of a cluster
# Rationale: Development took place on a large cluster which required
# requesting GPU access and resources. Hence this step is sectioned off
# in order to curtail the need to use the GPU for the entire pipeline
#
# MacKenzie Institute for Early Diagnosis (2022)
# Author: Damilola Oresegun, Peter Thorpe
#

import os
import argparse
from matplotlib
#####################################################################################################
def get_args():
    parser = argparse.ArgumentParser(description="Pipeline for basecalling " +
                                     "using the guppy basecaller. Requires " +
                                     " a GPU to carry out basecalling", 
                                     add_help=False)
    file_directory = os.path.realpath(__file__).split("NanoMetaBasecall")[0]
    if not os.path.isfile(os.path.join(file_directory, "NanoMetaBasecall.py")):
        file_directory = os.path.join(file_directory, "NanoMetaBasecall")
    if not os.path.isfile(os.path.join(file_directory, "NanoMetaBasecall.py")):
        print("Can't locate the correct path to the basecalling script")
    #################################################################################################
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument("-p", "--program",
                               dest="Guppy_basecaller",
                               action="store",
                               type=str,
                               default="guppy_basecaller",
                               help="Full path to the guppy_basecaller " +
                               "e.g. path/to/ont_guppy_v6.0.1/bin/guppy_basecaller" +
                               "program or if the bin folder is in the PATH " +
                               "just guppy_baseller. Default: guppy_basecaller",
                               required=True)
    required_args.add_argument("-r", "--raw",
                               dest="Raw_reads",
                               action="store",
                               type=str,
                               help="Full path to the fast5 folder of the " +
                               "raw sequenced reads. No defaults.",
                               required=True)
    required_args.add_argument("-o", "--output_directory",
                               dest="Output_directory",
                               action="store",
                               type=str,
                               help="Full path to the folder to save outputs",
                               required=True)
    required_args.add_argument("-m", "--model", 
                                dest="Model", 
                                action="store",
                                choices=["hac", "sup"],
                                type=str,
                                default="hac",
                                help="The model to use for basecalling. HAC is "+
                                "the high accuracy model that has been in used " +
                                "a few years now. It is the default model in " +
                                "in both guppy and this pipeline. SUP is the " +
                                "super accuracy model which is more recent and " +
                                "marginally more accurate than HAC however will " +
                                "take considerably longer to complete.",
                                required=True)
    #################################################################################################
    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument("-k", "--kit", 
                               dest="kit",
                               action="store",
                               type=str,
                               default="SQK-LSK109",
                               help="The sequencing kit used for this experiment")
    optional_args.add_argument("-f", "--flowcell",
                               dest="Flowcell",
                               action="store",
                               type=str,
                               default="FLO-MIN106",
                               help="The flowcell used for the experiment")
    args = parser.parse_args()
    return args, file_directory
#####################################################################################################
# set global variables
args, FILE_DIRECTORY = get_args()
RAWPATH = args.Raw_reads
SAVEPATH = args.Output_directory
ONT = args.Guppy_basecaller
KIT = args.kit
FLCLL = args.Flowcell
MODEL = args.Model

def basecalling_function(raw, out, flow, kit, guppy):
    try:
        Base_OUT = os.path.join(out, "Basecalled")
        if os.path.exists(Base_OUT):
            pass
        else:
            os.makedirs(Base_OUT)
        if MODEL == "hac":
            baseScp = (guppy, "-i", raw, "--save_path", Base_OUT, "--flowcell", flow,
                    "--kit", kit, "--disable_pings -r -v -q 0 -x auto",
                    "--trim_adapters --gpu_runners_per_device 12",
                    "--chunks_per_runner 1024 --num_callers 8")
            runBaseScp = ' '.join(baseScp)
            print(runBaseScp)
            #subprocess.call(runBaseScp, shell=True)
        elif MODEL == "sup":
            baseScp = (guppy, "-i", raw, "--save_path", Base_OUT, "--config", "dna_r9.4.1_450bps_sup.cfg", "--flowcell", flow,
                    "--kit", kit, "--disable_pings -r -v -q 0 -x auto",
                    "--trim_adapters --gpu_runners_per_device 12",
                    "--chunks_per_runner 1024 --num_callers 8")
            runBaseScp = ' '.join(baseScp)
            print(runBaseScp)
            #subprocess.call(runBaseScp, shell=True)
        else:
            print('An appropriate model was not provided. Please try again and' +
            'choose between hac and sup.')
        print('Basecalling is done')
        print("The results are be saved in: " + Base_OUT)
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

