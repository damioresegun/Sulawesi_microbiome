#!/usr/bin/env python3
import os
import subprocess

##########################################################################
''' Check the stats of the assemblies. Starts with using assembly-stats
 for descriptive statistics.  '''
##########################################################################
# to run stats on assemblies
def run_AssemStats(isolate,in_file,out_file):
    try:
        asSt = ("assembly-stats", in_file, ">", out_file)
        runAsSt = ' '.join(asSt)
        print(runAsSt)
        subprocess.call(runAsSt, shell=True)
    except FileNotFoundError:
        print('File not found error. Check your input file again')

def raw_Quast(assembly,output_dir,reference,reference_gff,threads):
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
                print('No Quast files found in the ' +
                      'existing folder. Proceeding ' +
                      'with Quast')
                rawQ = ("quast -t", threads, assembly, "-o", out_dir, "-r", reference, "-g", reference_gff, "--large -f --circos")
                print(rawQ)
                runRawQ = ' '.join(rawQ)
                print(runRawQ)
                subprocess.call(runRawQ, shell=True)
        else:
            os.makedirs(out_dir)
            rawQ = ("quast -t", threads, assembly, "-o", out_dir, "-r", reference, "-g", reference_gff, "--large -f --circos")
            print(rawQ)
            runRawQ = ' '.join(rawQ)
            print(runRawQ)
            subprocess.call(runRawQ, shell=True)
    except FileNotFoundError:
        print('Could not find the file')
    return out_dir

