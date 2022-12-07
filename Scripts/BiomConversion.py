import os
from pathlib import Path
import shutil
import subprocess
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

def makeBiom(inpt,outptr, seqtyp):
    '''
    Function to take in a folder, search through the folder to find all
    instances of a folder named 'Bracken'. If found, copy the report
    into a new folder and rename using the isolate name. Then convert 
    each report into a biom file AND make a combined biom file of all 
    isolates in the input folder
    Input: tester
    Output: tester/bracken
    Returns: String: The name of the combined Biom file
    '''
    print("Converting files to biom")
    # make the output directory if it doesnt exist
    makeDirectory(outptr)
    inptTru = os.path.join(inpt, "Kraken")
    for folder in Path(os.path.abspath(inptTru)).glob('*'):
        fnamer = os.path.basename(inpt)
        fname = fnamer.split("_")[0]
        orgTyp = os.path.basename(folder)
        outpt = os.path.join(outptr, orgTyp)
        makeDirectory(outpt)
        brakFol = os.path.join(folder, fname+"_"+seqtyp, "bracken_KrakenReport.txt")
        if os.path.exists(os.path.abspath(brakFol)):
            runErr = "Good"
            print("Copying " + brakFol + " to " + outpt)
            cpfile = os.path.join(outpt, fname+"_"+seqtyp + ".txt")
            shutil.copy2(brakFol, cpfile)
            outfile = os.path.join(outpt, fname+"_"+seqtyp + ".biom")
            # convert the txt files to biom files
            runBiom = ' '.join(["kraken-biom", cpfile, "-o", outfile])
            # print for log file
            print(runBiom)
            subprocess.call(runBiom, shell = True)
            print("Conversion completed")
        else:
            runErr = "Failed"
            print(str(folder) + " did not contain Bracken files")
            continue
    # now make a combined biom file
    if runErr == "Good":
        print("Now making combined biom file")
        allTfiles = os.path.join(outpt, "*.txt")
        allBiom = os.path.join(outpt, "CombinedIsolate.biom")
        runComBiom = ' '.join(["kraken-biom", allTfiles, "-o", allBiom])
        print(runComBiom)
        subprocess.call(runComBiom, shell = True)
    else:
        print(str(folder) + " did not contain Bracken files")
    return allBiom, runErr

# run
OUTDIR="/home/doresegu/scratch/private/JCS_MetaGenome_Project/MFMRCFS0322_WuM010"
biomOut = os.path.join(OUTDIR, "BiomFiles")
seqtyp = "DNA"
comBiom, runErr = makeBiom(OUTDIR, biomOut, seqtyp)