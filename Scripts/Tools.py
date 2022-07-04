#!/usr/bin/env python3
'''Script to hold various small tools and recurring commands'''

from asyncio import threads
import subprocess
import os

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


def zipFiles(directory):
    '''Function to zip files using multi-threads. Function is set up
    to take in only files, not folder. Input is the path to the 
    directory holding the files to be zipped. For each file in that
    folder, zipping will be done using pigz
    Input: path to the directory holding the files to be zipped
    Output: None
    Usage: zipFiles(path/to/directory)
    '''
    for file in os.listdir(directory,threads):
        zipDir = os.path.join(directory, file)
        runZip = ' '.join("pigz --best", zipDir, "-p", str(threads))
        print(runZip)
        subprocess.call(runZip, shell=True)
        print('FastQ files have been zipped')