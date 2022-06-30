#!/usr/bin/env python3
'''
Name: PreChecks
Purpose: Script to hold functions to use for pre-checks of the 
            NanoMetaPipe pipeline script. 
Author: Damilola Oresegun
'''
import logging
import logging.handlers
import sys

def isolateList(IsoList):
    DNA_ISOLATE = []
    CDNA_ISOLATE = []
    # check the isolate names given fit the needed format
    for isolate in IsoList:
        iso = isolate.lower()
        # if the isolate has DNA in its name
        if (iso.__contains__("_dna")):
            print('You have provided DNA sequences')
            print(isolate)
            # add the isolate to the DNA list
            DNA_ISOLATE.append(isolate)
        # if the isolate has cDNA in its name
        elif (iso.__contains__("cdna")):
            print('You have provided cDNA sequences')
            print(isolate)
            # add the isolate to the cDNA list
            CDNA_ISOLATE.append(isolate)
        else:
            print('You have not provided the isolates in a satisfactory format')
            print('Do all your isolate names have _dna or _cdna or _dscdna?')
            print('Please look at the help and try again')
            sys.exit(1)
