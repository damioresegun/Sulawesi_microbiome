#!/usr/bin/env python3
'''
Name: PreChecks
Purpose: Script to hold functions to use for pre-checks of the 
            NanoMetaPipe pipeline script. 
Author: Damilola Oresegun
'''
import os

def isolateList(IsoList):
    DNA_ISOLATE = []
    CDNA_ISOLATE = []
    # check the isolate names given fit the needed format
    for isolate in IsoList:
        iso = isolate.lower()
        # if the isolate has DNA in its name
        if (iso.__contains__("_dna")):
            #print('You have provided DNA sequences')
            #print(isolate)
            # add the isolate to the DNA list
            DNA_ISOLATE.append(isolate)
            DNAPRES = True
        # if the isolate has cDNA in its name
        elif (iso.__contains__("cdna")):
            #print('You have provided cDNA sequences')
            #print(isolate)
            # add the isolate to the cDNA list
            CDNA_ISOLATE.append(isolate)
            CDNAPRES = True
        else:
            DNAPRES = False
            CDNAPRES = False
            #print('You have not provided the isolates in a satisfactory format')
            #print('Do all your isolate names have _dna or _cdna or _dscdna?')
            #print('Please look at the help and try again')
            #sys.exit(1)
    return DNA_ISOLATE, CDNA_ISOLATE, DNAPRES, CDNAPRES


def seqCheck(seqType, makcref, creads, reference, cref, cadap):
    noCref = False 
    CrefMake = False
    Cref = False
    noRef = False
    noAdap = False
    Adap = False
    makCref = True
    AllCheck = False
    NoSeqType = False
    if seqType == "dna":
        pass
    elif seqType == "cdna":
    # if the sequence type is cdna, check if the user wants to make a transcriptome
        if makcref is True:
            # if the user wants to make a transcriptome, check if they provide reads
            if not creads:
                noCref = True
            else:
                Cref = True
        else:
            # if the user is not making a transcriptome, make the given reference the cdna transcriptome
            CrefMake = True
            pass
    # if the user choose both dna and cda, check inputs
    elif seqType == "both":
        # check if the dna reference is given
        if not reference:
            noRef = True
        # check if the cdna reference is given
        elif not cref:
            # check if the transcriptome assembly is to be made
            if makcref is True:
                if not creads:
                    noCref = True
                else:
                    Cref = True
            # check if the adaptors are given for the cdna
                if not cadap:
                    noAdap = True
                else:
                    Adap = True
            else:
                makCref = False
        else:
            AllCheck = True
    else:
        NoSeqType = True
    return noCref, Cref, CrefMake, noRef, noAdap, Adap, makCref, AllCheck, NoSeqType


def filterOptions(demulp,args,filter):
    # if its qcat
    if demulp == "qcat":
        G_KIT = args.Sequencing_Kit
        # set the correct 
        Q_KIT = "NBD103/" + args.Expansion_Kit
    # if its guppy
    else:
        G_KIT = "SQK-" + args.Sequencing_Kit
        Q_KIT = "EXP-" + args.Expansion_Kit
    # check if the user chose to carry out filtering
    if filter is True:
        # set the dna filter length
        DNA_FILT_LENGTH = args.dna_Filter_Length
        CDNA_FILT_LENGTH = args.cDNA_Filter_Length
        # set the quality
        FILT_QUAL = args.Filter_Quality
        # set bracken length to the filter length
        DBRACK_LENGTH = DNA_FILT_LENGTH
        CBRACK_LENGTH = CDNA_FILT_LENGTH
    else:
        # set bracken length to the default filter length
        DBRACK_LENGTH = args.dna_Filter_Length
        CBRACK_LENGTH = args.cDNA_Filter_Length
        pass
    return G_KIT, Q_KIT, DBRACK_LENGTH, CBRACK_LENGTH, FILT_QUAL