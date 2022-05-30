#!/usr/bin/env python3
''' Script to carry out rna-seq trimming, alignment and assembly
Rationale: the NanoMetaPipe pipeline may require some cDNA
processing and this means carrying out removal of host 
transcriptome. So to do this, the user can either provide a
complete transcriptome or the paired reads of the host RNA-seq
data to generate a transcriptome. 
Requirements:
- trinity
- trimmomatic
- fastqc
- star aligner
- bwa-mem2
Note: Trimmomatic is set up for paired-end trimming.'''
''' Author: Damilola Oresegun '''
#
from ast import parse
from cProfile import run
import os
import argparse
import sys
import configparser
import subprocess
import shutil
from pathlib import Path
#
def get_args():
    parser = argparse.ArgumentParser(description="A small script to carry "+
                                     "out the processing, alignment and " +
                                     "assembly of RNA-seq data, specifically " +
                                     "cDNA reads generated from metagenomic "+
                                     "samples.", add_help=False)
    ############################################################################
    required_args = parser.add_argument_group('Required arguments')
    required_args.add_argument("-s", "--seq_type",
                                dest="Sequence_TYPE",
                                type=str, action="store",
                                choices=['reads','genome'],
                                help="Input sequence type "+
                                "can be either RNA-seq raw " +
                                "reads or can be the " +
                                "transcriptome already done",
                                required=True)
    required_args.add_argument("-o", "--output", dest="Output_DIR", 
                               default=None,
                               type=str, action="store",
                               help="Full path to Folder to place outputs",
                               required=True)
    required_args.add_argument("-ur", "--user_reads",
                                dest="User_reads",
                                action="store",
                                nargs='+',
                                help="Path to the reads to be aligned " +
                                "against the uploaded or generated " +
                                "transcriptome. Must be entered like " +
                                "-ur readsForward.fastq readsReverse.fastq " +
                                "or as a single FASTQ file",
                                required=True)
    ############################################################################
    optional_args = parser.add_argument_group('Optional arguments')
    optional_args.add_argument("-t", "--transcriptome", 
                               dest="Transcriptome_Assembly",
                               default=None,
                               type=str, action="store",
                               help="Input transcriptome assembly " +
                               "file that has already been assembled. "+
                               "Only works with [genome] sequence " +
                               "type.")
    optional_args.add_argument("-hr", "--host_reference",
                                dest="Reference",
                                action="store",
                                type=str,
                                help="Path to the reference to align " +
                                "raw reads to make the genome guided assembly")
    optional_args.add_argument("-hg", "--host_gff",
                                dest="Host_GFF",
                                action="store",
                                type=str,
                                help="Path to the GFF file of the " +
                                "host reference genome. It is not " +
                                "necessary, but will improve alignment " +
                                "against the reference genome.")
    optional_args.add_argument("-r", "--reads",
                                dest="Input_reads",
                                action="store",
                                nargs='+',
                                help="Path to forward and reverse reads "+
                                "Only works with [reads] sequence " +
                                "type e.g -r reads1.fastq reads2.fastq." +
                                "Only for generating a transcriptome assembly")
    optional_args.add_argument("-p", "--processes",
                                dest="Processes",
                                action="store",
                                type=int,
                                default=8,
                                help="Number of threads to use")
    optional_args.add_argument("-m", "--max_memory",
                                dest="Max_Memory",
                                action="store",
                                type=str,
                                default="24G",
                                help="Maximum amount of memory to use" +
                                " for this process. e.g. 32G " +
                                "Default is 24G.")
    optional_args.add_argument("-mi", "--max_intron_size",
                                dest="Max_Intron_Size",
                                action="store",
                                type=int,
                                default="50000",
                                help="Maximum intron size to use" +
                                " for genome guided assembly. Must " +
                                "correlate to your desired reference " +
                                "genome provided. Default is 50000")
    optional_args.add_argument("-i", "--isolate",
                                dest="Isolate_name",
                                action="store",
                                type=str,
                                help="The name of the isolate")
    optional_args.add_argument("-a", "--adapters",
                                dest="Adapters",
                                action="store",
                                type=str,
                                help="Full path to the adapters used for " +
                                "Illumina sequencing. Only used for 'reads' " +
                                "input sequence type")
    optional_args.add_argument("-d", "--phred",
                                dest="Phred",
                                choices=["33", "64"],
                                default="33",
                                action="store",
                                type=str,
                                help="Phred encoding used for the " +
                                "sequencing. Default is 33. Only used for 'reads' " +
                                "input sequence type")
    optional_args.add_argument("-h", "--help",
                               action="help",
                               default=argparse.SUPPRESS,
                               help="Displays this help message")
    args = parser.parse_args()
    return args
##########################################################################
# set global variables
args = get_args()
SEQ_TYP = args.Sequence_TYPE
INP_READS = args.Input_reads
INP_TRANSCRIPT = args.Transcriptome_Assembly
OUTPUTDIR = args.Output_DIR
THREADS = args.Processes
ISOLATE = args.Isolate_name
ADAPTERS = args.Adapters
PHRED = args.Phred
REFERENCE = args.Reference
HGFF = args.Host_GFF
MAXMEM = args.Max_Memory
MAXINT = args.Max_Intron_Size
UREADS = args.User_reads
# check inputs
if not INP_READS and not INP_TRANSCRIPT:
    ''' checks if the reads files or transcriptome files are provided '''
    print('You have not given any reads or transcriptome')
    print('Please check again and try')
    sys.exit(1)
else:
    # if they are present, continue
    pass
if not REFERENCE:
    '''check if reference is given to get the reference location
    if no reference then continue '''
    pass
else:
    ''' if reference is given, get the location and make a new 
    path for the index '''
    REFDIR = os.path.abspath(REFERENCE).split(".")[0] + "_STAR_index"
    print(REFDIR)
if not UREADS:
    print('You have not provided reads for alignment.')
    print('This is necessary. Please try again')
    sys.exit(1)
else:
    pass
# define a global function
def assess_Align(transcript, outputDir, u_reads, threads):
    # make assessment output
    """ transOut = os.path.join(outputDir, "Transrate")
    if os.path.exists(transOut):
        pass
    else:
        os.makedirs(transOut)
    tranR = ("transrate --assembly", transcript,
                "--left", inp_reads[0], "--right", inp_reads[1],
                "--threads", str(threads), "--output",
                transOut)
    runTranR = ' '.join(tranR)
    print(runTranR)
    subprocess.call(runTranR, shell=True) """
    # now align against the user's reads
    bwaOut = os.path.join(outputDir, "Alignment")
    if os.path.exists(bwaOut):
        pass
    else:
        os.makedirs(bwaOut)
    fname = os.path.basename(u_reads[0]).split("_")[0]
    print(fname)
    fna = os.path.join(bwaOut, fname)
    print(fna)
    # index the reference first
    bwaI = ("bwa-mem2 index", transcript)
    runBwaI = ' '.join(bwaI)
    print(runBwaI)
    subprocess.call(runBwaI, shell=True)
    bwaR = ("bwa-mem2 mem -t", str(threads), transcript,
            u_reads[0], "|", 
            "samtools view -@", str(threads), "-b -",
            "| samtools sort -o", fna+"_cDNA.bam", "-")
    runBwaR = ' '.join(bwaR)
    print(runBwaR)
    subprocess.call(runBwaR, shell=True)
    print('Alignment done')
    # removing mapped and keeping unmapped
    isolOut = os.path.join(fna, fname)
    samSt = ("samtools flagstat --threads", str(threads), 
                fna+"_cDNA.bam", ">", fna+"_mappedStats.txt")
    alnEx = fna + "_cDNA_VsTranscriptome_unmapped.bam"
    alnOut = fna + "_cDNA_VsTranscriptome.fastq"
    samEx = ("samtools view --threads", str(threads), "-f 4 -b", fna+"_cDNA.bam", ">", alnEx)
    bamFq = ("bedtools bamtofastq -i", alnEx, "-fq", alnOut)
    runSamSt = ' '.join(samSt)
    runSamEx = ' '.join(samEx)
    runBamFq = ' '.join(bamFq)
    print(runSamSt)
    subprocess.call(runSamSt, shell=True)
    print(runSamEx)
    subprocess.call(runSamEx, shell=True)
    print(runBamFq)
    subprocess.call(runBamFq, shell=True)
    print('Unmapped reads have been extracted and saved in ' + fna)
    return fna
# check if the user chose reads or genome
if SEQ_TYP == "reads":
    ''' if reads chosen, go through trimming with trimmomatic
     after trimmomatic, go through alignment and then assembly
     make output folder '''
    fastOut = os.path.join(OUTPUTDIR, "FastQC")
    if os.path.exists(fastOut):
        pass
    else:
        os.makedirs(fastOut)
    # make trimmomatic output folder
    trimOut = os.path.join(OUTPUTDIR, "Trimmomatic")
    if os.path.exists(trimOut):
        pass
    else:
        os.makedirs(trimOut)
    ISOLATES = []
    # loop through the forward and reverse reads
    for i in INP_READS:
        f_name = os.path.basename(i).split(".")[0]
        print(i)
        print(f_name)
        # run fastqc first
        fatq = ("fastqc", "-t", str(THREADS), "-o", fastOut, i)
        runFatq = ' '.join(fatq)
        print(runFatq)
        # run the fastqc command
        subprocess.call(runFatq, shell=True)
        print('FastQC complete')
        ISOLATES.append(f_name)
    # next trimmomatic - make trimmomatic output
    trimOF = os.path.join(trimOut, "isolate_summary.txt")
    fname = os.path.basename(INP_READS[0]).split("_")[0]
    trimR = os.path.join(trimOut, fname)
    if not ADAPTERS:
        print('No adapters given. General trimming will be done')
        phrr = "-phred" + PHRED
        trimq = ("trimmomatic PE -summary", trimOF,
                "-threads", str(THREADS), phrr, INP_READS[0], 
                INP_READS[1], "-baseout", trimR+".fastq", 
                "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75")
    else: 
        print('Adapter file provided. Using adapters for trimming')
        illclip = "ILLUMINACLIP:" + ADAPTERS + ":2:30:10"
        phrr = "-phred" + PHRED
        trimq = ("trimmomatic PE -summary", trimOF,
                "-threads", str(THREADS), phrr, INP_READS[0], 
                INP_READS[1], "-baseout", trimR+".fastq", illclip,
                "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:75")
    runTrimq = ' '.join(trimq)
    print(runTrimq)
    subprocess.call(runTrimq, shell=True)
    print('Trimmomatic complete')
    # rename the trimmomatic outputs
    print(trimR)
    trimmedIso = []
    for i in Path(trimOut).glob('*'):
        if str(i).endswith('_1P.fastq') or str(i).endswith('_2P.fastq'):
            print(i)
            trimmedIso.append(i)
        else:
            pass
    print('Proceeding to assembly')
    print(trimmedIso[0])
    # align against the reference genome
    ''' index the genome,
     STAR does not like gzipped entries so will need to change that here '''
    if REFERENCE.endswith('.gz'):
        print('Reference is gzipped. Will unzip to proceed')
        refG = ("pigz -d", REFERENCE, "-p"+str(THREADS))
        runRefG = ' '.join(refG)
        print(runRefG)
        subprocess.call(runRefG, shell=True)
        REFERENCE = os.path.abspath(REFERENCE).split(".gz")[0]
    else:
        pass
    if HGFF.endswith('.gz'):
        print('GFF is gzipped. Will unzip to proceed')
        refGf = ("pigz -d", HGFF, "-p"+str(THREADS))
        runRefGf = ' '.join(refGf)
        print(runRefGf)
        subprocess.call(runRefGf, shell=True)
        HGFF = os.path.abspath(HGFF).split(".gz")[0]
    else:
        pass
    if os.path.exists(REFDIR):
        pass
        if not HGFF:
            print('You have not provided a GFF file so alignment without gff will be used')
            stRu = ("STAR --runThreadN", str(THREADS), "--runMode genomeGenerate --genomeDir",
                REFDIR, "--genomeFastaFiles", REFERENCE)
        else:
            print('You have provided a GFF file so alignment will be done with gff')
            stRu = ("STAR --runThreadN", str(THREADS), "--runMode genomeGenerate --genomeDir",
                REFDIR, "--genomeFastaFiles", REFERENCE, "--sjdbGTFtagExonParentTranscript", HGFF)
    else:
        os.makedirs(REFDIR)
        if not HGFF:
            stRu = ("STAR --runThreadN", str(THREADS), "--runMode genomeGenerate --genomeDir",
                REFDIR, "--genomeFastaFiles", REFERENCE)
        else:
            print('You have provided a GFF file so alignment will be done with gff')
            stRu = ("STAR --runThreadN", str(THREADS), "--runMode genomeGenerate --genomeDir",
                REFDIR, "--genomeFastaFiles", REFERENCE, "--sjdbGTFtagExonParentTranscript", HGFF)
    runStRu = ' '.join(stRu)
    print(runStRu)
    subprocess.call(runStRu, shell=True)
    # do alignment
    print('Starting Alignment')
    pref = os.path.join(OUTPUTDIR, "Assembly")
    prefix = pref + "/" + fname
    print(prefix)
    if not HGFF:
        print('You have not provided a GFF file so alignment without gff will be used')
        stRAl = ("STAR --genomeDir", REFDIR, "--runThreadN", str(THREADS), "--readFilesIn",
                str(trimmedIso[0]), str(trimmedIso[1]), "--outFileNamePrefix", prefix, "--outSAMtype BAM",
                "SortedByCoordinate", "--outBAMsortingThreadN", str(THREADS))
    else:
        print('You have provided a GFF file so alignment will be done with gff')
        stRAl = ("STAR --genomeDir", REFDIR, "--runThreadN", str(THREADS), "--readFilesIn",
                str(trimmedIso[0]), str(trimmedIso[1]), "--outFileNamePrefix", prefix, "--outSAMtype BAM",
                "SortedByCoordinate", "--sjdbGTFtagExonParentTranscript", HGFF,
                "--outBAMsortingThreadN", str(THREADS))
    runStRA1 = ' '.join(stRAl)
    print(runStRA1)
    subprocess.call(runStRA1, shell=True)
    print('Alignment against the host genome completed')
    print('Returning reference and GFF back to zipped format')
    if REFERENCE.endswith('.gz'):
        pass
    else:
        print('zipping host reference')
        zref = ("pigz", os.path.abspath(REFERENCE).split(".gz")[0], "-p"+str(THREADS))
        runZref = ' '.join(zref)
        subprocess.call(runZref, shell=True)
    if HGFF.endswith('.gz'):
        pass
    else:
        print('zipping host gff')
        zrefg = ("pigz", os.path.abspath(HGFF).split(".gz")[0], "-p"+str(THREADS))
        runZrefg = ' '.join(zrefg)
        subprocess.call(runZrefg, shell=True)
    print('Moving on assembly')
    TrinAssem = ("Trinity --genome_guided_bam", prefix+"Aligned.sortedByCoord.out.bam", 
                "--CPU", str(THREADS), "--max_memory", MAXMEM,
                "--output", pref+"/trinity_out", 
                "--genome_guided_max_intron", str(MAXINT))
    runTrinAssem = ' '.join(TrinAssem)
    print(runTrinAssem)
    subprocess.call(runTrinAssem, shell=True)
    print('Trinity assembly complete')
    INP_TRANSCRIPT = os.path.join(pref, "trinity_out", "Trinity-GG.fasta")
    # assess the output and align against the user reads
    hf_reads = assess_Align(INP_TRANSCRIPT, OUTPUTDIR, UREADS, THREADS)
    print('cDNA alignment complete')
    #print('Host free reads are saved in ' + hf_reads)
#########################################################################################################################
elif SEQ_TYP == "genome":
    '''This is if they choose genome. Should just align genome to the reads provided'''
    print('You have only provided the transcriptome')
    print('Your generated reads will be aligned against the transcriptome')
    hf_reads = assess_Align(INP_TRANSCRIPT, OUTPUTDIR, UREADS, THREADS)
    print('cDNA alignment complete')
    #print('Host free reads are saved in ' + hf_reads)
else:
    print('Something is wrong. Please check your sequence type again')