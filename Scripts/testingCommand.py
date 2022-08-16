#!/usr/bin/env python3

itn = "/home/doresegu/scratch/private/JCS_MetaGenome_Project/MFMRCFS0322_WuM010/Host_Free_Reads/MFMRCFS0322_dscDNAVsRef_unmapped.fastq"
out = "/home/doresegu/scratch/private/JCS_MetaGenome_Project/MFMRCFS0322_WuM010/Host_Free_Reads/MFMRCFS0322_dscDNAVsRef_unmapped_renamed.fastq"

def filter_fastq_file(in_fastq, out_fastq):
    ''' The function uses Biopython to iterate through a fastq
        file and remove duplicate names before carrying out
        de novo genome assembly on Flye. Function is adapted from
        Peter Thorpe's script. The function parses the fastq file, 
        collect names in a set, if the name is not in the set, 
        write out the entry 
        Input: FastQ file to filter
        Output: Filtered FastQ file
    '''
    inp_file = open(in_fastq)
    out_file = open(out_fastq, "w")
    count = 0
    name_set = set([])
    # use enumerate to count iterations
    for i, (title, seq, qual) in enumerate(FastqGeneralIterator(inp_file)):
        # print the title of the read
        if title.rstrip("\n") not in name_set:
            # add the title to the name set
            name_set.add(title.rstrip("\n"))
            # write to file
            out_file.write("@%s\n%s\n+\n%s\n" % (title, seq.upper(), qual))
        else:
            # change the title name by adding counter
            title = title.rstrip("\n") + "_" + str(count)
            out_file.write("@%s\n%s\n+\n%s\n" % (title, seq.upper(), qual))
            #print(title + " is a duplicate. It has been changed")
            count +=1
    print("The total number of duplicated reads is: " + str(count))
    out_file.close()
    inp_file.close()

    filter_fastq_file(itn, out)