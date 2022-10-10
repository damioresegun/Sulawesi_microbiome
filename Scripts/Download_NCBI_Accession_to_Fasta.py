#!/usr/bin/env python3
# Script to return the sequences from a list of 
# accession IDs
#
# Author: Damilola R Oresegun
# Adapted from: https://github.com/peterthorpe5/public_scripts/blob/master/generate_ITS1_database/bin/generate_database_direct_download.py
from Bio import SeqIO
from optparse import OptionParser
from Bio import Entrez, SeqIO

def download_accessions(accessions):
    """func takes in a list of accessions and return
    records. """
    # Acquire the NCBI records for each accession as SeqRecord objects
    # Set your email
    Entrez.email = "dro@st-andrews.ac.uk"
    idlist = ",".join(accessions)
    handle = Entrez.efetch(db="nucleotide", id=idlist, rettype="gb",
                       retmode="text")
    records = list(SeqIO.parse(handle, "genbank"))
    print("Downloaded %d GenBank records from NCBI" % len(records))
    return records


usage = """Use as follows:
$ Downloaded_NCBIA_Accession_to_Fasta.py -t tabfile.txt -o out.fasta
"""
# set the input options to look for
parser = OptionParser(usage=usage)

parser.add_option("-o", "--output", dest="out_file", default=None,
                  help="Output filename",
                  metavar="FILE")
parser.add_option("-t",  dest="tab_file", default=None,
                  help="tab_file containing the database informtation")


(options, args) = parser.parse_args()

out_file = options.out_file
accession_file = options.tab_file

(options, args) = parser.parse_args()

if __name__ == '__main__':
    accession_list = []
    for line in open(accession_file, 'r'):
        # remove the open spaces and the \n the TSV lines
        line = line.replace("\n", "").replace(" ","")
        # append the line to the list
        accession_list.append(line)
    f_out = open(out_file, 'w')
    wanted_set = set([])
    records = download_accessions(accession_list)
    print(records)
    for entry in accession_list:
        wanted_set.add(entry.split(".")[0])
    print(wanted_set)
    name_set = set([])
    for seq_record in records:
        # take the first 3 element of the description. Usaully 2 are the
        # species name. But some time 3.
        genbank_species = ("_".join(seq_record.description.split()[:3]))
        fasta_accession = seq_record.id.split(".")[0]
        if fasta_accession.split(".")[0] in wanted_set:
            name_set.add(fasta_accession.split(".")[0])
            SeqIO.write(seq_record, f_out, "fasta")
    not_found = wanted_set.difference(name_set)
    print("#DOWNLOAD INFO")
    print("Wanted number = %d" % (len(wanted_set)))
    print("we found %d" % (len(wanted_set.intersection(name_set))))