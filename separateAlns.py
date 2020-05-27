#!/usr/bin/env python

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC,Gapped
from operator import itemgetter
import glob

# for seq_record in SeqIO.parse(glob.glob('./pan_genome_sequences/%s.fa.aln.fasta'), "fasta") % gene:

# make list of mass isolate ids
# make list of absc isolate ids

MAM_ids = open('mam_ids_FINAL.txt','r')
MAA_ids = open('maa_ids_FINAL.txt','r')


def get_ids(idFile):
    idList = []
    for index,line in enumerate(idFile):
        line = line.strip()
        isolate = line.split("\t")[0]
        idList.append(isolate)
    return idList

mass = get_ids(MAM_ids)
absc = get_ids(MAA_ids)

genes = open("absc_mass_in_at_least_1_isolate_of_both_NOCORE.tsv",'r')

for index,line in enumerate(genes):
    line = line.strip()
    line = line.split()
    if index==0:
        isolateList = line
        print line
    else:
        gene = str(line[0])
        alnName = "./pan_genome_sequences/"+gene+".fa.aln"
#        print(alnName)
        newAlnName_mass = "./mass_genes/"+gene+".fasta"
        newAlnName_absc = "./absc_genes/"+gene+".fasta"
        mass_out = open(newAlnName_mass,"w")
        absc_out = open(newAlnName_absc,"w")
        mass_sequences = []
        absc_sequences = []
        for seq_record in SeqIO.parse(alnName, "fasta"):
            #print(seq_record.id.split("_")[0])
            if str(seq_record.id.split("_")[0]) in mass:
                mass_sequences.append(seq_record)
            else:
                absc_sequences.append(seq_record)
        SeqIO.write(mass_sequences, mass_out, "fasta")
        SeqIO.write(absc_sequences, absc_out, "fasta")
        mass_out.close()
        absc_out.close()

# parse gene aln
# if isolate id in mass list
#    write that seq to a new aln named geneName_mass.fasta
# else
#    write that seq to a new aln named geneName_absc.fasta

