#!/usr/bin/env python

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC,Gapped
from operator import itemgetter
import glob


if len(sys.argv) !=4:
    print("Usage: script.py <geneList_intermed_NOCORE.tsv> <absc_ids.txt> <mass_ids.txt>")
    sys.exit(0)

intermedNOCORE_genes = open(sys.argv[1],"r")
absc_ids = open(sys.argv[2],"r")
mass_ids = open(sys.argv[3],"r")

# for seq_record in SeqIO.parse(glob.glob('./pan_genome_sequences/%s.fa.aln.fasta'), "fasta") % gene:

def get_ids(idFile):
    idList = []
    for index,line in enumerate(idFile):
        line = line.strip()
        isolate = line.split("_")[0]
        idList.append(isolate)
    return idList

###################################################
# iterate through each sequence from each gene aln
# 	store each absc in a list
# 	store each mass in a list
# 	for each mass in the list:
# 		pair w each absc of the list
#		write fasta (with two seq)
###################################################

def make_alns(genes, absc, mass):
    for index,line in enumerate(genes):
        line = line.strip()
        line = line.split()
        if index==0:
            print line
        else:
            gene = str(line[0])
            alnName = "./pan_genome_sequences/"+gene+".fa.aln"
            newAlnName_mass = "./MAM_genes/"+gene+".fasta"
            newAlnName_absc = "./MAA_genes/"+gene+".fasta"
            mass_out = open(newAlnName_mass,"w")
            absc_out = open(newAlnName_absc,"w")
            mass_sequences = []
            absc_sequences = []
            for seq_record in SeqIO.parse(alnName, "fasta"):
                #print(seq_record.id.split("_")[0])
                if str(seq_record.id.split("_")[0]) in mass:
                    mass_sequences.append(seq_record)
                #	print(seq_record.id)
                else:
                    absc_sequences.append(seq_record)
       # print(mass_sequences)
       # print(absc_sequences)
            for massISO in mass_sequences:
                massName = str(massISO.id.split("_")[0])
                for abscISO in absc_sequences:
                    abscName = str(abscISO.id.split("_")[0])
                    tempPair = []
                    tempPair.append(massISO)
                    tempPair.append(abscISO)
                    tempAlnName = "./between_species_aln_pairs/"+gene+"-"+massName+"_"+abscName+".fasta"
                    temp_out = open(tempAlnName,"w")
                    SeqIO.write(tempPair, tempAlnName, "fasta")
                    temp_out.close()
               
            SeqIO.write(mass_sequences, mass_out, "fasta")
            SeqIO.write(absc_sequences, absc_out, "fasta")
            mass_out.close()
            absc_out.close()

absc_List = get_ids(absc_ids)
mass_List = get_ids(mass_ids)
make_alns(intermedNOCORE_genes,absc_List,mass_List)

