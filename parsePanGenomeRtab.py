#!/usr/bin/env python

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC,Gapped
from operator import itemgetter

###############
# This script parses roary pangenome output to calculate gene
# frequencies to be plotted in R
##################

if len(sys.argv) !=4:
    print("Usage: script.py <gene_presence_absence.Rtab> <clade1_ids.txt> <clade2_ids.txt>")
    sys.exit(0)

Rtab = open(sys.argv[1],"r")
clade1_ids = open(sys.argv[2],"r")
clade2_ids = open(sys.argv[3],"r")

outFile = open("clade1_clade2_pangenome_gene_freq.tsv",'w')
outClade1 = open("clade1_unique_gene_freq.tsv",'w')
outClade2 = open("clade2_unique_gene_freq.tsv",'w')
outintermed = open("clade1_clade2_in_at_least_1_isolate_of_both.tsv",'w')
outNOCORE = open("clade1_clade2_in_at_least_1_isolate_of_both_NOCORE.tsv",'w')

outFile.write("gene\tclade1Freq\tclade2Freq\n")
outClade1.write("gene\tclade1Freq\tclade2Freq\n")
outClade2.write("gene\tclade1Freq\tclade2Freq\n")
outintermed.write("gene\tclade1Freq\tclade2Freq\n")
outNOCORE.write("gene\tclade1Freq\tclade2Freq\n")

def get_ids(idFile):
    idList = []
    for index,line in enumerate(idFile):
        line = line.strip()
        idList.append(line) 
    return idList

def parse(Rtab, absc, mass):
    for index,line in enumerate(Rtab):
        line = line.strip()
        line = line.split()
        if index==0:
            # initialize idx starts and stops to 0
            clade1Idx = []
            clade2Idx = []
            isolateList = line
            for item in isolateList:
                if item in absc:
                    clade1Idx.append(isolateList.index(item))
                elif item in mass:
                    clade2Idx.append(isolateList.index(item))
            print line
            print clade1Idx
            print clade2Idx
            numClade1 = len(clade1Idx)
            numClade2 = len(clade2Idx)
            print numClade1
            print numClade2
        else:
            clade1Count = 0
            clade2Count = 0
            gene = line[0]
            for index,item in enumerate(line):
                if item == "1" and index in clade2Idx:
                    clade2Count += 1
                elif item == "1" and index in clade1Idx:
                    clade1Count += 1
            geneClade1freq = float(clade1Count / float(numClade1))
            geneClade2freq = float(clade2Count / float(numClade2))
            outFile.write(gene+"\t"+str(geneClade1freq)+"\t"+str(geneClade2freq)+"\n") 
            if geneClade1freq == 0:
                outClade2.write(gene+"\t"+str(geneClade1freq)+"\t"+str(geneClade2freq)+"\n")
            elif geneClade2freq == 0:
                outClade1.write(gene+"\t"+str(geneClade1freq)+"\t"+str(geneClade2freq)+"\n")
            else:
                outintermed.write(gene+"\t"+str(geneClade1freq)+"\t"+str(geneClade2freq)+"\n")
                if geneClade1freq != 1 and geneClade2freq != 1:
                    outNOCORE.write(gene+"\t"+str(geneClade1freq)+"\t"+str(geneClade2freq)+"\n")

clade1_List = get_ids(clade1_ids)
clade2_List = get_ids(clade2_ids)

parse(Rtab, clade1_List, clade2_List)

