#!/usr/bin/env python

import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC,Gapped
from operator import itemgetter
import glob
import numpy

if len(sys.argv) !=2:
    print("Usage: script.py <pangenome_btw_spp_pairs_selectionStats.txt>")
    sys.exit(0)

piFile = open(sys.argv[1],"r")

outFile = open("pangenome_between_spp_per_gene_avg_selectionStats.txt",'w')
outFile.write('Gene\tavgpi\n')

###################################################
# iterate through selectionstats file of all diff
# paired aln of each subsp. 
###################################################

geneDict = {}
geneList = []

for index,line in enumerate(piFile):
    line = line.strip()
    line = line.split()
    if index > 0:
        gene = str(line[0]).split("-")[0]
#        print(gene)
	if str(line[2]) == "None":
            pass
        else:
            pi = float(line[2]) #double check index of pi
        if gene not in geneList:
            geneList.append(gene)
            geneDict[gene] = [pi]
        else:
            geneDict[gene].append(pi)

for gene in geneDict:
    piValues = geneDict[gene]
   # print piValues
    #calculate per gene average
    perGeneAvg = numpy.mean(piValues)
    outFile.write(str(gene) + '\t' + str (perGeneAvg) + '\n')     

piFile.close()
outFile.close()
