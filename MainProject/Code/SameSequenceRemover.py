#! /usr/bin/env python
#Author: Grant Nickles

import os
import sys
import pdb
from Bio import SeqIO

oldFastaPath = "/Users/gnickles/Documents/GN_Botany563/MainProject/data/PF05141.fa"

RecordsToSave = []
SavedRecords = {}
for record in SeqIO.parse(oldFastaPath, "fasta"):
    species = str(record.description).split(",")[6]
    if species not in SavedRecords: #checks if key is already in the dictionary
        SavedRecords[species] = [record.seq]
        RecordsToSave.append(record)
    else:
        if str(record.seq) in SavedRecords[species]:
            continue
        else:
            SavedRecords[species].append(record.seq)
            RecordsToSave.append(record)

outputPath = r"/Users/gnickles/Documents/GN_Botany563/MainProject/data/PF05141_NoSameSeq.fa"
SeqIO.write(RecordsToSave, outputPath, 'fasta')
