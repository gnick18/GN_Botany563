#! /usr/bin/env python
#Author: Grant Nickles
import os
import sys
import pdb
from Bio import SeqIO

def SeqLength(fastaPath):
    editedRecords = []
    for record in SeqIO.parse(fastaPath, "fasta"):



fastaPath = r'/Users/gnickles/Documents/GN_Botany563/MainProject/data/PF05141_NoSameSeq.fa'
