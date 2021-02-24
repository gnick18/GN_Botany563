#! /usr/bin/env python
#Author: Grant Nickles
import os
from pandas import read_csv
import sys
import pdb
from Bio import SeqIO

def PF05141Filterer(rootDirectory, backboneFile, annotation):
    filteringRangesAll = {}
    filteringRanges = {}
    editedRecords = []
    for genomeDirs in os.listdir(rootDirectory):
        if os.path.isdir(os.path.join(rootDirectory, genomeDirs)):
            if genomeDirs.startswith("G"):
                genomeDir = os.path.join(rootDirectory, genomeDirs)
                genome = str(genomeDirs)[:15] #this removes all but the unique identifier of the genome
                for file in os.listdir(genomeDir):
                    if os.path.isfile(os.path.join(genomeDir, file)):
                        if file.endswith('.tsv'):
                            tsvPath = os.path.join(genomeDir, file)
                            protTSV = read_csv(tsvPath, sep='\t', names=['Protein Accession', 'Sequence MD5 digest', 'Sequence length', 'Analysis', 'Signature accession', 'Signature description', 'Start location', 'Stop location', 'Score', 'Status', 'Date', 'InterPro annotations - accession', 'InterPro annotations - description', 'GO annotations', 'Pathways annotations'])
                            rows = protTSV.loc[protTSV['Signature accession'] == annotation] #this selects all rows with PF05141
                            names = rows['Protein Accession'].tolist()
                            #startStop = [rows['Protein Accession'],rows['Start location'],rows['Stop location']]
                            filteringRangesAll[genome] = []
                            filteringRanges[genome] = []
                            for index, row in rows.iterrows():
                                filteringRangesAll[genome].append([row['Protein Accession'],row['Start location'],row['Stop location']]) #this line works
                            for f in filteringRangesAll[genome]: #genomes
                                if f[0] not in (item for sublist in filteringRanges[genome] for item in sublist):
                                    filteringRanges[genome].append(f) #checks if it is anywhere in the list of list
                            if len(names) != 0:
                                for file in os.listdir(genomeDir): #find the correct order of the files in the GIO folder
                                    if os.path.isfile(os.path.join(genomeDir, file)):
                                        if file.endswith('GIO_new.txt'):
                                            faPath = os.path.join(genomeDir, file)
                                            with open(faPath) as GIO:
                                                names = [i for n, i in enumerate(names) if i not in names[:n]] #removes redundancy in the name list
                                                lines = GIO.readlines()
                                                nameIndex = {}
                                                counter = 1
                                                genomeGIO = "_".join(file.split("_")[0:2])
                                                #will give the line number and the line itself
                                                for index, line in enumerate(lines):
                                                    line=line.rstrip()
                                                    for name in names:
                                                        if name in line: #for each occurance of the annotation
                                                            nameIndex[name] = counter
                                                            counter += 1
                                                filteringRanges[genomeGIO].sort(key=lambda x: nameIndex[x[0]])

    orderedDictionary = {}
    for g in filteringRanges: #genome
        counter = 1
        for entry in filteringRanges[g]: #list of lists
            newName = g + "_" + str(counter)
            counter += 1
            orderedDictionary[newName] = [entry[1], entry[2]]
    for record in SeqIO.parse(backboneFile, "fasta"):
        #the id should be the same as the genome that is in the filteringRanges dictionary
        try:
            start = orderedDictionary[record.id][0] #start
            stop = orderedDictionary[record.id][-1] #stop
            editedRecords.append(record[start-1:stop])
        except:
            editedRecords.append(record)
            print(record.id + " needs to be filtered by hand.")
    outputPath = os.path.join("/Users/gnickles/Documents/GN_Botany563/MainProject/data", "PF05141_Euro.fa")
    SeqIO.write(editedRecords,"/Users/gnickles/Documents/GN_Botany563/MainProject/data/PF05141_Euro.fa", "fasta")


rootDirectory = r'/Volumes/HardDrive/ICSProject/Eurotiales_all'
backboneFile = r'/Users/gnickles/Documents/GN_Botany563/MainProject/data/ICS_BackboneGenes.fa'
annotation = "PF05141"
PF05141Filterer(rootDirectory, backboneFile, annotation)
