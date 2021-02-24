#! /usr/bin/env python
#Author: Grant Nickles
from Bio import SeqIO
import os
import sys
import pdb
import pandas as pd
#This program will take all the fasta files and combine them to one file, adding on the GCA name as it does so.
def BackboneFastaCombiner(backboneDirectory, accession):
    #loading in the two files from their location on the computer
    accessionDF = pd.read_csv(accession, sep="\t", header=None)
    SpeciesInfo = accessionDF.iloc[:, -1] #seleting the last column in the DF and storing as series
    TempName = accessionDF.iloc[:,0]
    TempAndSpecies = pd.merge(TempName, SpeciesInfo, right_index = True, left_index = True)
    TempAndSpecies.columns = ["TempName", "SpeciesName"]
    accessionDict = {}
    for index, row in TempAndSpecies.iterrows(): #saving the accession data to a dictionary for easy access later on
        accessionDict[row["TempName"]] = row["SpeciesName"].strip()
    #going through all of the fasta files
    allFasta = []
    for fasta in os.listdir(backboneDirectory):
        if os.path.isfile(os.path.join(backboneDirectory, fasta)) and str(fasta).endswith(".fna"):
            fnaPath = os.path.join(backboneDirectory, fasta)
            genome = ".".join("_".join(fnaPath.split("_")[-2:]).split(".")[0:2]) #extracting only the genome name
            annotationNumber = fasta.split("_")[0]
            for record in SeqIO.parse(fnaPath, "fasta"):
                record.id = genome + "_" + annotationNumber #sets the things immediatly following the > to be the temp name
                record.description = accessionDict[genome] #this sets the description in the fasta to be the species information
                allFasta.append(record)
    outputPath = os.path.join(backboneDirectory, "ICS_BackboneGenes.fa")
    SeqIO.write(allFasta, outputPath, "fasta")

bbDir = r"/Volumes/HardDrive/ICSProject/Eurotiales_all/BackboneGenes"
accession = r"/Volumes/HardDrive/ICSProject/final_accession_sorted"
BackboneFastaCombiner(bbDir, accession)
