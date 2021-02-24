#! /usr/bin/env python
#Author: Grant Nickles

import os
import sys
import pdb
#using this program for the Botany 563 class. I want a dataset of all the backbone genes to practice making trees with.
def BackboneCopier(rootDirectory):
    #making all of the window folders in the root directory
    backbonePath = os.path.join(rootDirectory, "BackboneGenes")
    os.mkdir(backbonePath)
    for genomeFolders in os.listdir(rootDirectory): #stepping through the genome folders
        if os.path.isdir(os.path.join(rootDirectory, genomeFolders)) and str(genomeFolders).endswith("folder"):
            genomeFolder = os.path.join(rootDirectory, genomeFolders)
            genome = str(genomeFolders)[:15] #this removes all but the unique identifier of the genome
            for annotationFolders in os.listdir(genomeFolder): #stepping into the Annotation_Windows folder inside the genome folder
                if os.path.isdir(os.path.join(genomeFolder, annotationFolders)):
                    annoFolder = os.path.join(genomeFolder, annotationFolders)
                    for fasta in os.listdir(annoFolder):
                        if os.path.isfile(os.path.join(annoFolder, fasta)) and str(fasta).endswith(".fna"): #finding all of the zero kb files
                            kbNumber = int(fasta.split("_")[1][:-2]) #taking the number and not the kb ex 10kb = 10
                            hitNumber = fasta.split("_")[0]
                            if kbNumber == 0: #if it only has the backbone gene included then...
                                fastaPath = os.path.join(annoFolder, fasta) #save the path
                                pathToCopy = os.path.join(backbonePath, hitNumber + "_" + str(kbNumber) + "_" + genome + ".fna")
                                cmd = "cp " + fastaPath + " " + pathToCopy
                                os.system(cmd)

BackboneCopier(r'/Volumes/HardDrive/ICSProject/Eurotiales_all')
