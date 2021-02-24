# Notes and troubleshooting 

#### Botany 563 Main Project Pipeline

Grant Nickles; gnickles@wisc.edu

Spring of 2021 Semester: Botany 563

___

## Getting the Data

In order to get the data I had to make two scripts that consolidated everything from my main project. For this I used my GBKCreator.py pipeline grabbing 0kb around the backbone genes. I then added some of the metadata to the description park of the fasta files for later use. I'll include full copies of these two scripts that are debugged on my device bellow.

```python
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
                            if kbNumber == 0: #if it only has the backbone gene included then...
                                fastaPath = os.path.join(annoFolder, fasta) #save the path
                                pathToCopy = os.path.join(backbonePath, str(kbNumber) + "_" + genome + ".fna")
                                cmd = "cp " + fastaPath + " " + pathToCopy
                                os.system(cmd)

BackboneCopier(r'/Volumes/HardDrive/ICSProject/Eurotiales_all')
```

```python
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
            for record in SeqIO.parse(fnaPath, "fasta"):
                record.id = genome #sets the things immediatly following the > to be the temp name
                record.description = accessionDict[genome] #this sets the description in the fasta to be the species information
                allFasta.append(record)
    outputPath = os.path.join(backboneDirectory, "ICS_BackboneGenes.fa")
    SeqIO.write(allFasta, outputPath, "fasta")

bbDir = r"/Volumes/HardDrive/ICSProject/Eurotiales_all/BackboneGenes"
accession = r"/Volumes/HardDrive/ICSProject/final_accession_sorted"
BackboneFastaCombiner(bbDir, accession)
```

This gave me the final output file that is titled ICS_BackboneGenes.fa



### Big update on 02/23/2021

I thought really long and hard about the input data I was using. This is a custom dataset that only I have generated with Eurotiales. If I want to make an alignment I should be building it off the PF05141 domain only and not the other elements. Thus I need to copy and modify my GBKCreator pipeline to do this. I will copy in the code into this MainProject folder and make the needed edits.

Final output will be: all of the fasta locations that had PF05141 domain. Will NOT include the entire gene.

___

## Running alignment programs on the data

Date started: 02/22/2021

For this I will be downloading and testing out three different programs that can do alignments. I will keep all of the outputs, commands I used to run them, and will indicate roughly how long each took (min, hours? etc.) 

*NOTES: Adding things to you path*

```shell
#listing your path variable
echo $PATH 
#adding a directory to the $PATH
export PATH=[directory to add]:$PATH
#THIS ABOVE COMMAND DOES NOT MAKE IT A PERMANENT CHANGE
#To make it permanent you need to edit the .bashrc file
nano ~/.bash_profile #call from home directory
#add the export path command here
```



___

### Progam 1: ClustalW

link used to download: http://www.clustal.org/download/current/

- I downloaded the clustalw-2.1-macosx.dmg and installed it
  - all it had was the unix command `clustalw2` and the help file `clustalw_help`
- put the unix command into my bioinformatics software folder on box
- adding this to my $PATH

*Running the program: Notes from their documentation*







___

### Program 2: T-COFFEE

link used to download: http://www.tcoffee.org/Packages/Stable/Latest/

- I downloaded T-COFFEE_installer_Version_<version>_linux_x64.tar.gz
  - this said it had pre-compiled binaries for linux
- Added this to the bioinformatics software folder and put the bin in my $PATH

___

### Program 3: Muscle

link used to download: https://www.drive5.com/muscle/downloads.htm

- downloaded Mac OS X Intel i86 64 bit tar.gz file
  - When I untarred this it only gave me a single file? same thing happend on the 32 bit file

Looking at the documentation on Muscle this is correct. Just need to add the Bioinformatics folder into my $PATH so it is able to find it.

> ~*MUSCLE is a stand-alone binary*~
> Muscle is distributed as one file, known as the binary file or executable file. It is completely self-contained: it does not require configuration files, environment variables, third-party libraries or other external dependencies. There is no setup script or installer because they're not needed. To install it, all you do is download or copy the binary to a directory that is accessible from the computer where you want to run the code. For convenience, you may want to rename the binary file to muscle to avoid typing long names like muscle3.8.98_i86linux32.

___

*Running the program: Notes from their documentation*

Make an alignment and save to a file in FASTA format:

```muscle -in seqs.fa -out seqs.afa```

Write alignment to the console in CLUSTALW format (more readable than FASTA):

 ```muscle -in seqs.fa -clw```

I am opting to save the alignment in fasta format. If I end up needing it in the other format I can change that later on.

> ~*Issue that I had*~
>
> I wasn't finding the muscle command. I restarted the terminal and everything was deleted from the path??? I realized this was because I needed to add the command to the bashrc file. I will do that for all of them.
>
> I also had to allow the program to run in my security and privacy settings.

This gave me a ```Bus error: 10``` which from my research is from the amount of memory I have on my device. Likley I have some sequences that are simply too long.

- according to the program in order to fix this I need to run the `muscle -in seqs.fa -out seqs.afa -maxiters 1 -diags1 -sv` command 
  - This also might be an issue with my code 



**Command Run:**

```shell
muscle -in ./ICS_BackboneGenes.fa -out MuscleBackbone.afa
```

