from Bio import SeqIO
import pdb
fastaPath1 = r'/Users/gnickles/Documents/GN_Botany563/MainProject/data/MafftAlignment/PF05141_aln_trimmed10.fas'
fastaPath2 = r'/Users/gnickles/Documents/GN_Botany563/MainProject/data/MafftAlignment/PF05141_aln_trimmedGappyout.fas'
recordToSave1 = []
recordToSave2 = []
for record in SeqIO.parse(fastaPath1, "fasta"):
    name = record.id
    nameSplit = name.split("__")
    if len(nameSplit) == 4:
        part1 = nameSplit[0]
        part2 = nameSplit[1]
        part3 = nameSplit[2]
        newName = part1 + ":" + part2 + "(" + part3 + ")"
        record.description = ""
        record.id = newName
    recordToSave1.append(record)

for record in SeqIO.parse(fastaPath2, "fasta"):
    name = record.id
    nameSplit = name.split("__")
    if len(nameSplit) == 4:
        part1 = nameSplit[0]
        part2 = nameSplit[1]
        part3 = nameSplit[2]
        newName = part1 + ":" + part2 + "(" + part3 + ")"
        record.description = ""
        record.id = newName
    recordToSave2.append(record)

SeqIO.write(recordToSave1, "/Users/gnickles/Documents/GN_Botany563/MainProject/data/MafftAlignment/PF05141_aln_editedTrimmed10.fas", 'fasta')
SeqIO.write(recordToSave2, "/Users/gnickles/Documents/GN_Botany563/MainProject/data/MafftAlignment/PF05141_aln_editedTrimmedGappyout.fas", 'fasta')
