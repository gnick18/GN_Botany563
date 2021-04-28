from Bio import SeqIO
import pdb
fastaPath = r'/Users/gnickles/Documents/GN_Botany563/MainProject/data/MafftAlignment/PF05141_aln.fas'
recordToSave = []
for record in SeqIO.parse(fastaPath, "fasta"):
    name = record.description
    gene = record.id
    species = name.split(",")[-1]
    if species == "":
        species = name.split(",")[-2]
    tempName = name.split()[1]
    newName = species + "__" + gene + "__" + tempName + "__"
    record.id = newName
    record.description = ""
    recordToSave.append(record)

SeqIO.write(recordToSave, "/Users/gnickles/Documents/GN_Botany563/MainProject/data/MafftAlignment/PF05141_aln_edited2.fas", 'fasta')
