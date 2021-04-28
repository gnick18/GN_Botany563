from Bio import SeqIO
import pdb
import pandas as pd

def MakeTempNames(fastaPath):
    recordToSave = []
    tn = []
    sn = []
    gn = []
    GC = []
    counter = 1
    for record in SeqIO.parse(fastaPath, "fasta"):
        #adding the info to the key file
        tempName = str(counter) + "zzz"
        tn.append(tempName)
        counter += 1
        recordID = record.description
        geneName = recordID.split()[0]
        genome = recordID.split()[1]
        speciesName = recordID.split(",")[-2]
        sn.append(speciesName)
        gn.append(geneName)
        GC.append(genome)
        #changing to records ID to the tempName
        record.description = ""
        record.id = tempName
        recordToSave.append(record)
    df = pd.DataFrame({"TempName" : tn, "SpeciesName" : sn, "GeneName" : gn, "GenomeName" : GC})
    pdb.set_trace()
    SeqIO.write(recordToSave, "/Users/gnickles/Documents/GN_Botany563/MainProject/data/MafftAlignment/PF05141_aln_TN.fas", 'fasta')
    df.to_csv("/Users/gnickles/Documents/GN_Botany563/MainProject/data/MafftAlignment/PF05141_key.tsv", sep="\t", index=False, header = False)

fastaPath = r'/Users/gnickles/Documents/GN_Botany563/MainProject/data/MafftAlignment/PF05141_aln.fas'
MakeTempNames(fastaPath)
