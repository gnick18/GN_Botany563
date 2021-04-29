import pdb
import pandas as pd
import massedit

#assumes the key files has headers, code can be modified to the index of the columns if needed
def TempNamesToKeyNames(tempNameColumn, keyNameColumn, fileChanging, keyFile):
    df = pd.read_csv(keyFile, sep='\t', header=0) #reading in the key file
    #looping over the rows in the DataFrame from the bottom up (prevents 1zzz chaning 12(1zzz) for example
    for i, row in df[::-1].iterrows():
        tempName = row[tempNameColumn]
        keyName = row[keyNameColumn]
        print("Chaning " + tempName + " to " + keyName)
        #in order to pass variables in a reg expression you need to add the rf to the begining of the string
        #you also need to put the variable in '{}' with quotes around it
        massedit.edit_files(fileChanging, [rf"re.sub('{tempName}', '{keyName}', line)"], dry_run=False)


tempNameColumn = "TempName"
keyNameColumn = "SpeciesName"
keyFile = r'/Users/gnickles/Documents/GN_Botany563/MainProject/data/iqtree/PF05141_key.tsv'
#fileChanging = [r'/Users/gnickles/Documents/GN_Botany563/MainProject/data/iqtree/gappyout/PF05141_trimGappyout.contree', r'/Users/gnickles/Documents/GN_Botany563/MainProject/data/iqtree/gappyout/PF05141_trimGappyout.treefile']
fileChanging = [r'/Users/gnickles/Documents/GN_Botany563/MainProject/data/iqtree/PF05141_aln_trimmedGappyout.fas']
TempNamesToKeyNames(tempNameColumn, keyNameColumn, fileChanging, keyFile)
