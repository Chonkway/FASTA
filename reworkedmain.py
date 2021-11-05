import pathlib
import re
from types import GeneratorType
from typing import Sequence
from Bio import Seq
from numpy import float64
import pandas as pd
from Bio import SeqIO
from Bio.Seq import translate
import os
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import json
from pandas.core.frame import DataFrame
from pandas.core.indexes.base import Index



# This FASTA sequencer by default will take in the .csv file pointed to under the fasta_file variable.
# You will need to point to your own file, and modify the columns variable according to the names you want.
# This is mostly intended to handle two columns, one with the ID and the RPKM value for gene-expression adjustments.



#Modifiable Values
csv_file = input("Enter csv filename:\n")
fasta_file = input("Enter fasta file name:\n")
csv_file_path =(os.path.join(pathlib.Path(__file__).parent.resolve(), csv_file))
fasta_file_path =(os.path.join(pathlib.Path(__file__).parent.resolve(), fasta_file))
columns = ['ProteinID', 'BE_D_RPKM-relative']
excel_df = pd.read_csv(csv_file_path, sep =",", usecols=columns).replace('#VALUE!', 1) #Stores dataframe
#----------------------------------------------------------------------------------------------

finalseqcount = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}
nucleobasecount = {'A': 0, 'G':0, 'C': 0, 'T': 0, 'U': 0}

with open(os.path.join(pathlib.Path(__file__).parent.resolve(),'AminoAcids.json')) as json_file: #Opens a JSON of amino acids since the translation returns protien sequence
    amino_acid_data = json.load(json_file)
with open(os.path.join(pathlib.Path(__file__).parent.resolve(),'Nucleobases.json')) as json_file: #Opens a JSON of amino acids since the translation returns protien sequence
    nucleobase_data = json.load(json_file)

#------------------------------------------------------------------------------------------------


def sequence_handler_DNARNA(entry) -> int:
    """
    Meant to take an already parsed SeqRecord entry from a fasta. Will take in RNA or DNA fasta files.

    Returns PCount, NCount and CCount.
    """
    PCount = 0 #// Initialize PCount
    PCount = len(str(entry.seq))

    for i in entry.seq: #// Logs nucleobases
        if i in nucleobasecount: #// nucleobases for RNA/DNA
            nucleobasecount[i] = nucleobasecount[i] + 1

    NCount = 0
    for i in nucleobasecount.keys():
        for key in amino_acid_data:
            if i == key in amino_acid_data:
                NCount = NCount + (nucleobasecount[i]*nucleobase_data[key][2]) #Index 2 is the Nitrogen count of each Acid
    CCount = 0
    for i in nucleobasecount.keys():
        for key in nucleobase_data:
            if i == key in nucleobase_data:
                CCount = CCount + (nucleobasecount[i]*nucleobase_data[key][0]) #Index 0 is Carbon

    return CCount, NCount, PCount


def sequence_handler_AA(entry) -> int:
    """
    Meant to take an already parsed SeqRecord entry from a fasta.

    Returns NCount and CCount.
    """
    for i in entry.seq: #Logs nucleobases
        if i in finalseqcount:
            finalseqcount[i] = finalseqcount[i] + 1

    NCount = 0
    for i in finalseqcount.keys(): ##/ finalseqcount is for AA values
        for key in amino_acid_data:
            if i == key in amino_acid_data:
                NCount = NCount + (finalseqcount[i]*amino_acid_data[key][2]) #Index 2 is the Nitrogen count of each Acid
    CCount = 0
    for i in finalseqcount.keys():
        for key in amino_acid_data:
            if i == key in amino_acid_data:
                CCount = CCount + (finalseqcount[i]*amino_acid_data[key][0]) #Index 0 is Carbon

    return NCount, CCount

Seq_type = input("Specify if your sequence type is\n[1]DNA\n[2]RNA\n[3]AA\n[4]RNA -> AA (will run after translating)\n")

Final_Counts = {}

for entry in SeqIO.parse(fasta_file,"fasta"):
    # //This loop handles retrieving the correct return values as variable and also resetting dict values because im too lazy to fix my bad function loops
    if Seq_type == "1" or "2":
        return_values = sequence_handler_DNARNA(entry)
        nucleobasecount = dict.fromkeys(nucleobasecount, 0) #Gimpy way to reset values to 0? I guess?
    elif Seq_type == "3":
        return_values = sequence_handler_AA(entry)
        finalseqcount = dict.fromkeys(finalseqcount, 0)
    elif Seq_type == "4":
        translated_sequence = entry.translate()
        return_values = sequence_handler_DNARNA(translated_sequence)
        finalseqcount = dict.fromkeys(finalseqcount, 0)
    else:
        raise Exception("Invalid option. Please select 1-4.")

    for ProteinID in excel_df["ProteinID"]: 
        filtered_row = excel_df.loc[excel_df["ProteinID"]==ProteinID].to_dict('list')
        RPKM = filtered_row["BE_D_RPKM-relative"][0] #//Stores sequences RPKM value
        if Seq_type == "1":
            Final_Counts[str(ProteinID)] = list(return_values)
        else:
            adjust = [float64(RPKM)*float64(i) for i in list(return_values)] #// Maybe overkill but float64 probably prevents what i assume was overflowing?
            Final_Counts[str(ProteinID)] = adjust


#CNP Summative
C, N, P = 0,0,0 #Initialize values

if Seq_type == "1" or "2":
    for i in Final_Counts.values():
        C += i[0]
        N += i[1]
        P += i[2]
else:
    for i in Final_Counts.values():
        C += i[0]
        N += i[1]

Final_Counts["Total CNP"] = [C, N, P]

with open('Results_{}.json'.format(fasta_file.replace(".", "")), 'a',encoding="utf-8") as file:
    json.dump(Final_Counts, file)

print("Finished. Result file dumped.")

print(Final_Counts)
