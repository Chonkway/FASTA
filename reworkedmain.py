import pathlib
import re
from typing import Sequence
from numpy import float64
import pandas as pd
from Bio import SeqIO
import os
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import json
from pandas.core.frame import DataFrame



# This FASTA sequencer by default will take in the .csv file pointed to under the fasta_file variable.
# You will need to point to your own file, and modify the columns variable according to the names you want.
# This is mostly intended to handle two columns, one with the ID and the RPKM value for gene-expression adjustments.



#Modifiable Values
csv_file_path = r""
fasta_file_path = r""
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



Seq_type = input("Specify if your sequence type is\n[1]DNA\n[2]RNA\n[3]AA\n")

Adjusted_Values = {} #Values stored here will be ready for RPKM adjustment
for seq_entry in SeqIO.parse(fasta_file_path, "fasta"):
    for ProteinID in excel_df["ProteinID"]: # //Iterates over column ProtienID
        if str(ProteinID).endswith(seq_entry.id): # //Check if iterated value matches the file
            filtered_row = excel_df.loc[excel_df["ProteinID"]==ProteinID].to_dict('list')
            RPKM = filtered_row["BE_D_RPKM-relative"][0] #//Stores sequences RPKM value
            return_values = sequence_handler_DNARNA(seq_entry) #// returns tuple of values

            if Seq_type == "1": #// No RPKM adjustment needed for DNA
                Adjusted_Values[str(ProteinID)] = list(return_values)
                
            elif Seq_type == "2" or "3": # //apply RPKM adjustment
                adjust = [float64(RPKM)*float64(i) for i in list(return_values)] #// Maybe overkill but float64 probably prevents what i assume was overflowing?
                Adjusted_Values[str(ProteinID)] = adjust

with open('Results.json', 'a',encoding="utf-8") as file:
    json.dump(Adjusted_Values, file)
print("Finished. Result file dumped.")
