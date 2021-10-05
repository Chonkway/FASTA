import re
from typing import Sequence
from Bio import SeqIO
import os
from Bio.Seq import Seq, translate
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import json
from alive_progress import alive_bar; import time, logging


print("--------------") #
filename = input("Enter your filename(including extension). Ensure it is in the root directory.") #Sets target file
fileparse = SeqIO.parse(filename, "fasta")
bytecount = os.stat(filename).st_size #Byte size for alive_bar

PCount = 0 # Phosphorous count will be sum of sequence length after translating, since this is when the phosphate groups show up

finalseqcount = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}
nucleobasecount = {'A': 0, 'G':0, 'C': 0, 'T': 0, 'U': 0}

with alive_bar(bytecount,force_tty=True) as bar:

    for entry in fileparse:
        PCount = PCount + len(str(entry.seq)) #Takes the length of the entry before it it translated
        sequences = Seq(str(entry.seq)).translate() #Translates each entry into amino acids
        analyzed_seq = ProteinAnalysis(str(sequences))
        analyzed_seq.count_amino_acids() #Stores amino acid count in dict

        for i in entry.seq: #Logs nucleobases
            if i in nucleobasecount:
                nucleobasecount[i] = nucleobasecount[i] + 1

        for key in finalseqcount: #Logs amino acids
            if key in analyzed_seq.amino_acids_content: #Dictionary content from .count_amino_acids()
                finalseqcount[key] = finalseqcount[key] + analyzed_seq.amino_acids_content[key]
bar()


print("Finished.")
print("\n")

with open('AminoAcids.json') as json_file: #Opens a JSON of amino acids since the translation returns protien sequence
    data = json.load(json_file)

print("Your final amino acid counts are: \n")
print(finalseqcount)

print("Your final nucleobase counts are: \n") #Nucleobases for RNA/DNA
print(nucleobasecount)


NCount = 0 #Counts nitrogens
for i in finalseqcount.keys():
    for key in data:
        if i == key in data:
            NCount = NCount + (finalseqcount[i]*data[key][2]) #Index 2 is the Nitrogen count of each Acid
CCount = 0

with open("results_'{}'.txt".format(filename), "x") as results: #Logs results because I'm sick of waiting for console
    results.write("Amino acid count:\n")
    results.write(str(finalseqcount))
    results.write("\n")
    results.write("Nucleobase count:\n")
    results.write(str(nucleobasecount))
    results.write("\n")
    results.write("Phosphorous count: \n")
    results.write(str(PCount))
    results.write("\n")
    results.write("Nitrogen count: \n")
    results.write(str(NCount))
    results.write("\n")
    results.write("Carbon count: \n")
    results.write(str(CCount))