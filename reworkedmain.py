import pathlib
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
#Global Variables

PCount = 0 # Phosphorous count will be sum of sequence length after translating, since this is when the phosphate groups show up
finalseqcount = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}
nucleobasecount = {'A': 0, 'G':0, 'C': 0, 'T': 0, 'U': 0}

with open(os.path.join(pathlib.Path(__file__).parent.resolve(),'AminoAcids.json')) as json_file: #Opens a JSON of amino acids since the translation returns protien sequence
    data = json.load(json_file)
with open(os.path.join(pathlib.Path(__file__).parent.resolve(),'Nucleobases.json')) as json_file: #Opens a JSON of amino acids since the translation returns protien sequence
    nucleobase_data = json.load(json_file)

def sequence_handler_Translate_RNA(files, filename):
    fileparse = SeqIO.parse(files, "fasta")
    bytecount = os.stat(files).st_size #Byte size for alive_bar
    with alive_bar(bytecount) as bar:
        for entry in fileparse:
            PCount = PCount + len(str(entry.seq)) #Takes the length of the entry before it it translated
            sequences = Seq(str(entry.seq)).translate() #Translates each entry into amino acids
            analyzed_seq = ProteinAnalysis(str(sequences))
            analyzed_seq.count_amino_acids() #Stores amino acid count in dict

            bar()

            for i in entry.seq: #Logs nucleobases
                if i in nucleobasecount:
                    nucleobasecount[i] = nucleobasecount[i] + 1
            bar()

            for key in finalseqcount: #Logs amino acids
                if key in analyzed_seq.amino_acids_content: #Dictionary content from .count_amino_acids()
                    finalseqcount[key] = finalseqcount[key] + analyzed_seq.amino_acids_content[key]
            bar()
    bar()


    print("Finished.")
    print("\n")

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
    for i in nucleobasecount.keys():
        for key in nucleobase_data:
            if i == key in nucleobase_data:
                CCount = CCount + (nucleobasecount[i]*6) + (nucleobasecount[i]*nucleobase_data[key][0])




def sequence_handler_RNA(files, filename):
    fileparse = SeqIO.parse(files, "fasta")
    bytecount = os.stat(files).st_size #Byte size for alive_bar
    with alive_bar(bytecount) as bar:
        for entry in fileparse:
            PCount = PCount + len(str(entry.seq)) #Takes the length of the entry before it it translated
            sequences = Seq(str(entry.seq)) #Leaves sequence as-is
            analyzed_seq = ProteinAnalysis(str(sequences))
            analyzed_seq.count_amino_acids() #Stores amino acid count in dict

            bar()

            for i in entry.seq: #Logs nucleobases
                if i in nucleobasecount:
                    nucleobasecount[i] = nucleobasecount[i] + 1
            bar()

            for key in finalseqcount: #Logs amino acids
                if key in analyzed_seq.amino_acids_content: #Dictionary content from .count_amino_acids()
                    finalseqcount[key] = finalseqcount[key] + analyzed_seq.amino_acids_content[key]
            bar()
    bar()


    print("Finished.")
    print("\n")

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
    for i in nucleobasecount.keys():
        for key in nucleobase_data:
            if i == key in nucleobase_data:
                CCount = CCount + (nucleobasecount[i]*6) + (nucleobasecount[i]*nucleobase_data[key][0])
    with open("rna_results_'{}'.txt".format(filename), "x") as results:
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



def sequence_handler_DNA(files, filename):
    fileparse = SeqIO.parse(files, "fasta")
    bytecount = os.stat(files).st_size #Byte size for alive_bar
    with alive_bar(bytecount) as bar:
        for entry in fileparse:
            PCount = PCount + len(str(entry.seq)) #Takes the length of the entry before it it translated
            sequences = Seq(str(entry.seq)) #Leaves sequence as-is
            analyzed_seq = ProteinAnalysis(str(sequences))
            analyzed_seq.count_amino_acids() #Stores amino acid count in dict

            bar()

            for i in entry.seq: #Logs nucleobases
                if i in nucleobasecount:
                    nucleobasecount[i] = nucleobasecount[i] + 1
            bar()

            for key in finalseqcount: #Logs amino acids
                if key in analyzed_seq.amino_acids_content: #Dictionary content from .count_amino_acids()
                    finalseqcount[key] = finalseqcount[key] + analyzed_seq.amino_acids_content[key]
            bar()
    bar()


    print("Finished.")
    print("\n")

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
    for i in nucleobasecount.keys():
        for key in nucleobase_data:
            if i == key in nucleobase_data:
                CCount = CCount + (nucleobasecount[i]*6) + (nucleobasecount[i]*nucleobase_data[key][0])
    with open("dna_results_'{}'.txt".format(filename), "x") as results:
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



def sequence_handler_AA(files, filename):
    fileparse = SeqIO.parse(files, "fasta")
    bytecount = os.stat(files).st_size #Byte size for alive_bar
    with alive_bar(bytecount) as bar:
        for entry in fileparse:
            sequences = Seq(str(entry.seq)) #Leaves sequence as-is
            analyzed_seq = ProteinAnalysis(str(sequences))
            analyzed_seq.count_amino_acids() #Stores amino acid count in dict

            bar()

            for i in entry.seq: #Logs nucleobases
                if i in nucleobasecount:
                    nucleobasecount[i] = nucleobasecount[i] + 1
            bar()

            for key in finalseqcount: #Logs amino acids
                if key in analyzed_seq.amino_acids_content: #Dictionary content from .count_amino_acids()
                    finalseqcount[key] = finalseqcount[key] + analyzed_seq.amino_acids_content[key]
            bar()
    bar()


    print("Finished.")
    print("\n")

    with open('AminoAcids.json') as json_file: #Opens a JSON of amino acids since the translation returns protien sequence
        data = json.load(json_file)
    with open('Nucleobases.json') as json_file: #Opens a JSON of amino acids since the translation returns protien sequence
        nucleobase_data = json.load(json_file)

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
    for i in nucleobasecount.keys():
        for key in nucleobase_data:
            if i == key in nucleobase_data:
                CCount = CCount + (nucleobasecount[i]*6) + (nucleobasecount[i]*nucleobase_data[key][0])
    with open("AA_results_'{}'.txt".format(filename), "x") as results:
        results.write("Amino Acid count:\n")
        results.write(str(finalseqcount))
        results.write("\n")
        results.write("Nitrogen count: \n")
        results.write(str(NCount))
        results.write("\n")
        results.write("Carbon count: \n")
        results.write(str(CCount))



seq_type = input("Select your sequence type:\n[1]DNA\n[2]RNA\n[3]Count AA Sequence\n")

path = input("Enter the name of the folder that contains all the FASTA files:\n")
directory = pathlib.Path(__file__).parent.resolve()
for files in os.listdir(os.path.join(directory, path)):
    if seq_type == "1":
        sequence_handler_DNA(os.path.join(directory, path, files), files)
    if seq_type == "2":
        sequence_handler_RNA(os.path.join(directory, path, files), files)
    if seq_type == "3":
        sequence_handler_AA(os.path.join(directory, path, files), files)

