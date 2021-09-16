import re, cProfile
from typing import Sequence
from Bio import SeqIO
import os
from Bio.AlignIO import parse, write
from Bio.Seq import Seq, translate
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import json
from alive_progress import alive_bar; import time, logging
#Reimpliement batch iteration -

def batch_iterator(iterator, batch_size):
    """
    Returns multiple lists of length batch_size.
    """

    entry = True #Ensures the loop runs once
    while entry:
        batch = []
        while len(batch) < batch_size: 
            try:
                entry = next(iterator)
            except StopIteration:
                entry = None
            if entry is None:
                break#End of file
            batch.append(entry)
        if batch:
            yield batch

print("--------------")
filename = input("Enter your filename(including extension). Ensure it is in the root directory.") #Sets target file
query = input("If the file is significantly large, it is recommended that you split it. Would you like to split the file?(y/n)")
seqtype = input("Does your file need transcribed?(y/n)") #Used for a check below to use SeqIO's .translate() module

if query.lower() == "y": #Checks if you want/need to split the file
    record_iter = SeqIO.parse(open(filename), "fasta") #Uses SeqIO to parse the file as the iterator (subject to change to a different fasta parser)
    for i, batch in enumerate(batch_iterator(record_iter, 60000)):
        file = "group_%i.fasta" % (i + 1)
        with open(file, "w") as handle:
            count = SeqIO.write(batch, handle, "fasta")
        print("Wrote %i records to %s" % (count, file))
        
with alive_bar() as bar:
    print("---------------------")
    print("\n \n \n \n \n")

    #Scan directory for common fasta file extensions using OS module -
    ext = ('.fasta', '.fna', '.fnn', '.faa', '.frn' , '.fa')

    PCount = 0 # Phosphorous count will be sum of sequence length pre-translation, creates variable.
    finalseqcount = {'A': 0, 'C': 0, 'D': 0, 'E': 0, 'F': 0, 'G': 0, 'H': 0, 'I': 0, 'K': 0, 'L': 0, 'M': 0, 'N': 0, 'P': 0, 'Q': 0, 'R': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0}

    #Loop over all files, count their Amino Acids and create a total - 
    print("Now parsing files. Depending on the number of files and their size, this may take a while...")

    if query == 'y':
        for files in os.listdir(): #Scan for all files (applicable only if files are split)
            if files.endswith(ext):
                if files != filename:
                    fileparse = SeqIO.parse(files, "fasta") #Begins the file parse. Trying to keep this as the only parse step to avoid costly time
        
    else:
        fileparse = SeqIO.parse(filename, "fasta") #If there's only one file used, no need to scan directory

    for entry in fileparse: #Depending on if the file is large/a large series of a lot of smaller files, this can take a long time
        PCount = PCount + len(str(entry.seq))

        if seqtype.lower() == "y": # Check for mRNA
            mRNA_translate = Seq(str(entry.seq)).translate()
            analyzed_seq = ProteinAnalysis(str(mRNA_translate)) #Allows .count_amino_acids() to apply to the files
            


        else:
            analyzed_seq = ProteinAnalysis(str(entry.seq))

        aacount = analyzed_seq.count_amino_acids()

        for key in finalseqcount:
            if key in analyzed_seq.amino_acids_content: #Dictionary content from .count_amino_acids()
                finalseqcount[key] = finalseqcount[key] + analyzed_seq.amino_acids_content[key]
    else:
        pass
    bar()
print("Finished.")
print("\n")

print("Your final amino acid counts are: \n")
print(finalseqcount)
log = input("Would you like to write these to a log file?(y/n)")

with open('AminoAcids.json') as json_file:
    data = json.load(json_file)

NCount = 0
for i in finalseqcount.keys():
    for key in data:
        if i == key in data:
            NCount = NCount + (finalseqcount[i]*data[key][2])
            
if log.lower() == "y":   
    with open("results_'{}'.txt".format(filename), "x") as results: #Logs results because I'm sick of waiting for console
        results.write("Amino acid count:\n")
        results.write(str(finalseqcount))
        results.write("\n")
        results.write("Phosphorous count: \n")
        results.write(str(PCount))
        results.write("\n")
        results.write("Nitrogen count: \n")
        results.write(str(NCount))


