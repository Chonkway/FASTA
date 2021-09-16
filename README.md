## Dependencies

You will need to install [biopython](https://biopython.org) from the linked site or by using pip:

> \>pip install biopython

Along with biopython, you will need to ensure that `AminoAcids.json` is in the same directory with the main script.


# USAGE

Supported files are the standard fasta file format:
> .fasta, .fna, .fnn, .faa, .frn ,.fa


This program uses `Bio.SeqIO` to parse the files and obtain the following information:

**PCount** - The total Phosphorous count of the sequence. Obtained by taking the total sequence length (prior to translating if you have a transcriptome).

**NCount** - The Nitrogen content taken by multiplying the Amino Acid count by the values listed in AminoAcids.json. More information on the json file can be found on my[AACount repository](https://github.com/Chonkway/AACount) page.


*NOTE* : The logfile will write to a file named `results_{}` where the {} contains your parent file name. You will need to move, delete, or rename the first logfile it spits out if you want to rerun the program on that same filename.


----
### Future plans 

I hope to include the carbon count of organims eventually, although I am unsure when I'll ever get to that.