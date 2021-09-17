## Dependencies

You will need to install [biopython](https://biopython.org) from the linked site and [alive-progress](https://github.com/rsalmei/alive-progress) or by using pip:

> \>pip install biopython

> \> pip install alive-progress

Along with biopython, you will need to ensure that `AminoAcids.json` and `Nucleobases.json` is in the same directory with the main script.

If you want, you can delete the alive_progress wrapper and remove an indent level from everything within to avoid the need for the progress dependency. I just like it.
# USAGE

Supported files are the standard fasta file format:
> .fasta, .fna, .fnn, .faa, .frn ,.fa


This program uses `Bio.SeqIO` to parse the files and obtain the following information:

**PCount** - The total Phosphorous count of the sequence. Obtained by taking the total sequence length (prior to translating if you have a transcriptome).

**NCount** - The Nitrogen content taken by multiplying the Amino Acid count (or the Nucleobases) by the values listed in `AminoAcids.json` or `Nucleobases.json`. More information on the json file can be found on my[AACount repository](https://github.com/Chonkway/AACount) page.

**AminoAcid Count/Nucleobase Count** - The count of either Nucleobases or Amino Acids, depending on if you feed the script DNA or RNA and if you decide to translate the RNA

*NOTE* : The logfile will write to a file named `results_{}` where the {} contains your parent file name. You will need to move, delete, or rename the first logfile it spits out if you want to rerun the program on that same filename.

You can split files, although currently I am unaware how this may affect codon sequences. I might add a feature to split the files so they're always divisible by 3 but for now it might throw results off to a degree I'm uncertain of.

You can throw the program any fasta file, split the file into multiple files of any batch size (in bytes I believe it is), translate RNA to polypeptides, transcribe RNA to DNA or just leave the sequences as is.

----
### Future plans 

I hope to include the carbon count of organims eventually, although I am unsure when I'll ever get to that.

The code got messy along the way so if anything breaks I'll be sure to try and fix it, I don't know if I'll rewrite the program again though since I already did it once.
