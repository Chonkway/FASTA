## Dependencies

You will need to install [biopython](https://biopython.org) from the linked site) or by using pip:

> \>pip install biopython

Along with biopython, you will need to ensure that `AminoAcids.json` and `Nucleobases.json` is in the same directory with the main script.

# USAGE

Supported files are the standard fasta file format:
> .fasta, .fna, .fnn, .faa, .frn ,.fa

~~This program supports batch-file handling, although you will need to create folders that separate each fasta sequence type. 
> As an example, if you had 150 AminoAcid FASTA files you could create a folder `AA` inside the project directory and point to it for the AA option.~~

It used to do batch handling, though as it stands that functionality is gone. May be reimplemented. Really, the sequence handlers should work with whatever you want them to do, so you could rewrite the programs end loop to do what you want.

This program uses `Bio.SeqIO` to parse the files and obtain the following information:

**PCount** - The total Phosphorous count of the sequence. Obtained by taking the total sequence length. (Only for DNA/RNA) It assumes 1 phosphate to 1 base, hence the length.

**NCount** - The Nitrogen content taken by multiplying the Amino Acid count (or the Nucleobases) by the values listed in `AminoAcids.json` or `Nucleobases.json`. More information on the json file structure can be found on my [AACount repository](https://github.com/Chonkway/AACount) page.

**CCOunt** A rough count of the Carbon in the sequence pre-translation. It calculates this by taking the total # of each base, multiplying it by 6 (for the backbone ring + methyl group) and adding the native carbons of each base pair.

*NOTE* : The generated JSON logfile will not overwrite itself. The program will throw an error if you do not move the JSON file out of the programs root before running a new instance.
