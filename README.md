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

**PCount** - The total Phosphorous count of the sequence. Obtained by taking the total sequence length.

**NCount** - The Nitrogen content taken by multiplying the Amino Acid count (or the Nucleobases) by the values listed in `AminoAcids.json` or `Nucleobases.json`. More information on the json file structure can be found on my [AACount repository](https://github.com/Chonkway/AACount) page.

**CCOunt** A rough count of the Carbon in the sequence pre-translation. It counts the Carbon from each acid like NCount, but it adds the carbon from the sugar backbone.

*NOTE* : The logfile will write to a file named `results_{}` where the {} contains your parent file name. You will need to move, delete, or rename the first logfile it spits out if you want to rerun the program on that same filename.
