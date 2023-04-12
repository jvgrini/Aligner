# Aligner

Aligner is a general purpose sequence alignment tool for pairwise alignments. Aligner takes a FASTA file containing two protein or nucleotide sequences and aligns them using either the Needleman-Wunsch or Smith-Waterman algorithm.

## Run Locally

Clone the project

```bash
  git clone https://github.com/jvgrini/Aligner
```

Go to the project directory

```bash
  cd aligner
```

Install dependencies

```bash
  pip install -r requirements.txt
```

## Usage/Examples

Aligner is run from its command line interface. Two positional arguments are required, namely the path of the FASTA file containing the sequences and the alignment method to be used (either 'local' or 'global'). A list of positional and optional arguments are provided below.

Performing a local alignment without any optional arguments or parameters can be done by running the following command.

```bash
python align.py protseqs.fa local
```
If one wishes to specify gap, match and mismatch values, the following command is run. Adding '-s' saves the output as a txt file named 'output.txt'.

```bash
python align.py protseqs.fa local -2 -1 4 -s
```

To use a substitution matrix, add the optional argument coresponding to the matrix of choice

```bash
python align.py protseqs.fa global -s -b62
```

Mandetory positional arguments:

```bash
sequences: path to FASTA file
method: 'local' or 'global'
```
Optional positional arguments:
```bash
gap penalty: default -1
mismatch penalty: default -1
match score: default 3
```

Positional arguments:
```bash
-h, --help            show this help message and exit
-s, --save            save output
-n, --no_output       do not print the output in the terminal
-b45, --blosum45
-b50, --blosum50
-b62, --blosum62
-b80, --blosum80
-b90, --blosum90
```
