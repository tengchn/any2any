# Miscellaneous
**any2any.py** is a Python script that facilitates the conversion of sequence alignment files between different formats. It supports commonly used formats and can automatically guess the format of the input and output files based on their suffixes.
## Requirements
- Python 3
- Biopython

## Usage
**any2any.py**
```
python3 any2any.py -i INPUT -it INPUTTYPE -o OUTPUT -ot OUTTYPE

  Options:
      -i, --input: Specify the input file path.
      -it, --inputtype: Specify the input file format. If not provided, the script attempts to guess the format based on the file suffix.
      -o, --output: Specify the output file path.
      -ot, --outtype: Specify the output file format. If not provided, the script attempts to guess the format based on the file suffix.
  Example: python3 any2any.py -i example.fasta -o example.phy
