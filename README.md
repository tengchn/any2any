# Miscellaneous
**any2any.py** is a Python script that facilitates the conversion of sequence alignment files between different formats. It supports commonly used formats and can automatically guess the format of the input and output files based on their suffixes.

**any2any_multithreads.py** is a multithreads version of **any2any.py**, with an option "-t" or "--threads" to define the number of threads you want to use. It was designed for the large files conversion, e.g., fastq file.
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
```

**any2any_multithreads.py**
```
python3 any2any_multithreads.py -i INPUT -it INPUTTYPE -o OUTPUT -ot OUTTYPE -t THREADS

  Options:
      -i, --input: Specify the input file path.
      -it, --inputtype: Specify the input file format. If not provided, the script attempts to guess the format based on the file suffix.
      -o, --output: Specify the output file path.
      -ot, --outtype: Specify the output file format. If not provided, the script attempts to guess the format based on the file suffix.
      -t, --threads: Number of threads to use (default: 4)
  Example: python3 any2any_multithreads.py -i example.fasta -o example.phy -t 4
```
