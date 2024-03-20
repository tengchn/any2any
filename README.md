# Miscellaneous
**any2any.py** is a Python script that facilitates the conversion of sequence alignment files between different formats. It supports commonly used formats and can automatically guess the format of the input and output files based on their suffixes.

**any2any_parallel.py** and **any2any_multithreads.py** are multithreads version of **any2any.py**, with an option "-t" or "--threads" to define the number of threads you want to use. It was designed for the large files conversion, e.g., fastq file. **any2any_parallel.py** will process the file based on file size without regard to the content, then use the safe split points to ensure that each chunk starts and ends at a complete record boundary. If the input is not in fastq format and the file size is less than 5 GB, **any2any.py** will be used instead.

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

**any2any_parallel.py**
```
python3 any2any_parallel.py -i INPUT -it INPUTTYPE -o OUTPUT -ot OUTTYPE -t THREADS

  Options:
      -i, --input: Specify the input file path.
      -it, --inputtype: Specify the input file format. If not provided, the script attempts to guess the format based on the file suffix.
      -o, --output: Specify the output file path.
      -ot, --outtype: Specify the output file format. If not provided, the script attempts to guess the format based on the file suffix.
      -t, --threads: Number of threads to use (default: Number of CPU cores-2, prevent overload)
  Example: python3 any2any_parallel.py -i example.fasta -o example.phy
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
