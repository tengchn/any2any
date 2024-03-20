#/usr/bin/env python3
"""Transform frequently used format to any other alignment file format by parallel processing."""

import sys
import re 
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from joblib import Parallel, delayed
import multiprocessing

def any2any(args):
    infile = args.input
    outfile = args.output
    if not args.inputtype:
        it = Guess(args.input)
    else:  
        it = args.inputtype
    if not args.outtype:
        ot = Guess(args.output)
    else:
        ot = args.outtype   
    sequences = [] 
    for seq_re in SeqIO.parse(infile, it):
        sequences.append(seq_re)      	
    SeqIO.write(sequences, outfile, ot)

def Guess(fmt):
    guess=os.path.splitext(fmt)[1][1:]
    if guess == "fas" or guess == "fa":
        guess = "fasta"
    elif guess == "fq":
        guess = "fastq"
    elif guess == "phy":
        guess = "phylip"
    elif guess == "nex":
        guess = "nexus"
    elif guess == "xml":
        guess = "seqxml"
    return guess

def argparser():
    parser = argparse.ArgumentParser(description="""
                                     Transform frequently used format to any other alignment file format by parallel processing. It also can guess your input based on the file suffix.
                                     """)
    parser.add_argument("-i", "--input", help="Input file")
    parser.add_argument("-it", "--inputtype", help="Input file format (default: Guess based on the file suffix).")
    parser.add_argument("-o", "--output", help="Output file")
    parser.add_argument("-ot","--outtype", help="Output filetype (default: Guess based on the file suffix).")
    parser.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count()-2, help="Number of threads to use (default: Number of CPU cores-2, prevent overload)")    
    return parser.parse_args()

def find_safe_split_points(file_path, num_chunks):
    file_size = os.path.getsize(file_path)
    approximate_chunk_size = file_size // num_chunks
    split_points = [0]
    with open(file_path, 'rb') as file:
        for _ in range(1, num_chunks):
            file.seek(approximate_chunk_size * _)    
            # Read until the end of the current record
            while True:
                line = file.readline()
                if line.startswith(b'@') and b' ' in line and file.tell() > split_points[-1]:
                    file.seek(-len(line), 1)
                    break
            split_points.append(file.tell())
    split_points.append(file_size)
    return split_points

def read_and_process_chunk(input_file, input_format, start_offset, end_offset):
    """Read and process a chunk of a sequence file from start_offset to end_offset."""
    processed_sequences = []
    with open(input_file, 'rb') as file:
        # Move to the start of the chunk
        file.seek(start_offset)
        # Read the chunk into memory
        chunk_data = file.read(end_offset - start_offset)
        # Assuming the file is text, decode the binary data
        chunk_text = chunk_data.decode('utf-8')
        # Use StringIO to mimic a file object for SeqIO
        from io import StringIO
        chunk_file = StringIO(chunk_text)
        for record in SeqIO.parse(chunk_file, input_format):
            processed_sequences.append(record)
    return processed_sequences

def split_and_parallel_process(input_file, input_format, num_chunks):
    """Split the input file into chunks and process them in parallel."""
    split_points = find_safe_split_points(input_file, num_chunks)
    tasks = []
    for i in range(len(split_points) - 1):
        start_pos = split_points[i]
        end_pos = split_points[i + 1]
        tasks.append((input_file, input_format, start_pos, end_pos))
    
    results = Parallel(n_jobs=num_chunks)(
        delayed(read_and_process_chunk)(*task) for task in tasks
    )
    # Combine results from all chunks
    all_sequences = [seq for result in results for seq in result]
    return all_sequences

def main():
    args = argparser()
    input_format = args.inputtype if args.inputtype else Guess(args.input)
    output_format = args.outtype if args.outtype else Guess(args.output)
    if input_format == "fastq" and os.path.getsize(args.input) > 5000000000:
        all_sequences = split_and_parallel_process(args.input, input_format, args.threads)
        SeqIO.write(all_sequences, args.output, output_format)
    else: any2any(args)
    
if __name__ == '__main__':
    main()
