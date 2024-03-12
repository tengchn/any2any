#'/usr/bin/env python3
"""Transform frequently used format to any other alignment file format by using multithreads."""

import sys
import re 
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from concurrent.futures import ThreadPoolExecutor, as_completed

def process_chunk(input_chunk, input_format):
    # Read sequences from the chunk
    sequences = [seq for seq in SeqIO.parse(input_chunk, input_format)]
    return sequences

def split_sequences(input_file, input_format, chunks):
    sequences = list(SeqIO.parse(input_file, input_format))
    chunk_size = len(sequences) // chunks + (len(sequences) % chunks > 0)
    for i in range(0, len(sequences), chunk_size):
        yield sequences[i:i+chunk_size]

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
    parser = argparse.ArgumentParser(description=
                                     """
                                     Transform frequently used format to any other alignment file format by using multithreads. It also can guess your input based on the file suffix.
                                     """)
    parser.add_argument("-i", "--input", help="Input file")
    parser.add_argument("-it", "--inputtype", help="Input file format (default: Guess based on the file suffix).")
    parser.add_argument("-o", "--output", help="Output file")
    parser.add_argument("-ot","--outtype", help="Output filetype (default: Guess based on the file suffix).")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use (default: 4)")
	
    return parser.parse_args()

def main():
    args = argparser()
    input_format = args.inputtype if args.inputtype else Guess(args.input)
    output_format = args.outtype if args.outtype else Guess(args.output)
    
    # Split sequences and create temporary chunk files
    chunks = list(split_sequences(args.input, input_format, args.threads))
    chunk_files = []
    for i, chunk in enumerate(chunks):
        chunk_file = f"temp_chunk_{i}.fasta"  # Adjust temporary file format as necessary
        SeqIO.write(chunk, chunk_file, input_format)
        chunk_files.append(chunk_file)
    
    all_sequences = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(process_chunk, chunk_file, input_format): chunk_file for chunk_file in chunk_files}
        for future in as_completed(futures):
            try:
                result = future.result()
                print(f"Processed {len(result)} sequences from {chunk_file}.")
                all_sequences.extend(result)
            except Exception as exc:
                print(f"Chunk processing generated an exception: {exc}")
    
    # Write all sequences to the final output file
    SeqIO.write(all_sequences, args.output, output_format)

    # Clean up temporary chunk files
    for chunk_file in chunk_files:
        os.remove(chunk_file)

if __name__ == '__main__':
    main()

