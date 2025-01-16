#'/usr/bin/env python3
"""Transform frequently used alignment format to any other alignment file format."""

import sys
import re 
import os
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

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
    elif guess == "gbf" or guess == "gb":
        guess = "genbank"
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
                                     Transform frequently used alignment format to any other alignment file format. It also can guess your input based on the file suffix.
                                     """)
    parser.add_argument("-i", "--input", help="Input file")
    parser.add_argument("-it", "--inputtype", help="Input file format (default: Guess based on the file suffix).")
    parser.add_argument("-o", "--output", help="Output file")
    parser.add_argument("-ot","--outtype", help="Output filetype (default: Guess based on the file suffix).")
    args = parser.parse_args()
    if not args.input: 
        print("Please check the input file") 
        exit()
    if not args.output: 
        print("Please check the output file") 
        exit()
    return args


if __name__ == '__main__':
    args = argparser()
    any2any(args)
