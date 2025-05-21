#!/usr/bin/env python3

from Bio.Seq import Seq
import argparse
import numpy as np
from tqdm.auto import tqdm
import os

try:
    from utils import MUMdata, get_seq_paths, get_sequence_lengths
except ImportError:
    from mumemto.utils import MUMdata, get_seq_paths, get_sequence_lengths

def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description='Extract the MUM sequences')
    parser.add_argument('-m', '--mumfile', type=str, help='bumbl file to process')
    parser.add_argument('-i', '--index', type=int, default=0, help='The index of the file in the corresponding filelist')
    parser.add_argument('-o', '--output', type=str, help='The name of the output file')
    parser.add_argument('-f', '--filelist', type=str, help='filelist file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print verbose output')
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)
    
    if args.filelist == None:
        args.filelist = os.path.splitext(args.mumfile)[0] + '.lengths'
        if not os.path.exists(args.filelist):
            raise FileNotFoundError(f"Filelist {args.filelist} not found, and no filelist provided")
    
    if args.output is None:
        args.output = os.path.splitext(args.mumfile)[0] + '_mums.fa'
    if not args.output.endswith('.fa') and not args.output.endswith('.fasta'):
        args.output += '.fa'
    return args

def main(args):
    file = get_seq_paths(args.filelist)[args.index]
    lengths = get_sequence_lengths(args.filelist)
    if not os.path.exists(file):
        raise FileNotFoundError(f"File {file} not found. Ensure lengths file is formatted correctly. Perhaps the paths are relative, not absolute?")
    seq = Seq(''.join([l for l in open(file, 'r').read().splitlines() if not l.startswith('>')]))
    assert len(seq) == lengths[args.index], f"Sequence length {len(seq)} does not match expected length {lengths[args.index]}. There may be Ns or other non-ACGT chars in the FASTA file."
    mums = MUMdata(args.mumfile, sort=False)
    mum_lengths, mum_starts, mum_strands = mums.lengths, mums.starts, mums.strands
    with open(args.output, 'w') as out:
        outlines = []
        for i in tqdm(range(len(mum_lengths)), desc="Extracting MUMs", disable=not args.verbose):
            outlines.append(f'>mum_{i}')
            cur = seq[mum_starts[i, args.index] : mum_starts[i, args.index] + mum_lengths[i]]  
            if mum_strands[i, args.index]:
                outlines.append(str(cur) + '#')
            else:
                outlines.append(str(cur.reverse_complement()) + '#')
        out.write('\n'.join(outlines))
if __name__ == "__main__":
    args = parse_arguments()
    main(args)

