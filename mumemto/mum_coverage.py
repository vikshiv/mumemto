#!/usr/bin/env python3

import numpy as np
from tqdm.auto import tqdm
import argparse
import os
import sys
from numba import njit
try:
    from utils import stream_mums, get_sequence_lengths
except ImportError:
    from mumemto.utils import stream_mums, get_sequence_lengths

def parse_arguments(args=None):    
    parser = argparse.ArgumentParser(description="Aggregates MUM coverage from mumemto output.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--input-prefix', '-i', dest='prefix', help='prefix for filelist, mums, and lengths files')
    group.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto')
    
    parser.add_argument('--lengths','-l', dest='lens', help='lengths file, first column is seq length in order of filelist')
    parser.add_argument('--len-filter','-L', dest='lenfilter', help='only consider MUMs longer than threshold', default=0, type=int)
    parser.add_argument('--seq-idx', '-s', dest='seq_idx', help='sequence index to compute coverage for', default=0, type=int)
    parser.add_argument('--verbose','-v', dest='verbose', help='verbose mode', action='store_true', default=False)
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)
    
    if args.mumfile:
        args.prefix = os.path.splitext(args.mumfile)[0]
    elif args.prefix:
        # Check if both .bumbl and .mums exist
        bumbl_exists = os.path.exists(args.prefix + '.bumbl')
        mums_exists = os.path.exists(args.prefix + '.mums')
        
        if bumbl_exists and mums_exists:
            print(f"Error: Both {args.prefix + '.bumbl'} and {args.prefix + '.mums'} exist. Please specify the file explicitly with --mums.", file=sys.stderr)
            sys.exit(1)
        elif bumbl_exists:
            args.mumfile = args.prefix + '.bumbl'
        elif mums_exists:
            args.mumfile = args.prefix + '.mums'
        elif args.prefix.endswith('.bumbl') or args.prefix.endswith('.mums'):
            args.mumfile = args.prefix
            args.prefix = os.path.splitext(args.prefix)[0]
        else:
            args.mumfile = args.prefix + '.mums'
    
    if args.lens is None:
        args.lens = args.prefix + '.lengths'
        
    return args


@njit
def update_coverage(coverage, start, length, lenfilter=0):
    if start != -1 and length >= lenfilter: 
        coverage[start:start+length] = True

def main(args):
    # Read sequence lengths
    seq_lengths = get_sequence_lengths(args.lens)
    
    # Get target sequence length
    if args.seq_idx >= len(seq_lengths) or args.seq_idx < 0:
        print(f'Error: sequence index {args.seq_idx} is out of range (0-{len(seq_lengths)-1})', file=sys.stderr)
        sys.exit(1)
    
    target_length = seq_lengths[args.seq_idx]
    
    # Initialize coverage array for single sequence
    coverage = np.zeros(target_length, dtype=bool)
    
    # compute coverage over mums
    mum_gen = stream_mums(args.mumfile, seq_idx=args.seq_idx, verbose=args.verbose)
    for mum in mum_gen:
        length, start, _ = mum
        update_coverage(coverage, start, length, args.lenfilter)
        
    # Print results
    coverage_pct = np.count_nonzero(coverage) * 100 / target_length
    print(f'seq{args.seq_idx}: {coverage_pct:.3f}%', file=sys.stderr)

if __name__ == "__main__":
    args = parse_arguments()
    main(args)
