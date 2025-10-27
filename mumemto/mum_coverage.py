#!/usr/bin/env python3

import numpy as np
from tqdm.auto import tqdm
import argparse
import os
import sys
from numba import njit
try:
    from utils import parse_mums_generator, parse_bumbl_generator, get_sequence_lengths
except ImportError:
    from mumemto.utils import parse_mums_generator, parse_bumbl_generator, get_sequence_lengths

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
        elif args.prefix.endswith('.bumbl'):
            args.mumfile = args.prefix
            args.prefix = os.path.splitext(args.prefix)[0]
        elif args.prefix.endswith('.mums'):
            args.prefix = args.prefix[:-5]
            args.mumfile = args.prefix + '.mums'
        else:
            args.mumfile = args.prefix + '.mums'
    
    if args.lens is None:
        args.lens = args.prefix + '.lengths'
        
    return args

@njit
def update_coverage_ranges(coverage, starts, lengths):
    """Update coverage array for multiple MUM ranges using Numba acceleration"""
    n = len(starts)
    for i in range(n):
        start = starts[i]
        if start != -1:  # Skip if MUM is not present in this sequence
            length = lengths[i]
            coverage[start:start+length] = True

def update_coverage(coverage, mum_gen, lenfilter=0, verbose=False):
    """Update coverage array based on MUM positions"""
    for mum in tqdm(mum_gen, disable=not verbose, desc='Updating coverage'):
        length, start, _ = mum
        if start != -1 and length >= lenfilter: 
            coverage[start:start+length] = True
    return coverage

def update_coverage_chunked(coverage, mum_gen, lenfilter=0, verbose=False):
    """Update coverage array based on MUM positions from chunks"""
    for chunk in tqdm(mum_gen, disable=not verbose, desc='Updating coverage'):
        lengths, starts, _ = chunk
        # Filter out entries where MUM is not present and apply length filter
        valid_mask = (starts != -1) & (lengths >= lenfilter)
        if valid_mask.any():
            valid_starts = starts[valid_mask]
            valid_lengths = lengths[valid_mask]
            update_coverage_ranges(coverage, valid_starts, valid_lengths)
    return coverage

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
    
    # Determine file type and use appropriate generator
    if args.mumfile.endswith('.bumbl'):
        mum_gen = parse_bumbl_generator(args.mumfile, seq_idx=args.seq_idx, verbose=args.verbose, return_chunk=True)
        # Update coverage using chunked processing with numpy vectorization
        coverage = update_coverage_chunked(coverage, mum_gen, lenfilter=args.lenfilter, verbose=args.verbose)
    else:
        mum_gen = parse_mums_generator(args.mumfile, seq_idx=args.seq_idx, verbose=args.verbose)
        # Update coverage for individual MUMs
        coverage = update_coverage(coverage, mum_gen, lenfilter=args.lenfilter, verbose=args.verbose)
    
    # Print results
    coverage_pct = np.count_nonzero(coverage) * 100 / target_length
    print(f'seq{args.seq_idx}: {coverage_pct:.3f}%')

if __name__ == "__main__":
    args = parse_arguments()
    main(args)
