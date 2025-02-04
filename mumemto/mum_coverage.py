#!/usr/bin/env python3

import numpy as np
from tqdm.auto import tqdm
import argparse
import os
import sys
from numba import njit
try:
    from utils import parse_mums_generator, get_sequence_lengths
except ImportError:
    from mumemto.utils import parse_mums_generator, get_sequence_lengths

def parse_arguments(args=None):    
    parser = argparse.ArgumentParser(description="Aggregates MUM coverage from mumemto output.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--input-prefix', '-i', dest='prefix', help='prefix for filelist, mums, and lengths files')
    group.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto')
    
    parser.add_argument('--lengths','-l', dest='lens', help='lengths file, first column is seq length in order of filelist')
    parser.add_argument('--len-filter','-L', dest='lenfilter', help='only consider MUMs longer than threshold', default=0, type=int)
    parser.add_argument('--verbose','-v', dest='verbose', help='verbose mode', action='store_true', default=False)
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)
    
    if args.mumfile:
        args.prefix = os.path.splitext(args.mumfile)[0]
    elif args.prefix:
        if args.prefix.endswith('.mums'):
            args.prefix = args.prefix[:-5]
        args.mumfile = args.prefix + '.mums'
    
    if args.lens is None:
        args.lens = args.prefix + '.lengths'
        
    return args

@njit
def update_coverage_single(coverage, start, length):
    """Update coverage array for a single MUM position"""
    coverage[start:start+length] = True

def update_coverage(coverage, mum_gen, verbose=False):
    """Update coverage array based on MUM positions"""
    for mum in tqdm(mum_gen, disable=not verbose, desc='Updating coverage'):
        length, starts, _ = mum
        for idx, start in enumerate(starts):
            if start != -1:  # Skip if MUM is not present in this sequence
                update_coverage_single(coverage[idx], start, length)
    return coverage

def main(args):
    # Read sequence lengths
    seq_lengths = get_sequence_lengths(args.lens)
    
    # Initialize coverage array
    max_len = max(seq_lengths)
    coverage = np.zeros((len(seq_lengths), max_len), dtype=bool)
    
    # Process MUMs using generator
    if args.verbose:
        print(f'Processing MUMs from {args.mumfile}...', file=sys.stderr)
    mum_gen = parse_mums_generator(args.mumfile, seq_lengths=seq_lengths, 
                                 lenfilter=args.lenfilter)
    
    # Update coverage
    coverage = update_coverage(coverage, mum_gen, verbose=args.verbose)
    
    # Print results
    print('coverages:')
    for i in range(len(seq_lengths)):
        print('seq%d: %.3f%%' % (i, np.count_nonzero(coverage[i, :seq_lengths[i]]) * 100 / seq_lengths[i]))

if __name__ == "__main__":
    args = parse_arguments()
    main(args)
