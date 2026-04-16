#!/usr/bin/env python3

import argparse
import numpy as np
import os
import sys

try:
    from utils import MUMdata, get_sequence_lengths
except ImportError:
    from mumemto.utils import MUMdata, get_sequence_lengths

def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description='Trim MUMs that extend beyond sequence boundaries')
    parser.add_argument('-m', '--mumfile', type=str, help='MUM or bumbl file to process', required=True)
    parser.add_argument('-l', '--lengths', type=str, help='Path to lengths file')
    parser.add_argument('-o', '--output', type=str, help='Path to output file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print verbose output')
    parser.add_argument('--min-length', type=int, default=20, help='Minimum MUM length after trimming (default: 20)')
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)
    
    if not args.mumfile.endswith('.mums') and not args.mumfile.endswith('.bumbl'):
        if not os.path.exists(args.mumfile + '.mums') and not os.path.exists(args.mumfile + '.bumbl'):
            print(f"MUM file {args.mumfile} not found.", file=sys.stderr)
            sys.exit(1)
        elif os.path.exists(args.mumfile + '.mums'):
            args.mumfile = args.mumfile + '.mums'
        elif os.path.exists(args.mumfile + '.bumbl'):
            args.mumfile = args.mumfile + '.bumbl'
    elif not os.path.exists(args.mumfile):
        print(f"MUM file {args.mumfile} not found.", file=sys.stderr)
        sys.exit(1)
    
    if args.lengths is None:
        args.lengths = os.path.splitext(args.mumfile)[0] + '.lengths'
    if not os.path.exists(args.lengths):
        print(f"Lengths file {args.lengths} not found.", file=sys.stderr)
        sys.exit(1)
    
    if args.output is None:
        base = os.path.splitext(args.mumfile)[0]
        ext = '.mums' if args.mumfile.endswith('.mums') else '.bumbl'
        args.output = base + '_trimmed' + ext
    
    return args

def trim_mums(mums, seq_lengths, min_length=20, verbose=False):
    """
    Trim MUMs that extend beyond sequence boundaries using vectorized numpy operations.
    """    
    seq_lengths_arr = np.array(seq_lengths, dtype=np.int64)
    valid_mask = mums.starts != -1
    lengths_2d = mums.lengths[:, np.newaxis]
    
    excess = mums.starts + lengths_2d - seq_lengths_arr
    excess[~valid_mask] = np.iinfo(np.int64).min
    max_excess = np.max(excess, axis=1)
    needs_trimming = max_excess > 0
    
    new_lengths = mums.lengths.copy()
    new_lengths[needs_trimming] = mums.lengths[needs_trimming] - max_excess[needs_trimming]
    
    trimmed_count = np.sum(needs_trimming)
    
    below_min = new_lengths < min_length
    removed_count = np.sum(below_min)
    new_lengths[below_min] = 0
    valid_mask = new_lengths > 0
    if not np.all(valid_mask):
        new_lengths = new_lengths[valid_mask]
        new_starts = mums.starts[valid_mask]
        new_strands = mums.strands[valid_mask]
    else:
        new_starts = mums.starts
        new_strands = mums.strands
    
    if verbose:
        print(f"Trimmed {trimmed_count} MUMs", file=sys.stderr)
        if removed_count > 0:
            print(f"Removed {removed_count} MUMs that would be below minimum length ({min_length}bp)", file=sys.stderr)
    
    return MUMdata.from_arrays(new_lengths, new_starts, new_strands)

def main(args):
    seq_lengths = get_sequence_lengths(args.lengths)
    mums = MUMdata(args.mumfile, sort=False, verbose=args.verbose)
    if args.verbose:
        print(f"Loaded {mums.num_mums} MUMs across {mums.num_seqs} sequences", file=sys.stderr)
        print(f"Sequence lengths: {seq_lengths}", file=sys.stderr)
    
    trimmed_mums = trim_mums(mums, seq_lengths, min_length=args.min_length, verbose=args.verbose)
    
    if args.verbose:
        print(f"Output: {trimmed_mums.num_mums} MUMs", file=sys.stderr)
    
    if args.output.endswith('.bumbl'):
        trimmed_mums.write_bums(args.output)
    else:
        trimmed_mums.write_mums(args.output, blocks=trimmed_mums.blocks)
    
    if args.verbose:
        print(f"Written to {args.output}", file=sys.stderr)

if __name__ == "__main__":
    args = parse_arguments()
    main(args)

