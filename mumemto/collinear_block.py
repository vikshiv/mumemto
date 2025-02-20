#!/usr/bin/env python3

import os, sys
import numpy as np
import argparse
try:
    from utils import MUMdata, find_coll_blocks, get_sequence_lengths
except ImportError:
    from mumemto.utils import MUMdata, find_coll_blocks, get_sequence_lengths

def parse_arguments(args=None):    
    parser = argparse.ArgumentParser(description="Plots a synteny plot of MUMs from mumemto")
    # parser.add_argument('--filelist', '-f', dest='filelist', help='path to filelist from mumemto', required=True)
    # parser.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto', required=True)
    # parser.add_argument('--lengths','-l', dest='lens', help='lengths file, first column is seq length in order of filelist', required=True)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--input-prefix', '-i', dest='prefix', help='prefix for filelist, mums, and lengths files')
    group.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto')
    
    parser.add_argument('--lengths','-l', dest='lens', help='lengths file, first column is seq length in order of filelist')
    
    parser.add_argument('--fout','-o', dest='filename', help='plot fname (default: input_prefix + _sorted)')
    parser.add_argument('--max-gap-len','-g', dest='max_break', help='maximum break between collinear mums within a collinear block (default: 1kbp)', default=1000, type=int)
    parser.add_argument('--verbose','-v', dest='verbose', help='verbose mode', action='store_true', default=False)
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)
    
    if args.mumfile:
        args.prefix = os.path.splitext(args.mumfile)[0]
    elif args.prefix:
        if args.prefix.endswith('.mums') or args.prefix.endswith('.bums'):
            args.prefix = args.prefix[:-5]
        args.mumfile = args.prefix + '.mums'
    else:
        parser.error("Either --mums or --prefix must be provided")
        
    if args.lens is None:
        args.lens = args.prefix + '.lengths'
    
    if args.filename is None:
        args.filename = args.prefix + '_sorted'
    
    return args

def main(args):
    seq_lengths = get_sequence_lengths(args.lens)
    max_length = max(seq_lengths)
    
    mums = MUMdata(args.mumfile, verbose=args.verbose)
    if args.verbose:
        print(f'Found {mums.num_mums} MUMs', file=sys.stderr)
            
    ### filter out pmums for collinear blocks
    mums.filter_pmums()
    if len(mums) == 0:
        print('No strict MUMs found after filtering partial MUMs.', file=sys.stderr)
        return
    collinear_blocks = find_coll_blocks(mums, max_break=args.max_break, verbose=args.verbose)
    if args.verbose:
        print(f'found {len(collinear_blocks)} collinear blocks', file=sys.stderr)
    
    mums.write_bums(args.filename)
    with open(args.filename + '.blocks', 'w') as f:
        for l, r in collinear_blocks:
            f.write(f'{l}\t{r}\n')
    
if __name__ == "__main__":
    args = parse_arguments()
    main(args)
