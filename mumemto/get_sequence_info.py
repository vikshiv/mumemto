#!/usr/bin/env python3

import argparse
import numpy as np
import os, sys

try:
    from utils import MUMdata, get_sequence_lengths
except ImportError:
    from mumemto.utils import MUMdata, get_sequence_lengths

def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description='Extract the MUM sequences')
    parser.add_argument('-m', '--mumfile', type=str, help='bumbl file to process', required=True)
    parser.add_argument('-o', '--output', type=str, help='Path to the output file')
    parser.add_argument('-l', '--lengths', type=str, help='Path to alternate lengths file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print verbose output')
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)
        
    if not args.mumfile.endswith('.mums') and not os.path.exists(args.mumfile + '.mums'):
        print(f"MUM file {args.mumfile} not found.", file=sys.stderr)
        sys.exit(1)
    elif not args.mumfile.endswith('.mums'):
        args.mumfile = args.mumfile + '.mums'
    if args.lengths is None:
        args.lengths = os.path.splitext(args.mumfile)[0] + '.lengths'
    if not os.path.exists(args.lengths):
        print(f"Lengths file {args.lengths} not found.", file=sys.stderr)
        sys.exit(1)
    if args.output is None:
        args.output = os.path.splitext(args.mumfile)[0] + '_labeled.mums'
    return args


def offset_mums(mums, lengths):
    print('Transforming MUMs to contig-relative coordinates...', file=sys.stderr)
    NUM_SEQS = len(lengths)
    ### offset the mums by the contig locations
    offsets = np.cumsum(lengths, axis=1)
    ### label each mum with the contig it belongs to
    # breaks = np.hstack((np.array([[0]*NUM_SEQS]).transpose(), offsets))
    contig_idx = np.array([np.searchsorted(offsets[idx], mums.starts[:,idx], side='right') for idx in range(NUM_SEQS)]).transpose()
    ### get the relative offset of each mum to the start of its contig
    left_start = np.hstack((np.zeros((offsets.shape[0],1), dtype=int), offsets[:,:-1]))
    rel_offsets = mums.starts - left_start[np.arange(NUM_SEQS), contig_idx]
    partial_mask = mums.starts != -1
    mums.starts[partial_mask] = rel_offsets[partial_mask]
    return contig_idx, rel_offsets

def main(args):
    try:
        lengths = get_sequence_lengths(args.lengths, multilengths=True)
    except ValueError as e:
        print("Multi-FASTA input required for contig ID annotation.", file=sys.stderr)
        sys.exit(1)
    mums = MUMdata(args.mumfile, sort=False, verbose=args.verbose)
    mum_lengths, mum_strands = mums.lengths, mums.strands
    
    contig_idx, rel_offsets = offset_mums(mums, lengths)
    mums.starts = rel_offsets
    with open(args.output, 'w') as out:
        for i in range(mums.num_mums):
            strands_str = ['+' if s else '-' for s in mum_strands[i]]
            cur_contig_idx = ','.join(map(str, contig_idx[i]))
            out.write(f"{mum_lengths[i]}\t{','.join(map(str, rel_offsets[i]))}\t{','.join(strands_str)}\t{cur_contig_idx}\n")
    
if __name__ == "__main__":
    args = parse_arguments()
    main(args)

