#!/usr/bin/env python3

import argparse
import numpy as np
import os, sys

try:
    from utils import MUMdata, get_sequence_lengths, serialize_coll_blocks, get_contig_names
except ImportError:
    from mumemto.utils import MUMdata, get_sequence_lengths, serialize_coll_blocks, get_contig_names

def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description='Extract the MUM sequences')
    parser.add_argument('-m', '--mumfile', type=str, help='bumbl file to process', required=True)
    parser.add_argument('-o', '--output', type=str, help='Path to the output file')
    parser.add_argument('-l', '--lengths', type=str, help='Path to alternate lengths file')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print verbose output')
    parser.add_argument('-n', '--contig-names', dest='contig_names', action='store_true', help='Print contig names/sequence IDs instead of numerical indexes')
    
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
        args.output = os.path.splitext(args.mumfile)[0] + '_labeled.mums'
    return args


def offset_mums(args, mums, lengths):
    if args.verbose:
        print('Transforming MUMs to contig-relative coordinates...', file=sys.stderr)
    NUM_SEQS = len(lengths)
    contig_idx = []
    rel_offsets = []
    for i in range(NUM_SEQS):
        ### offset the mums by the contig locations
        offsets = np.cumsum(lengths[i])
        ### label each mum with the contig it belongs to
        contig_idx.append(np.searchsorted(offsets, mums.starts[:,i], side='right'))
        ### get the relative offset of each mum to the start of its contig
        left_start = np.hstack((0, offsets[:-1]))
        rel_offsets.append(mums.starts[:,i] - left_start[contig_idx[i]])
    contig_idx = np.vstack(contig_idx).transpose()
    rel_offsets = np.vstack(rel_offsets).transpose()
    partial_mask = mums.starts == -1
    rel_offsets[partial_mask] = -1
    return contig_idx, rel_offsets

def main(args):
    try:
        lengths = get_sequence_lengths(args.lengths, multilengths=True)
    except ValueError as e:
        print("Multi-FASTA input required for contig ID annotation.", file=sys.stderr)
        sys.exit(1)
    if args.contig_names:
        names = get_contig_names(args.lengths)
    
    mums = MUMdata(args.mumfile, sort=False, verbose=args.verbose)
    mum_lengths, mum_starts, mum_strands = mums.lengths, mums.starts, mums.strands
    is_blocked = mums.blocks is not None
    if is_blocked:
        blocks = serialize_coll_blocks(mums.blocks, mums.num_mums)
        
    contig_idx, rel_offsets = offset_mums(args, mums, lengths)
    mums.starts = rel_offsets
    with open(args.output, 'w') as out:
        for i in range(mums.num_mums):
            strands_str = ['+' if s else '-' for s in mum_strands[i]]
            if args.contig_names:
                cur_contig_idx = ','.join([names[idx][i] for idx, i in enumerate(contig_idx[i])])
            else:
                cur_contig_idx = ','.join(map(str, contig_idx[i]))
            if is_blocked:
                out.write(f"{mum_lengths[i]}\t{','.join(map(str, mum_starts[i]))}\t{','.join(strands_str)}\t{blocks[i]}\t{cur_contig_idx}\t{','.join(map(str, rel_offsets[i]))}\n")
            else:
                out.write(f"{mum_lengths[i]}\t{','.join(map(str, mum_starts[i]))}\t{','.join(strands_str)}\t*\t{cur_contig_idx}\t{','.join(map(str, rel_offsets[i]))}\n")
    
if __name__ == "__main__":
    args = parse_arguments()
    main(args)

