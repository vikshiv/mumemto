#!/usr/bin/env python3

import os
import sys
import argparse
import numpy as np
from mumemto.utils import MUMdata, find_coll_blocks, get_block_order
from tqdm import tqdm

def parse_arguments():
    parser = argparse.ArgumentParser(description="Detect inversions from MUMs. Optionally checks if inversions are flanked by scaffold breaks when AGP files are provided.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--input-prefix', '-i', dest='prefix', help='prefix for filelist, mums, and lengths files')
    group.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto')
    parser.add_argument('--lengths', '-l', dest='lens', help='Path to lengths file')
    
    parser.add_argument('--agp-filelist', '-a', dest='agp_filelist', help='Path to filelist containing AGP files. Each line should contain the path to an AGP file, in the same order as sequences in the mumemto filelist/lengths file. Assume first sequence is reference.')
    parser.add_argument('--filelist', '-f', dest='filelist', help='Path to filelist for sequence names')
    parser.add_argument('--chr', '-c', help='Chromosome number (required if using --agp-filelist)')
    parser.add_argument('--margin', '-d', dest='margin', type=float, default=0.01, help='Proximity margin for scaffold break detection (default: within 1%% of inversion length)')
    parser.add_argument('--max-length', '-L', dest='max_length', type=int, help='Maximum inversion length to detect')
    parser.add_argument('--max-block-gap-len','-g', dest='max_block_gap', help='maximum break between collinear mums within a collinear block (default: 1000)', default=1000, type=int)

    parser.add_argument('--verbose', '-v', action='store_true', help='Print progress updates')

    args = parser.parse_args()
    
    # Validate AGP and chromosome arguments
    if bool(args.agp_filelist) ^ bool(args.chr):
        parser.error("--agp-filelist and --chr must be provided together")
    args.scaffold = bool(args.agp_filelist) and bool(args.chr)
    
    if args.mumfile:
        args.prefix = os.path.splitext(args.mumfile)[0]
    elif args.prefix:
        if args.prefix.endswith('.mums'):
            args.prefix = args.prefix[:-5]
        args.mumfile = args.prefix + '.mums'
    else:
        parser.error("Either --mums or --prefix must be provided")
    
    if args.lens is None:
        args.lens = args.prefix + '.lengths'
        
    if not args.filelist:
        args.filelist = args.prefix + '_filelist.txt'
        
    return args

def get_sequence_info(args):
    """Load sequence lengths and names"""
    seq_lengths = np.array([int(l.split()[1]) for l in open(args.lens).readlines()])
    
    if args.filelist:
        if args.chr:
            hap_ids = [os.path.basename(l.split()[0]).split(f'_chr{args.chr}')[0] 
                      for l in open(args.filelist).readlines()]
        else:
            hap_ids = [os.path.basename(l.split()[0]) for l in open(args.filelist).readlines()]
    else:
        hap_ids = [f'seq_{i}' for i in range(len(seq_lengths))]
        
    return seq_lengths, hap_ids

def get_scaffold_breaks(args, hap_ids):
    """Load contig break positions from AGP files"""
    breaks = {}
    contig_names = {}
        
    # Read AGP filelist and create mapping of hap_id to AGP file
    agp_files = {}
    with open(args.agp_filelist) as f:
        for idx, line in enumerate(f):
            if idx >= len(hap_ids):  # Skip if we have more AGP files than sequences
                break
            agp_path = line.strip()
            if agp_path:  # Skip empty lines
                agp_files[hap_ids[idx+1]] = agp_path
    
    # Process AGP files with progress bar if verbose
    haps_iter = tqdm(hap_ids[1:], desc="Processing AGP files", disable=not args.verbose)
    for hap in haps_iter:
        if hap not in agp_files:
            continue
        agp_file = agp_files[hap]
        if not os.path.exists(agp_file):
            if args.verbose:
                print(f"Warning: AGP file not found: {agp_file}", file=sys.stderr)
            continue
            
        with open(agp_file) as f:
            lines = [l for l in f.read().splitlines() if l.startswith('chr' + str(args.chr))]
            lengths = [l.split() for l in lines if l.split()[4] == 'W']
            breaks[hap] = np.array([int(l[2]) - int(l[1]) + 1 for l in lengths])
            contig_names[hap] = [l[5] for l in lengths if l[4] == 'W']
    return breaks, contig_names

def find_reversals(coll_block_order, mums, blocks):
    # output ranges are inclusive coll_block_order[s:e+1] ranges
    stretches = []
    for i in range(1, len(coll_block_order)):
        decreases = np.where(np.diff(coll_block_order[i]) == -1)[0]
        ranges = np.split(decreases, np.where(np.diff(decreases) != 1)[0] + 1)
        for r in ranges:
            if len(r) == 0:
                continue
            if np.all([not mums[blocks[x][0]][2][i] for x in coll_block_order[i][r[0]:r[-1]+2]]):
                stretches.append((i, r[0], r[-1]+1))
    return stretches

def inversion_coords(coll_block_order, mums, blocks, i, s, e):
    block_range = coll_block_order[i][s:e+1]
    first, last = block_range[0], block_range[-1]
    # Get coordinates in sequence i
    seq_start = mums[blocks[first][1]][1][i]
    seq_end = mums[blocks[last][0]][1][i] + mums[blocks[last][0]][0]
    # Get coordinates in reference (sequence 0)
    ref_start = mums[blocks[first][1]][1][0]
    ref_end = mums[blocks[last][0]][1][0] + mums[blocks[last][0]][0]
    return (i, seq_start, seq_end, ref_start, ref_end)

def main():
    args = parse_arguments()
    
    if args.verbose:
        print("hahaLoading sequence information...", file=sys.stderr)
    
    # Load data
    seq_lengths, hap_ids = get_sequence_info(args)
    
    # Only get scaffold breaks if both AGP and chr are provided
    if args.scaffold:
        breaks, contig_names = get_scaffold_breaks(args, hap_ids)
        
    # Load and process MUMs
    mums = MUMdata(args.mumfile, seq_lengths=seq_lengths, verbose=args.verbose)
    
    # Find collinear blocks and inversions
    _, small_blocks, _ = find_coll_blocks(mums, max_break=args.max_block_gap, verbose=args.verbose)
    block_orders = get_block_order(mums, small_blocks)
    if args.verbose:
        print("Finding inversions...", file=sys.stderr)
    
    # Find inversions using find_reversals and inversion_coords
    reversed_stretches = find_reversals(block_orders, mums, small_blocks)
    reversed_ranges = []
    for i, s, e in reversed_stretches:
        ranges = inversion_coords(block_orders, mums, small_blocks, i, s, e)
        if args.max_length is None or np.abs(ranges[2] - ranges[1]) <= args.max_length:
            reversed_ranges.append(ranges)
    
    if args.verbose:
        print(f"Found {len(reversed_ranges)} inversions", file=sys.stderr)
        print("Writing results...", file=sys.stderr)
    
    # Output results to stdout
    print(f"hap_id\tstart\tend\tref_start\tref_end" + ("\tscaffold_break\tcontig" if args.scaffold else ""))
    for seq_idx, start, end, ref_start, ref_end in reversed_ranges:
        hap = hap_ids[seq_idx]
        
        # Only check scaffold breaks if AGP files were provided
        if args.scaffold and hap in breaks:
            diffs_start = np.abs(np.cumsum(breaks[hap]) - start)
            diffs_end = np.abs(np.cumsum(breaks[hap]) - end)
            margin = (end - start) * args.margin  # within X% of inversion size
            contig_id = []
            if diffs_start.min() < margin:
                contig_id.extend([contig_names[hap][x] for x in np.where(diffs_start < margin)[0]])
            if diffs_end.min() < margin:
                contig_id.extend([contig_names[hap][x] for x in np.where(diffs_end < margin)[0]])
            # Output with scaffold break info
            print(f"{hap}\t{start}\t{end}\t{ref_start}\t{ref_end}\t{True if contig_id else False}\t{','.join(contig_id) if contig_id else 'NA'}")
        else:
            # Basic output without scaffold break info
            print(f"{hap}\t{start}\t{end}\t{ref_start}\t{ref_end}")

if __name__ == '__main__':
    main()

