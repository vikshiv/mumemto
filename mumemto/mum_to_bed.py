import numpy as np
import os
import argparse
from mumemto.utils import get_sequence_lengths, stream_mums
import sys


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(description='Convert MUMs file to BED file')
    parser.add_argument('mums_file', help='Path to the MUMs file')
    parser.add_argument('--lengths-file', '-l', help='Path to the lengths file (optional)')
    parser.add_argument('-v', '--verbose', action='store_true', 
                       help='Enable verbose output with progress bars')
    parser.add_argument('--min-singleton-length', '-L', type=int, default=100,
                        help='Minimum length of singleton blocks to include')
    parser.add_argument('--seq-idx', '-s', type=int, default=0,
                        help='Sequence to output BED coordinates (default: first sequence)')
    parser.add_argument('--output', '-o', help='Path to the output file (default: stdout)', default=None)
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    if args.lengths_file is None:
        args.lengths_file = get_lengths_file(args.mums_file)
    
    return args

def get_lengths_file(mums_file):
    """Get the corresponding lengths file path"""
    base = os.path.splitext(mums_file)[0]
    lengths_file = f"{base}.lengths"
    if not os.path.exists(lengths_file):
        raise FileNotFoundError(f"Lengths file {lengths_file} not found")
    return lengths_file
    
def process_mums_file(mums_file, seq_idx=0, verbose=False, min_singleton_length=100):
    """Process MUMs file and calculate statistics"""
    with open(mums_file, 'r') as f:
        line = f.readline().split()
        if (len(line) < 4 or line[3] == '*'):
            print('MUMs file does not contain blocks. Try running mumemto collinear first.', file=sys.stderr)
            sys.exit(1)
    last_block = '-'
    last_start, last_end = None, None
    last_strand = None
    intervals = []
    mum_idx = 0
    has_blocks = False
    for l, start, strand, block in stream_mums(mums_file, seq_idx=seq_idx, verbose=verbose, return_blocks=True):
        if block is None:
            print('No collinear blocks found. Only writing mums to BED intervals.', file=sys.stderr)
        else:
            has_blocks = True
        if has_blocks:
            # Handle block intervals
            if last_block != '-' and block == last_block:  # n,n
                if strand:
                    last_end = start + l
                else:
                    last_start = start
            elif last_block != '-' and block != last_block:  # n,- or n,n+1
                intervals.append((last_start, last_end, last_strand, f'block_{last_block}'))
                if block != '-':  # n,n+1
                    if strand:
                        last_start = start
                    else:
                        last_end = start + l
            elif block != '-':  # -,n
                if strand:
                    last_start = start
                else:
                    last_end = start + l
            if block == '-' and l >= min_singleton_length:
                intervals.append((start, start+l, strand, f'mum_{mum_idx}'))
        elif l >= min_singleton_length: ## no blocks
            intervals.append((start, start+l, strand, f'mum_{mum_idx}'))
            
        last_block = block
        last_strand = strand
        mum_idx += 1
    return intervals

def find_chr(intervals, lengths):
    offsets = np.cumsum(lengths)
    starts = np.array([i[0] for i in intervals])
    ### label each mum with the contig it belongs to
    contig_idx = np.searchsorted(offsets, starts, side='right')
    ### get the relative offset of each mum to the start of its contig
    left_start = np.hstack((0, offsets[:-1]))
    rel_offsets = starts - left_start[contig_idx]
    return contig_idx, rel_offsets


def get_contig_names(lengths_file):
    ### assumes lengths_file is formatted as multilengths
    names = []
    cur_name = []
    first_line = True
    for l in open(lengths_file, 'r').readlines():
        l = l.strip().split()
        if first_line and l[1] != '*':
            print('Lengths file must be formatted as multilengths.', file=sys.stderr)
            sys.exit(1)
        first_line = False
        if l[1] == '*':
            if cur_name:
                names.append(cur_name)
            cur_name = []
            continue
        cur_name.append(l[1])
    names.append(cur_name)
    return names

def main(args):
    lengths = get_sequence_lengths(args.lengths_file, multilengths=True)
    if args.seq_idx >= len(lengths):
        print(f"Sequence index {args.seq_idx} too large for dataset with {len(lengths)} sequences.", file=sys.stderr)
        sys.exit(1)
    lengths = lengths[args.seq_idx]
    contig_names = get_contig_names(args.lengths_file)[args.seq_idx]

    intervals = process_mums_file(args.mums_file, args.seq_idx, args.verbose, args.min_singleton_length)
    contig_idx, rel_offsets = find_chr(intervals, lengths)

    if args.output is None:
        outfile = sys.stdout
    else:
        outfile = open(args.output, 'w')
    for i in range(len(intervals)):
        interval_len = intervals[i][1] - intervals[i][0]
        outfile.write(f'{contig_names[contig_idx[i]]}\t{rel_offsets[i]}\t{rel_offsets[i] + interval_len}\t{intervals[i][3]}\t{'+' if intervals[i][2] else '-'}\n')
    
    outfile.close()
    
    


if __name__ == "__main__":
    args = parse_arguments()
    main(args)
