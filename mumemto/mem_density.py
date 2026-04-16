import numpy as np
from tqdm.auto import tqdm
from numba import njit
import argparse
import os
try:
    from utils import get_sequence_lengths
except ImportError:
    from mumemto.utils import get_sequence_lengths

def parse_arguments():    
    parser = argparse.ArgumentParser(description="Aggregates MEM density for downstream plotting")
    parser.add_argument('--mems', '-m', dest='mems', help='path to mems from mumemto', required=True)
    parser.add_argument('--lengths', '-l', dest='lengths', help='lengths of each sequence', required=True)
    parser.add_argument('--verbose', '-v', dest='verbose', action='store_true', help='enable verbose mode')
    args = parser.parse_args()
    return args


@njit
def update_coverage(coverage, idxs, starts, l):
    for start, idx in zip(starts, idxs):
        coverage[idx, start:start+l] += 1

def main(args):
    file = open(args.mems, 'r')
    args.lengths = get_sequence_lengths(args.lengths)
    args.size = max(args.lengths)
    args.num = len(args.lengths)
    coverage = np.zeros((args.num, args.size))
    for m in tqdm(file, disable=not args.verbose, desc='Parsing '+args.mems):
        m = m.strip().split()
        l = int(m[0])
        idxs = np.fromstring(m[2], sep=',', dtype=int)
        starts = np.fromstring(m[1], sep=',', dtype=int)
        update_coverage(coverage, idxs, starts, l)
    filename = os.path.splitext(args.mems)[0]
    coverage_file = filename + '_coverage.npy'
    np.save(coverage_file, coverage)

if __name__ == "__main__":
    args = parse_arguments()
    main(args)