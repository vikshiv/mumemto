import numpy as np
from tqdm.auto import tqdm
from collections import namedtuple

MUM = namedtuple('MUM', ['length', 'starts', 'strands'])

def find_coll_blocks(mums, max_break=100000, verbose=False):
    starts = mums.starts
    strands = mums.strands
    lengths = mums.lengths
    mum_orders = starts.transpose().argsort()
    mum_gaps = []
    flips = set([])
    for i in tqdm(range(mum_orders.shape[0]), desc='Finding collinear blocks', disable=not verbose):
        cur = []
        for l in range(1, mum_orders.shape[1]):
            left, right = mum_orders[i][l-1], mum_orders[i][l]
            if strands[left][i] == strands[right][i]:
                if strands[left][i]:
                    cur.append((left, right))
                else:
                    cur.append((right, left))
                    flips.add((right, left))
        mum_gaps.append(cur)
    common_gaps = set.intersection(*map(set, mum_gaps))
    left, right = zip(*common_gaps)
    left, right = set(list(left)), set(list(right))
    true_collinear_mums = sorted(list(left.intersection(right)))
    right_coll_mums = sorted(list(left.difference(set(true_collinear_mums)))) # have a right pair, but not a left
    left_coll_mums = sorted(list(right.difference(set(true_collinear_mums)))) # have a left pair, but not a right
    large_blocks = list(zip(right_coll_mums, left_coll_mums))
    ### find the longest stretches of collinear mums
    small_blocks = []
    for l, r in large_blocks:
        last = l
        for i in range(l, r):
            lens = np.full(len(starts[i]), lengths[i])
            lens[(starts[i+1] < starts[i])] = lengths[i+1] 
            gap_lens = np.abs(starts[i] - starts[i+1]) - lens
            if gap_lens.max() > max_break and last < i:
                small_blocks.append((last, i))
                last = i + 1
        if last != r:
            small_blocks.append((last, r))
    return large_blocks, small_blocks, mum_gaps

def get_block_order(mums, blocks):
    starts = mums.starts
    mum_orders = starts.transpose().argsort()
    ### get coll_block order
    coll_block_starts = [b[0] for b in blocks]
    coll_block_orders = []
    for i in range(mum_orders.shape[0]):
        poi = np.where(np.isin(mum_orders[i], coll_block_starts))[0]
        values = mum_orders[i][poi]
        poi = np.argsort(values)
        coll_block_orders.append(np.argsort(poi))
    return coll_block_orders

def parse_mums_generator(mumfile, seq_lengths, lenfilter=0, subsample=1):
    """Generator that streams MUMs from mumfile"""
    count = 0
    with open(mumfile, 'r') as f:
        for line in f:
            if subsample == 1 or count % subsample == 0:
                line = line.strip().split()
                length = int(line[0])
                if length >= lenfilter:
                    # Parse the line
                    strands = [s == '+' for s in line[2].split(',')]
                    starts = [int(pos) if pos != '' else -1 for pos in line[1].split(',')]
                    yield MUM(length, starts, strands)
            count += 1

class MUMdata:
    def __init__(self, mumfile, seq_lengths=None, lenfilter=0, subsample=1, verbose=False):
        # If seq_lengths not provided, try to load from .lengths file
        if seq_lengths is None:
            lengths_file = mumfile.replace('.mums', '.lengths')
            try:
                with open(lengths_file, 'r') as f:
                    seq_lengths = [int(l.strip()) for l in f.readlines()]
            except FileNotFoundError:
                raise ValueError("Either a *.lengths file or an input seq_lengths array is required")
                
        self.lengths, self.starts, self.strands = self.parse_mums(
            mumfile, 
            seq_lengths, 
            lenfilter, 
            subsample,
            verbose
        )
        self.num_mums = len(self.lengths)
        self.num_seqs = self.starts.shape[1] if self.num_mums > 0 else 0
        # sort by reference offset position
        order = self.starts[:,0].argsort()
        self.lengths = self.lengths[order]
        self.starts = self.starts[order]
        self.strands = self.strands[order]
        
    @staticmethod
    def parse_mums(mumfile, seq_lengths, lenfilter=0, subsample=1, verbose=False):
        lengths, starts, strands = [], [], []
        
        count = 0
        for l in tqdm(open(mumfile, 'r').readlines(), disable=not verbose, desc='Parsing MUMs'):
            if subsample == 1 or count % subsample == 0:
                l = l.strip().split()
                length = int(l[0])
                if length >= lenfilter:
                    # Parse the line
                    start_positions = [int(v) if v else -1 for v in l[1].split(',')]
                    strand_info = l[2].split(',')
                    
                    # Handle reverse strands
                    # for idx, (pos, strand) in enumerate(zip(start_positions, strand_info)):
                    #     if strand == '-':
                    #         start_positions[idx] = seq_lengths[idx] - pos - length
                    
                    lengths.append(length)
                    starts.append(start_positions)
                    strands.append([s == '+' for s in strand_info])
            count += 1
        
        # Convert to numpy arrays all at once
        return (
            np.array(lengths),
            np.array(starts),
            np.array(strands)
        )

    def filter_pmums(self):
        """Remove any MUMs that have -1 in their start positions"""
        valid_rows = ~np.any(self.starts == -1, axis=1)
        self.lengths = self.lengths[valid_rows]
        self.starts = self.starts[valid_rows]
        self.strands = self.strands[valid_rows]
        self.num_mums = len(self.lengths)
        return self

    def __iter__(self):
        """Iterate over MUMs, yielding (length, starts, strands) for each"""
        for i in range(self.num_mums):
            yield self[i]

    def __getitem__(self, idx):
        """Get a single MUM as (length, starts, strands)"""
        return MUM(self.lengths[idx], self.starts[idx], self.strands[idx])

    def __len__(self):
        """Return number of MUMs"""
        return self.num_mums