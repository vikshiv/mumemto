import numpy as np
from tqdm.auto import tqdm
from collections import namedtuple

MUM = namedtuple('MUM', ['length', 'starts', 'strands'])

def find_coll_blocks(mums, max_break=0, verbose=False, return_order=False):
    def find_blocks(coll_mums):
        diffs = np.diff(np.concatenate(([False], coll_mums, [False])).astype(int))
        starts = np.where(diffs == 1)[0]
        ends = np.where(diffs == -1)[0]
        blocks = list(zip(starts, ends))
        return blocks
    """
    The output bool array (quick_coll) is True in position i if mum[i] and mum[i+1] are collinear
    This function works by finding the ranks of each mum, then identifying consecutive MUMs with consecutive ranks.
    The strand orientiation of collinear MUMs must be identical, and MUMs in - strands should be reversed in rank.
    """
    starts = mums.starts
    strands = mums.strands
    lengths = mums.lengths
    mum_orders = starts.transpose().argsort()
    strand_changes = (~np.diff(strands, axis=0)).all(axis=1)
    strand_change_diff = np.where(strands, 1, -1)
    mum_order_pos = np.argsort(mum_orders, axis=1)
    quick_coll = (strand_change_diff.T[:, :-1] == np.diff(mum_order_pos, axis=1)).all(axis=0)
    np.logical_and(quick_coll, strand_changes, out = quick_coll)
    large_blocks = find_blocks(quick_coll)
    if max_break > 0:
        small_blocks = []
        for l, r in tqdm(large_blocks, desc='Truncating blocks (max gap length > {})'.format(max_break), disable=not verbose):
            last = l
            for i in range(l, r):
                lens = np.full(len(starts[i]), lengths[i])
                lens[(starts[i+1] < starts[i])] = lengths[i+1] 
                gap_lens = np.abs(starts[i] - starts[i+1]) - lens
                if (gap_lens.max() > max_break):
                    if last < i:
                        small_blocks.append((last, i))
                    last = i + 1
            if last != r:
                small_blocks.append((last, r))
        blocks = small_blocks
    else:
        blocks = large_blocks
    if return_order:
        order = mum_order_pos[:,[b[0] for b in blocks]].argsort(axis=1)
        return blocks, order
    return blocks

def parse_mums_generator(mumfile, lenfilter=0, subsample=1, verbose=False):
    """Generator that streams MUMs from mumfile"""
    count = 0
    with open(mumfile, 'r') as f:
        for line in tqdm(f, desc='parsing MUM file', disable=not verbose):
            if subsample == 1 or count % subsample == 0:
                line = line.strip().split()
                length = int(line[0])
                if length >= lenfilter:
                    # Parse the line
                    strands = [s == '+' for s in line[2].split(',')]
                    starts = [int(pos) if pos != '' else -1 for pos in line[1].split(',')]
                    yield MUM(length, starts, strands)
            count += 1

def get_sequence_lengths(lengths_file, multilengths=False):
    def get_lengths(lengths_file):
        return [int(l.split()[1]) for l in open(lengths_file, 'r').read().splitlines()]
    def get_multilengths(lengths_file): 
        offset = []
        cur_offset = []
        for l in open(lengths_file, 'r').readlines():
            l = l.strip().split()
            if l[1] == '*':
                if cur_offset:
                    offset.append(cur_offset)
                cur_offset = []
                continue
            cur_offset.append(int(l[2]))
        offset.append(cur_offset)
        offset = np.array(offset)
        return offset
    simple = True
    try:
        with open(lengths_file, 'r') as f:
            first_line = f.readline().strip().split()
            if len(first_line) > 1 and first_line[1] == '*':
                simple = False
    except FileNotFoundError:
        raise ValueError("Either a *.lengths file or an input seq_lengths array is required")
    if simple and multilengths:
        raise ValueError("Multi-FASTA lengths not available in ", lengths_file)
    if not simple:
        offsets = get_multilengths(lengths_file)
        return offsets if multilengths else offsets.sum(axis=1).tolist()
    else:
        return get_lengths(lengths_file)

def unpack_flags(flags):
    """
    Unpack flags from a single uint64 value
    """
    flag_labels = ['partial', 'coll_blocks', 'merge']
    flags = np.unpackbits(np.array([flags], dtype=np.uint8))[-len(flag_labels):]
    return {f : bool(b) for f, b in zip(flag_labels, flags)}

def pack_flags(flags):
    """
    Pack flags into a single uint64 value
    """
    flag_labels = ['partial', 'coll_blocks', 'merge']
    flags = ([0] * (64 - len(flag_labels))) + [int(flags[f]) for f in flag_labels]
    return np.packbits(flags)

class MUMdata:
    def __init__(self, mumfile, lenfilter=0, subsample=1, verbose=False):
        if mumfile.endswith('.bums'):
            self.lengths, self.starts, self.strands = self.parse_bums(
                mumfile, 
                lenfilter, 
                subsample
            )
        else:
            self.lengths, self.starts, self.strands = self.parse_mums(
                mumfile, 
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
        self.partial = -1 in self.starts
    
    @classmethod
    def from_arrays(cls, lengths, starts, strands):
        """Create a MUMdata object directly from arrays.
        
        Args:
            lengths: Array of MUM lengths
            starts: 2D array of start positions (num_mums x num_seqs)
            strands: 2D array of strand information (num_mums x num_seqs)
        """
        instance = cls.__new__(cls)
        instance.lengths = lengths
        instance.starts = starts 
        instance.strands = strands
        instance.num_mums = len(lengths)
        instance.num_seqs = starts.shape[1] if instance.num_mums > 0 else 0
        return instance
        
    @staticmethod
    def parse_mums(mumfile, lenfilter=0, subsample=1, verbose=False):
        count = 0
        lengths, starts, strands = [], [], []
        with open(mumfile, 'r') as f:
            for line in tqdm(f, desc='parsing MUM file', disable=not verbose):
                if subsample == 1 or count % subsample == 0:
                    line = line.strip().split()
                    length = int(line[0])
                    if length >= lenfilter:
                        # Parse the line
                        strand = [s == '+' for s in line[2].split(',')]
                        start = [int(pos) if pos != '' else -1 for pos in line[1].split(',')]
                        starts.append(start)
                        strands.append(strand)
                        lengths.append(length)
                count += 1
        try:
            lengths = np.array(lengths, dtype=np.uint16)
        except OverflowError:
            raise ValueError("MUM length must be less than 65,535bp")
        try:
            starts = np.array(starts, dtype=np.int64)
        except OverflowError:
            raise ValueError("MUM start position must be less than 2^63")
        
        return lengths, starts, np.array(strands, dtype=bool)
        
    @staticmethod
    def parse_bums(bumfile, lenfilter=0, subsample=1):
        with open(bumfile, 'rb') as f:
            flags, n_seqs, n_mums = np.fromfile(f, count = 3, dtype=np.uint64)
            flags = unpack_flags(flags)
            start_dtype = np.int64
            mum_lengths = np.fromfile(f, count = n_mums, dtype=np.uint16)
            mum_starts = np.fromfile(f, count = n_seqs * n_mums, dtype=start_dtype).reshape(n_mums, n_seqs)
            mum_strands = np.fromfile(f, dtype=np.uint8)
            mum_strands = np.unpackbits(mum_strands, count=n_mums * n_seqs).reshape(n_mums, n_seqs).astype(bool)
        
        # Create boolean mask for subsampling
        mask = np.zeros(n_mums, dtype=bool)
        if subsample == 1:
            mask[:] = True
        else:
            mask[::subsample] = True            

        mask &= mum_lengths >= lenfilter
            
        return mum_lengths[mask], mum_starts[mask], mum_strands[mask]

    def filter_pmums(self):
        """Remove any MUMs that have -1 in their start positions"""
        if self.partial:
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
    
    def write_mums(self, filename, blocks=None):
        with open(filename, 'w') as f:
            if blocks is None:
                for i in range(self.num_mums):
                    strands_str = ['+' if s else '-' for s in self.strands[i]]
                    f.write(f"{self.lengths[i]}\t{','.join(map(str, self.starts[i]))}\t{','.join(strands_str)}\n")
            else:
                for idx, (l, r) in enumerate(blocks):
                    for i in range(l, r + 1):
                        strands_str = ['+' if s else '-' for s in self.strands[i]]
                        f.write(f"{self.lengths[i]}\t{','.join(map(str, self.starts[i]))}\t{','.join(strands_str)}\t{idx}\n")
    
    def write_bums(self, filename, blocks=None):
        if blocks:
            block_idx = [x for l,r in blocks for x in range(l, r + 1)]
            self.lengths = self.lengths[block_idx]
            self.starts = self.starts[block_idx]
            self.strands = self.strands[block_idx]
            self.num_mums = len(self.lengths)
        with open(filename, 'wb') as f:
            f.write(pack_flags({'partial': self.partial, 'coll_blocks': True if blocks else False, 'merge': False}).tobytes())
            f.write(np.uint64(self.num_seqs).tobytes())
            f.write(np.uint64(self.num_mums).tobytes())
            f.write(self.lengths.tobytes())
            f.write(self.starts.tobytes())
            f.write(self.strands.tobytes())
            if blocks:
                block_idx = np.array([idx for idx, (l,r) in enumerate(blocks) for _ in range(l, r + 1)], dtype=np.uint32)
                f.write(block_idx.tobytes())
