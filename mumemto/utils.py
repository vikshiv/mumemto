import numpy as np
from tqdm.auto import tqdm
from collections import namedtuple
import sys, os

MUM = namedtuple('MUM', ['length', 'starts', 'strands'])
MUM_BLOCK = namedtuple('MUM_BLOCK', ['length', 'starts', 'strands', 'block'])

def find_coll_blocks(mums, max_break=0, verbose=False, return_order=False, min_singleton_length=None):
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
    if verbose:
        print(f'Finding collinear blocks (max gap = {None if max_break == 0 else max_break} bp)...', file=sys.stderr)
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
        
    if min_singleton_length is not None:
        ### find singleton mums that are long enough to be single blocks    
        is_coll = np.zeros(len(lengths), dtype=bool)
        for s, e in blocks:
            is_coll[s:e+1] = True
        singleton = np.where((~is_coll) & (lengths >= min_singleton_length))
        for i in singleton[0]:
            blocks.append((i, i))
    
    blocks = sorted(blocks, key=lambda x: x[0])
    if return_order:
        order = mum_order_pos[:,[b[0] for b in blocks]].argsort(axis=1)
        return blocks, order
    return blocks

def get_coll_block_order(mums, blocks):
    return mums.starts[[b[0] for b in blocks],:].transpose().argsort(axis=1)
    
def parse_mums_generator(mumfile, seq_idx=None, verbose=False, return_blocks=False):
    """Generator that streams MUMs from mumfile"""
    with open(mumfile, 'r') as f:
        for line in tqdm(f, desc='parsing MUM file', disable=not verbose):
            line = line.strip().split()
            length = int(line[0])
            block = None if (len(line) < 4 or line[3] == '*') else line[3]
            # parse the full line
            if seq_idx is None:
                strands = [s == '+' for s in line[2].split(',')]
                starts = [int(pos) if pos != '' else -1 for pos in line[1].split(',')]
                yield MUM_BLOCK(length, starts, strands, block) if return_blocks else MUM(length, starts, strands)
            else:
                start = line[1].split(',')[seq_idx]
                ### only yield if mum appears in seq
                if start: 
                    start = int(start) if start != '' else -1
                    strand = line[2].split(',')[seq_idx] == '+'
                    yield MUM_BLOCK(length, start, strand, block) if return_blocks else MUM(length, start, strand)

def parse_first_mum(mumfile, verbose=False):
    """Special case, optimized parser to get MUM positions in the first sequence"""
    with open(mumfile, 'r') as f:
        for line in tqdm(f, desc='parsing MUM file', disable=not verbose):
            line = line.split()
            length = int(line[0])
            start = line[1][:line[1].index(',')]
            strand = line[2][:line[2].index(',')] == '+'
            ### only yield if mum appears in seq
            if start:
                start = int(start)
                yield (length, start, strand)

def parse_bumbl_generator(mumfile, seq_idx=None, verbose=False, chunksize=1024, return_chunk=False, return_blocks=False):
    """Generator that streams MUMs from bumbl file"""
    start_size = 8
    length_size = 4
    length_handle = open(mumfile, "rb")
    starts_handle = open(mumfile, "rb")
    strands_handle = open(mumfile, "rb")
    
    # Read flags
    flags_bytes = np.fromfile(length_handle, count=1, dtype=np.uint16)
    flags = unpack_flags(flags_bytes)
    
    # Read header information
    n_seqs, n_mums = np.fromfile(length_handle, count=2, dtype=np.uint64)
    
    # Calculate positions
    lengths_pos = 2 + 8 + 8  # After flags, num_seqs, and num_mums
    offsets_pos = lengths_pos + (n_mums * length_size)  # After lengths data
    strands_pos = offsets_pos + (n_mums * n_seqs * start_size)  # After starts data

    starts_handle.seek(offsets_pos)
    
    # Read all strands
    strands_handle.seek(strands_pos)
    strands_bytes = np.fromfile(strands_handle, count=np.ceil(n_seqs * n_mums / 8).astype(int), dtype=np.uint8)
    all_strands = np.unpackbits(strands_bytes, count=n_mums * n_seqs).reshape((n_mums, n_seqs)).astype(bool)
    
    # Check if coll_blocks flag is set to determine if blocks exist
    if return_blocks and flags.get('coll_blocks', False):
        # Read blocks data
        num_blocks = int.from_bytes(strands_handle.read(8), byteorder='little')
        all_blocks = np.fromfile(strands_handle, count=num_blocks * 2, dtype=np.uint32).reshape((num_blocks, 2))
        strands_handle.close()
        all_blocks = serialize_coll_blocks(all_blocks, n_mums)    
    else: 
        all_blocks = [None] * n_mums
        
    chunk = chunksize
    for idx in tqdm(range(0, n_mums, chunk), desc='parsing bumbl file', disable=not verbose):
        if idx + chunk > n_mums:
            chunk = n_mums - idx
        lengths = np.fromfile(length_handle, count=chunk, dtype=np.uint32)
        starts = np.fromfile(starts_handle, count=chunk * n_seqs, dtype=np.int64).reshape((chunk, n_seqs))
        strands = all_strands[idx:idx+chunk]
        blocks = all_blocks[idx:idx+chunk]
        
        if return_chunk:
            if seq_idx is None:
                yield (lengths, starts, strands) if not return_blocks else (lengths, starts, strands, blocks)
            else:
                yield (lengths, starts[:, seq_idx], strands[:, seq_idx]) if not return_blocks else (lengths, starts[:, seq_idx], strands[:, seq_idx], blocks)
        else:
            for i in range(chunk):
                if seq_idx is None:  
                    yield MUM_BLOCK(lengths[i], starts[i], strands[i], blocks[i]) if return_blocks else MUM(lengths[i], starts[i], strands[i])
                else:
                    start = starts[i, seq_idx]
                    if start != -1:  # Only yield if MUM appears in sequence
                        strand = strands[i, seq_idx]
                        yield MUM_BLOCK(lengths[i], start, strand, blocks[i]) if return_blocks else MUM(lengths[i], start, strand)
                            
    length_handle.close()
    starts_handle.close()
    strands_handle.close()

def stream_mums(mumfile, seq_idx=None, verbose=False, return_blocks=False):
    if mumfile.endswith('.mums') and seq_idx == 0 and not return_blocks:
        yield from parse_first_mum(mumfile, verbose=verbose)
    if mumfile.endswith('.mums'):
        yield from parse_mums_generator(mumfile, seq_idx=seq_idx, verbose=verbose, return_blocks=return_blocks)
    elif mumfile.endswith('.bumbl'):
        yield from parse_bumbl_generator(mumfile, seq_idx=seq_idx, verbose=verbose, return_blocks=return_blocks)
    else:
        raise ValueError('mumfile arg does not end with .mums or .bumbl')
    
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
        return offset
    simple = True
    try:
        with open(lengths_file, 'r') as f:
            first_line = f.readline().strip().split()
            if len(first_line) > 1 and first_line[1] == '*':
                simple = False
    except FileNotFoundError:
        raise FileNotFoundError(f"File {lengths_file} not found.")
    if simple and multilengths:
        raise ValueError("Multi-FASTA lengths not available in ", lengths_file)
    if not simple:
        offsets = get_multilengths(lengths_file)
        return offsets if multilengths else [sum(o) for o in offsets]
    else:
        return get_lengths(lengths_file)

def get_seq_paths(lengths_file):
    with open(lengths_file, 'r') as f:
        first_line = f.readline().strip().split()
        if len(first_line) > 1 and first_line[1] == '*':
            simple = False
        else:
            simple = True
    if simple:
        return [l.split()[0] for l in open(lengths_file, 'r').read().splitlines()]
    else:
        return [l.split()[0] for l in open(lengths_file, 'r').read().splitlines() if l.split()[1] == '*']

def get_contig_names(lengths_file):
    """
    Get contig names from a multilengths file.
    Returns a list of lists where each inner list contains the contig names for one sequence.
    
    Args:
        lengths_file: Path to the lengths file in multilengths format
        
    Returns:
        List of lists where names[i] is the list of contig names for sequence i
    """
    names = []
    cur_name = []
    first_line = True
    for l in open(lengths_file, 'r').readlines():
        l = l.strip().split()
        if first_line and l[1] != '*':
            raise ValueError('Lengths file must be formatted as multilengths.')
        first_line = False
        if l[1] == '*':
            if cur_name:
                names.append(cur_name)
            cur_name = []
            continue
        cur_name.append(l[1])
    names.append(cur_name)
    return names

def unpack_flags(packed_value):
    """
    Unpack a uint16 value into individual flags.
    """
    flag_labels = ['partial', 'coll_blocks', 'length32']
    # Convert uint16 to 16-bit binary representation
    bits = np.unpackbits(np.array([packed_value], dtype=np.uint16).view(np.uint8), bitorder='little')
    # Extract the last `len(flag_labels)` flags
    flags = {label: bool(bits[-(len(flag_labels) - i)]) for i, label in enumerate(flag_labels)}
    return flags

def pack_flags(flags):
    """
    Pack flags into a single uint16 value.
    """
    flag_labels = ['partial', 'coll_blocks', 'length32']
    bits = ([0] * (16 - len(flag_labels))) + [int(flags[f]) for f in flag_labels]
    # Pack into bytes and convert to uint16
    packed = np.packbits(bits, bitorder='little')
    return np.frombuffer(packed, dtype=np.uint16)[0]

def deserialize_coll_blocks(coll_blocks):
    coll_blocks = np.array([-1 if x == '-' else int(x) for x in coll_blocks])
    change_points = np.where(np.diff(coll_blocks) != 0)[0] + 1
    l_vals = np.concatenate(([0], change_points))
    r_vals = np.concatenate((change_points - 1, [len(coll_blocks) - 1]))
    valid_ranges = [(l, r) for l, r in zip(l_vals, r_vals) if coll_blocks[l] != -1]
    return valid_ranges

def serialize_coll_blocks(coll_blocks, num_mums):
    idx = 0
    block_idx = []
    left_block, right_block = coll_blocks[idx]
    for i in range(num_mums):
        if i > right_block:
            idx += 1
            if idx < len(coll_blocks):
                left_block, right_block = coll_blocks[idx]
        if i < left_block or i > right_block:
            block_idx.append('-')
        else:
            block_idx.append(str(idx))
    return block_idx

class MUMdata:
    def __init__(self, mumfile, lenfilter=0, subsample=1, sort=True, verbose=False):
        ### set types for each data type
        self.length_dtype = np.uint32
        self.offset_dtype = np.int64
        self.max_length = np.iinfo(self.length_dtype).max
        self.max_offset = np.iinfo(self.offset_dtype).max
        
        ### parse input file
        if mumfile.endswith('.bumbl'):
            self.lengths, self.starts, self.strands, self.blocks = self.parse_bums(
                mumfile, 
                lenfilter, 
                subsample,
                length_dtype=self.length_dtype,
                offset_dtype=self.offset_dtype
            )
            self.extra_fields = None
        else:
            self.lengths, self.starts, self.strands, self.blocks, self.extra_fields = self.parse_mums(
                mumfile, 
                lenfilter, 
                subsample,
                verbose,
                length_dtype=self.length_dtype,
                offset_dtype=self.offset_dtype
            )
        if sort:
            sorted = np.all(np.diff(self.starts[:,0]) >= 0)
            if self.blocks is not None and not sorted:
                print("MUMs must be sorted by first column to store blocks; ignoring blocks and sorting.", file=sys.stderr)
                self.blocks = None
            if not sorted:
            # sort by reference offset position
                order = self.starts[:,0].argsort()
                self.lengths = self.lengths[order]
                self.starts = self.starts[order]
                self.strands = self.strands[order]
                if self.extra_fields is not None:
                    self.extra_fields = [self.extra_fields[i] for i in order]
    
    @property
    def num_mums(self):
        """Number of MUMs in the dataset"""
        return len(self.lengths)
    
    @property
    def num_seqs(self):
        """Number of sequences in the dataset"""
        return self.starts.shape[1] if self.num_mums > 0 else 0
    
    @classmethod
    def from_arrays(cls, lengths, starts, strands, blocks=None, extra_fields=None):
        """Create a MUMdata object directly from arrays.
        
        Args:
            lengths: Array of MUM lengths
            starts: 2D array of start positions (num_mums x num_seqs)
            strands: 2D array of strand information (num_mums x num_seqs)
        """
        instance = cls.__new__(cls)
        instance.lengths = lengths
        instance.starts = starts.astype(np.int64, copy=False) 
        instance.strands = strands
        instance.blocks = blocks
        instance.extra_fields = extra_fields
        return instance
    
    def copy(self):
        """Create a deep copy of this MUMdata object with independent arrays.
        
        Returns:
            MUMdata: New MUMdata object with copied arrays
        """
        # Create copies of all arrays
        new_lengths = self.lengths.copy()
        new_starts = self.starts.copy()
        new_strands = self.strands.copy()
        
        # Handle extra_fields if present
        new_extra_fields = None
        if self.extra_fields is not None:
            new_extra_fields = self.extra_fields.copy()
        
        # Handle blocks if present
        new_blocks = None
        if self.blocks is not None:
            new_blocks = self.blocks.copy()
        
        # Create new instance using from_arrays class method
        new_mumdata = MUMdata.from_arrays(new_lengths, new_starts, new_strands, 
                                         blocks=new_blocks, extra_fields=new_extra_fields)
        return new_mumdata
        
    @staticmethod
    def parse_mums(mumfile, lenfilter=0, subsample=1, verbose=False, length_dtype=np.uint32, offset_dtype=np.int64):
        count = 0
        lengths, starts, strands, coll_blocks, extra_fields = [], [], [], [], []
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
                        if len(line) > 3 and line[3] != '*':
                            coll_blocks.append(line[3])
                        if len(line) > 4:
                            extra_fields.append('\t'.join(line[4:]))
                count += 1
                
        try:
            lengths = np.array(lengths, dtype=length_dtype)
        except OverflowError:
            raise ValueError("MUM length must be less than 2^32")
        
        try:
            starts = np.array(starts, dtype=offset_dtype)
        except OverflowError:
            raise ValueError("MUM start position must be less than 2^63")

        if len(coll_blocks) > 0:
            blocks = deserialize_coll_blocks(coll_blocks)
        else:
            blocks = None
        
        if len(extra_fields) == 0:
            extra_fields = None
        
        return lengths, starts, np.array(strands, dtype=bool), blocks, extra_fields
    
    @staticmethod
    def parse_bums(bumfile, lenfilter=0, subsample=1, length_dtype=np.uint32, offset_dtype=np.int64):
        filesize = os.path.getsize(bumfile)
        with open(bumfile, 'rb') as f:
            flags = np.fromfile(f, count = 1, dtype=np.uint16)
            n_seqs, n_mums = np.fromfile(f, count = 2, dtype=np.uint64)
            flags = unpack_flags(flags)
            mum_lengths = np.fromfile(f, count = n_mums, dtype=length_dtype)
            mum_starts = np.fromfile(f, count = n_seqs * n_mums, dtype=offset_dtype).reshape(n_mums, n_seqs)
            mum_strands = np.fromfile(f, count=np.ceil(n_seqs*n_mums/8).astype(int), dtype=np.uint8)
            mum_strands = np.unpackbits(mum_strands, count=n_mums * n_seqs).reshape(n_mums, n_seqs).astype(bool)
            if flags['coll_blocks']:
                num_blocks = int.from_bytes(f.read(8), byteorder='little')
                blocks = np.fromfile(f, count=num_blocks * 2, dtype=np.uint32).reshape(num_blocks, 2)
            else:
                blocks = None
    
        # Create boolean mask for subsampling
        mask = np.zeros(n_mums, dtype=bool)
        if subsample == 1:
            mask[:] = True
        else:
            mask[::subsample] = True            

        mask &= mum_lengths >= lenfilter
            
        return mum_lengths[mask], mum_starts[mask], mum_strands[mask], blocks
    
    def filter_pmums(self):
        """Remove any MUMs that have -1 in their start positions"""
        if -1 in self.starts:
            valid_rows = ~np.any(self.starts == -1, axis=1)
            self.lengths = self.lengths[valid_rows]
            self.starts = self.starts[valid_rows]
            self.strands = self.strands[valid_rows]
            if self.extra_fields is not None:
                self.extra_fields = [self.extra_fields[i] for i in valid_rows]
        return self

    def slice(self, indices, copy=False):
        """Create a new MUMdata object with numpy-style slicing.
           Supports: mumdata[rows], mumdata[:, cols], mumdata[rows, cols]
           The new MUMdata object will not have blocks
        
        Args:
            indices: Single index/slice/mask for rows, or tuple (rows, cols) for 2D slicing
            copy: If True, create copies of the arrays instead of views
            
        Returns:
            MUMdata: New MUMdata object containing only the selected MUMs and sequences
        """
        # Handle tuple indexing (2D slicing)
        if isinstance(indices, tuple):
            if len(indices) == 2:
                row_indices, col_indices = indices
            elif len(indices) == 1:
                row_indices = indices[0]
                col_indices = slice(None)
            else:
                raise ValueError("Too many indices for 2D slicing")
        else:
            # Single index - slice only rows
            row_indices = indices
            col_indices = slice(None)
        
        # Apply indexing directly (numpy handles all the cases)
        new_lengths = self.lengths[row_indices]
        new_starts = self.starts[row_indices, col_indices]
        new_strands = self.strands[row_indices, col_indices]
        
        # Handle extra_fields if present
        new_extra_fields = None
        if self.extra_fields is not None:
            new_extra_fields = [self.extra_fields[i] for i in np.arange(self.num_mums)[row_indices]]
        
        # Create new instance using from_arrays class method
        new_mumdata = MUMdata.from_arrays(new_lengths, new_starts, new_strands, extra_fields=new_extra_fields)
        if copy:
            new_mumdata = new_mumdata.copy()
            
        return new_mumdata

    def __getitem__(self, idx):
        """Get MUM(s) by index, slice, or list of indices.
        Supports numpy-style slicing: mumdata[rows], mumdata[:, cols], mumdata[rows, cols]
        
        Args:
            idx: Integer, slice, list, array, boolean array, or tuple for 2D slicing
            
        Returns:
            MUM: Single MUM if idx is integer
            MUMdata: New MUMdata object if idx is slice/list/array/tuple
        """
        if isinstance(idx, (int, np.integer)):
            # Single MUM access - return MUM object
            return MUM(self.lengths[idx], self.starts[idx], self.strands[idx])
        else:
            # Multiple MUMs access - return new MUMdata object using slice method
            return self.slice(idx)

    def __repr__(self):
        """String representation"""
        return (f"MUMdata(num_mums={self.num_mums}, num_seqs={self.num_seqs}")
    
    def __str__(self):
        """Human-readable string representation"""
        return f"MUMdata with {self.num_mums} MUMs across {self.num_seqs} sequences"
    
    def __bool__(self):
        """Boolean evaluation - True if has MUMs, False if empty"""
        return self.num_mums > 0
    
    def __iter__(self):
        """Iterate over MUMs, yielding (length, starts, strands) for each"""
        for i in range(self.num_mums):
            yield self[i]
            
    def __add__(self, other):
        """Concatenate two MUMdata objects. Will ignore extra fields."""
        if not isinstance(other, MUMdata):
            raise TypeError("Can only concatenate MUMdata objects")
        
        if self.num_seqs != other.num_seqs:
            raise ValueError("Cannot concatenate MUMdata objects with different numbers of sequences")
        
        # Concatenate arrays
        new_lengths = np.concatenate([self.lengths, other.lengths])
        new_starts = np.vstack([self.starts, other.starts])
        new_strands = np.vstack([self.strands, other.strands])
        
        # Create new instance
        new_mumdata = MUMdata.from_arrays(new_lengths, new_starts, new_strands)
        return new_mumdata

    def __len__(self):
        """Return number of MUMs"""
        return self.num_mums
    
    def filter_length(self, length):
        """Filter MUMs < length threshold"""
        return self[self.lengths < length]
        
    def write_mums(self, filename, blocks=None):
        with open(filename, 'w') as f:
            if blocks is None:
                for i in range(self.num_mums):
                    strands_str = ['+' if s else '-' for s in self.strands[i]]
                    f.write(f"{self.lengths[i]}\t{','.join(map(str, self.starts[i]))}\t{','.join(strands_str)}\n")
            else:
                if not np.all(np.diff(self.starts[:,0]) >= 0):
                    print("MUMs must be sorted by first column to write blocks; ignoring blocks.", file=sys.stderr)
                else:
                    idx = 0
                    block_idx = 0
                    left_block, right_block = blocks[idx]
                    for i in range(self.num_mums):
                        if i > right_block:
                            idx += 1
                            if idx < len(blocks):
                                left_block, right_block = blocks[idx]
                        if i < left_block or i > right_block:
                            block_idx = '-'
                        else:
                            block_idx = idx
                        strands_str = ['+' if s else '-' for s in self.strands[i]]
                        if self.extra_fields is not None:
                            f.write(f"{self.lengths[i]}\t{','.join(map(str, self.starts[i]))}\t{','.join(strands_str)}\t{block_idx}\t{self.extra_fields[i]}\n")
                        else:
                            f.write(f"{self.lengths[i]}\t{','.join(map(str, self.starts[i]))}\t{','.join(strands_str)}\t{block_idx}\n")
    
    def write_bums(self, filename, blocks=None):
        with open(filename, 'wb') as f:
            f.write(pack_flags({'partial': -1 in self.starts, 
                                'coll_blocks': blocks is not None, 
                                'length32': self.lengths.dtype == np.uint32}).tobytes())
            f.write(np.uint64(self.num_seqs).tobytes())
            f.write(np.uint64(self.num_mums).tobytes())
            f.write(self.lengths.tobytes())
            f.write(self.starts.tobytes())
            strands_array = np.packbits(self.strands)
            f.write(strands_array.tobytes())
            if blocks is not None:
                if not np.all(np.diff(self.starts[:,0]) >= 0):
                    print("MUMs must be sorted by first column to write blocks; ignoring blocks.", file=sys.stderr)
                else:
                    f.write(np.uint64(len(blocks)).tobytes())
                    block_idx = np.array(blocks, dtype=np.uint32)
                    f.write(block_idx.tobytes())
