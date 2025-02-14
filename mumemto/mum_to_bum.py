#!/usr/bin/env python3

import os, sys
import numpy as np
import argparse
from tqdm.auto import tqdm
try:
    from utils import parse_mums_generator, unpack_flags, pack_flags
except ImportError:
    from mumemto.utils import parse_mums_generator, unpack_flags, pack_flags

def parse_arguments(args=None):    
    parser = argparse.ArgumentParser(description="Plots a synteny plot of MUMs from mumemto")
    parser.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto')
    parser.add_argument('--bums', '-b', dest='bumfile', help='path to *.bum file from mumemto')

    parser.add_argument('--fout','-o', dest='out', help='output fname')
    parser.add_argument('--verbose','-v', dest='verbose', help='verbose mode', action='store_true', default=False)
    parser.add_argument('--chunk-size','-c', dest='chunk_size', help='chunk size for writing output MUM file', default=8, type=int)
    parser.add_argument('--64bit', dest='bit64', help='use 64-bit integers for start positions', action='store_true', default=False)
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)
    assert args.chunk_size % 8 == 0, "Chunk size must be a multiple of 8"
    args.to_bum = True
    if args.mumfile and args.bumfile and os.path.exists(args.mumfile) and os.path.exists(args.bumfile):
        parser.error("Multiple input files provided, only one is allowed")
        
    
    if args.mumfile and os.path.exists(args.mumfile) and args.bumfile and not os.path.exists(args.bumfile):
        args.to_bum = True
        args.out = args.bumfile
    elif args.bumfile and os.path.exists(args.bumfile) and args.mumfile and not os.path.exists(args.mumfile):
        args.to_bum = False
        args.out = args.mumfile
    elif args.out == None and args.bumfile == None and args.mumfile:
        args.out = args.mumfile.replace('.mums', '.bums')
        args.to_bum = True
    elif args.out == None and args.mumfile == None and args.bumfile:
        args.out = "-"
        args.to_bum = False
        
    if not args.to_bum and args.bumfile:
        if not os.path.exists(args.bumfile):
            parser.error("BUM file does not exist")
    if args.to_bum and args.mumfile:
        if not os.path.exists(args.mumfile):
            parser.error("MUM file does not exist")            
    return args

def mum_to_bum(mumfile, outfile, verbose=False, bit64=False):
    length_dtype = np.uint16
    start_dtype = np.int64
    
    parser = parse_mums_generator(mumfile, verbose=verbose)
    mum_count = 0
    lengths_out = open(outfile + '.len', 'wb')
    starts_out = open(outfile + '.starts', 'wb')
    strands_out = open(outfile + '.strands', 'wb')
    strands_array = []
    is_partial = False
    for l, starts, strands in parser:
        # Write length as uint64
        if l > 65535:
            raise ValueError("MUM length must be less than 65535")
        lengths_out.write(length_dtype(l).tobytes())
        
        # Write starts as packed int64 array
        # Convert -1 to max value for the start_dtype
        starts = np.array(starts)
        if not is_partial and -1 in starts:
            is_partial = True
        starts_out.write(starts.astype(start_dtype).tobytes())
        
        # Write strands as packed bitvector
        strands_array.append(strands)
        
        mum_count += 1
    
    strands_array = np.array(strands_array, dtype=bool)
    strands_array = np.packbits(strands_array)
    strands_out.write(strands_array.tobytes())    
    
    lengths_out.close()
    starts_out.close()
    strands_out.close()
    
    flags = pack_flags({'partial': is_partial, 'coll_blocks': False, 'merge': False})
    
    with open(outfile, 'wb') as out:
        out.write(flags.tobytes())
        out.write(np.uint64(len(strands)).tobytes())
        out.write(np.uint64(mum_count).tobytes())
        with open(outfile + '.len', 'rb') as f:
            out.write(f.read())
        with open(outfile + '.starts', 'rb') as f:
            out.write(f.read())
        with open(outfile + '.strands', 'rb') as f:
            out.write(f.read())
        
    # Clean up temporary files
    os.remove(outfile + '.len')
    os.remove(outfile + '.starts') 
    os.remove(outfile + '.strands')

def bum_to_mum(bumfile, outfile, verbose=False, chunk_size=8):
    length_size = 2
    length_dtype = np.uint16
    if outfile == "-":
        outfile = sys.stdout
    else:
        outfile = open(outfile, 'w')
    try:
        with open(bumfile, "rb") as f:
            # Read first byte and unpack into 8 bools
            flags = int.from_bytes(f.read(8), byteorder='little')
            flags = unpack_flags(flags)
            start_size = 8
            start_dtype = np.int64
            # Read num_seqs (1 byte) and num_mums (1 byte)
            num_seqs = int.from_bytes(f.read(8), byteorder='little')
            num_mums = int.from_bytes(f.read(8), byteorder='little')
            # Compute starting positions of each section
            lengths_pos = 8 * 3  # Immediately after flags, num_seqs, and num_mums
            offsets_pos = lengths_pos + (num_mums * length_size)  # After num_mums bytes
            strands_pos = offsets_pos + (num_seqs * num_mums * start_size)  # After offsets
            # Compute bit length for strands
            strand_buffer = None
            for mum_index in tqdm(range(0, num_mums, chunk_size), desc="Processing MUMs", disable=not verbose):
                f.seek(lengths_pos + (mum_index * length_size))
                # Read length (1 byte per mum)
                lengths = np.fromfile(f, count = chunk_size, dtype=length_dtype)

                # Seek to the correct offset for this mum
                offset_start = offsets_pos + (mum_index * num_seqs * start_size)
                f.seek(offset_start)
                starts = np.fromfile(f, count = num_seqs * chunk_size, dtype=start_dtype).reshape(chunk_size, num_seqs)
                
                if mum_index % chunk_size == 0:
                    f.seek(strands_pos + ((mum_index // chunk_size) * (chunk_size // 8) * num_seqs))
                    strand_buffer = np.fromfile(f, count = num_seqs * (chunk_size // 8), dtype=np.uint8)
                    strand_buffer = np.unpackbits(strand_buffer, count=num_seqs * chunk_size).reshape(chunk_size, num_seqs).astype(bool)
                strands = strand_buffer
                # Output the current mum as human-readable text
                if mum_index + chunk_size >= num_mums:
                    chunk_size = num_mums - mum_index
                outfile.write('\n'.join([f"{lengths[i]}\t{','.join(starts[i,:].astype(str))}\t{','.join(np.where(strands[i, :], '+', '-'))}" for i in range(chunk_size)]) + '\n')
                
    except BrokenPipeError:
        # Python flushes standard streams on exit; redirect remaining output
        # to devnull to avoid another BrokenPipeError at shutdown
        sys.stdout = os.fdopen(sys.stdout.fileno(), 'wb', 0)
        sys.exit(0)

def main(args): 
    if args.to_bum:
        mum_to_bum(args.mumfile, args.out, verbose=args.verbose, bit64=args.bit64)
    else:
        bum_to_mum(args.bumfile, args.out, verbose=args.verbose, chunk_size=args.chunk_size)
    
if __name__ == "__main__":
    args = parse_arguments()
    main(args)
