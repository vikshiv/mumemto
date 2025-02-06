#!/usr/bin/env python3

import os, sys
import numpy as np
import argparse
from tqdm.auto import tqdm
import math
try:
    from utils import parse_mums_generator
except ImportError:
    from mumemto.utils import parse_mums_generator

def parse_arguments(args=None):    
    parser = argparse.ArgumentParser(description="Plots a synteny plot of MUMs from mumemto")
    # parser.add_argument('--filelist', '-f', dest='filelist', help='path to filelist from mumemto', required=True)
    # parser.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto', required=True)
    # parser.add_argument('--lengths','-l', dest='lens', help='lengths file, first column is seq length in order of filelist', required=True)
    # group = parser.add_mutually_exclusive_group(required=True)
    parser.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto')
    parser.add_argument('--bums', '-b', dest='bumfile', help='path to *.bum file from mumemto')

    parser.add_argument('--fout','-o', dest='out', help='output fname')
    parser.add_argument('--verbose','-v', dest='verbose', help='verbose mode', action='store_true', default=False)
    parser.add_argument('--len-filter','-L', dest='lenfilter', help='only plot MUMs longer than threshold', default=0, type=int)
    parser.add_argument('--subsample','-s', dest='subsample', help='subsample every Nth mum', default=1, type=int)
    parser.add_argument('--chunk-size','-c', dest='chunk_size', help='chunk size for writing output MUM file', default=8, type=int)
    
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

def mum_to_bum(mumfile, outfile, lenfilter=0, subsample=1, verbose=False):
    parser = parse_mums_generator(mumfile, lenfilter=lenfilter, subsample=subsample, verbose=verbose)
    mum_count = 0
    lengths_out = open(outfile + '.len', 'wb')
    starts_out = open(outfile + '.starts', 'wb')
    strands_out = open(outfile + '.strands', 'wb')
    strands_array = []
    for l, starts, strands in parser:
        # Write length as uint64
        lengths_out.write(np.uint64(l).tobytes())
        
        # Write starts as packed int64 array
        starts_array = np.array(starts, dtype=np.int64)
        starts_out.write(starts_array.tobytes())
        
        # Write strands as packed bitvector
        strands_array.append(strands)
        
        mum_count += 1
    
    strands_array = np.array(strands_array, dtype=bool)
    strands_array = np.packbits(strands_array)
    strands_out.write(strands_array.tobytes())    
    
    lengths_out.close()
    starts_out.close()
    strands_out.close()
    
    with open(outfile, 'wb') as out:
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
        

# def bum_to_mum(bumfile, outfile, lenfilter=0, subsample=1, verbose=False):
#     try:
#         with open(bumfile, "rb") as f:
#             # Read num_seqs (1 byte) and num_mums (1 byte)
#             num_seqs = int.from_bytes(f.read(8), byteorder='little')
#             num_mums = int.from_bytes(f.read(8), byteorder='little')
#             # Compute starting positions of each section
#             lengths_pos = 16  # Immediately after num_seqs and num_mums
#             offsets_pos = lengths_pos + (num_mums * 8)  # After num_mums bytes
#             strands_pos = offsets_pos + (num_seqs * num_mums * 8)  # After offsets
#             # Compute bit length for strands
#             strands_byte_size = math.ceil(num_seqs * num_mums / 8)
#             strand_buffer = None
#             for mum_index in range(num_mums):
#                 f.seek(lengths_pos + (mum_index * 8))
#                 # Read length (1 byte per mum)
#                 length = int.from_bytes(f.read(8), byteorder='little')

#                 # Seek to the correct offset for this mum
#                 offset_start = offsets_pos + (mum_index * num_seqs * 8)
#                 f.seek(offset_start)
#                 starts = np.fromfile(f, count = num_seqs, dtype=np.int64)
                
#                 if mum_index % 8 == 0:
#                     f.seek(strands_pos + ((mum_index // 8) * num_seqs))
#                     strand_buffer = np.fromfile(f, count = num_seqs, dtype=np.uint8)
#                     strand_buffer = np.unpackbits(strand_buffer, count=num_seqs * 8).reshape(8, num_seqs).astype(bool)
#                 strands = strand_buffer[mum_index % 8, :]
                
#                 # Output the current mum as human-readable text
#                 print(f"{length}\t{','.join(starts.astype(str))}\t{','.join(np.where(strands, '+', '-'))}")
                
#     except BrokenPipeError:
#         # Python flushes standard streams on exit; redirect remaining output
#         # to devnull to avoid another BrokenPipeError at shutdown
#         import os
#         import sys
#         sys.stdout = os.fdopen(sys.stdout.fileno(), 'wb', 0)
#         sys.exit(0)

def bum_to_mum(bumfile, outfile, lenfilter=0, subsample=1, verbose=False, chunk_size=8):
    length_size, start_size = 8, 8
    if outfile == "-":
        outfile = sys.stdout
    else:
        outfile = open(outfile, 'w')
    try:
        with open(bumfile, "rb") as f:
            # Read num_seqs (1 byte) and num_mums (1 byte)
            num_seqs = int.from_bytes(f.read(8), byteorder='little')
            num_mums = int.from_bytes(f.read(8), byteorder='little')
            # Compute starting positions of each section
            lengths_pos = 16  # Immediately after num_seqs and num_mums
            offsets_pos = lengths_pos + (num_mums * length_size)  # After num_mums bytes
            strands_pos = offsets_pos + (num_seqs * num_mums * start_size)  # After offsets
            # Compute bit length for strands
            strand_buffer = None
            for mum_index in tqdm(range(0, num_mums, chunk_size), desc="Processing MUMs", disable=not verbose):
                f.seek(lengths_pos + (mum_index * length_size))
                # Read length (1 byte per mum)
                lengths = np.fromfile(f, count = chunk_size, dtype=np.uint64)

                # Seek to the correct offset for this mum
                offset_start = offsets_pos + (mum_index * num_seqs * start_size)
                f.seek(offset_start)
                starts = np.fromfile(f, count = num_seqs * chunk_size, dtype=np.int64).reshape(chunk_size, num_seqs)
                
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

# def bum_to_mum(bumfile, outfile, lenfilter=0, subsample=1, verbose=False):
#     try:
#         with open(bumfile, "rb") as lengths_file, open(bumfile, "rb") as starts_file, open(bumfile, "rb") as strands_file:
#             num_seqs, num_mums = np.fromfile(lengths_file, count = 2, dtype=np.uint64)
#             # Compute starting positions of each section
#             lengths_pos = 2  # Immediately after num_seqs and num_mums
#             offsets_pos = lengths_pos + (num_mums * 8)  # After num_mums bytes
#             strands_pos = offsets_pos + (num_seqs * num_mums * 8)  # After offsets
            
#             # starts_file.seek(offsets_pos)
#             strands_file.seek(strands_pos)
#             strands_count = 8
#             strand_buffer = None
#             lengths = np.fromfile(lengths_file, count = num_mums, dtype=np.uint64)
#             for mum_index in range(num_mums):
#                 # Read length (1 byte per mum)
                
#                 # starts_file.seek(-1, 1)
#                 starts = np.fromfile(lengths_file, count = num_seqs, dtype=np.int64)
                
#                 if strands_count == 8:
#                     strand_buffer = np.fromfile(strands_file, count = num_seqs, dtype=np.uint8)
#                     strand_buffer = np.unpackbits(strand_buffer, count=num_seqs * 8).reshape(8, num_seqs).astype(bool)
#                     strands_count = 0
                
#                 strands = strand_buffer[strands_count, :]
#                 strands_count += 1

#                 # Output the current mum as human-readable text
#                 print(f"{lengths[mum_index]}\t{','.join(starts.astype(str))}\t{','.join(np.where(strands, '+', '-'))}")
#     except BrokenPipeError:
#         # Python flushes standard streams on exit; redirect remaining output
#         # to devnull to avoid another BrokenPipeError at shutdown
#         import os
#         import sys
#         sys.stdout = os.fdopen(sys.stdout.fileno(), 'wb', 0)
#         sys.exit(0)
def main(args): 
    if args.to_bum:
        mum_to_bum(args.mumfile, args.out, lenfilter=args.lenfilter, subsample=args.subsample, verbose=args.verbose)
    else:
        bum_to_mum(args.bumfile, args.out, lenfilter=args.lenfilter, subsample=args.subsample, verbose=args.verbose, chunk_size=args.chunk_size)
    
if __name__ == "__main__":
    args = parse_arguments()
    main(args)
