#!/usr/bin/env python3

import os, sys
import numpy as np
import argparse
from tqdm.auto import tqdm

from mumemto.utils import parse_bumbl_generator
try:
    from utils import stream_mums, unpack_flags, pack_flags, deserialize_coll_blocks, serialize_coll_blocks, MUMdata
except ImportError:
    from mumemto.utils import stream_mums, unpack_flags, pack_flags, deserialize_coll_blocks, serialize_coll_blocks, MUMdata

def parse_arguments(args=None):    
    parser = argparse.ArgumentParser(description="Converts Mumemto mum and bumbl formats")
    parser.add_argument('--mums', '-m', dest='mumfile', help='path to *.mum file from mumemto')
    parser.add_argument('--bumbl', '-b', dest='bumfile', help='path to *.bumbl file from mumemto')
    parser.add_argument('--length-upsize', '-l', dest='length_upsize', action='store_true', help='convert bumbl with 16 bit length to 32 bit (v1.3.1 onwards)')
    parser.add_argument('--fout','-o', dest='out', help='output fname')
    parser.add_argument('--verbose','-v', dest='verbose', help='verbose mode', action='store_true', default=False)
    parser.add_argument('--chunk-size','-c', dest='chunk_size', help='chunk size for writing output MUM file', default=8, type=int)
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)
    assert args.chunk_size % 8 == 0, "Chunk size must be a multiple of 8"
    args.to_bum = True
    if args.mumfile and args.bumfile and os.path.exists(args.mumfile) and os.path.exists(args.bumfile):
        parser.error("Multiple input files provided, only one is allowed")
        
    if args.out != None:
        if args.mumfile:
            args.to_bum = True
        else:
            args.to_bum = False
    elif args.mumfile and os.path.exists(args.mumfile) and args.bumfile and not os.path.exists(args.bumfile):
        args.to_bum = True
        args.out = args.bumfile
    elif args.bumfile and os.path.exists(args.bumfile) and args.mumfile and not os.path.exists(args.mumfile):
        args.to_bum = False
        args.out = args.mumfile
    elif args.out == None and args.bumfile == None and args.mumfile:
        args.out = args.mumfile.replace('.mums', '.bumbl')
        args.to_bum = True
    elif args.out == None and args.mumfile == None and args.bumfile:
        args.out = "-"
        args.to_bum = False
        
    if not args.to_bum and args.bumfile:
        if not os.path.exists(args.bumfile):
            parser.error("bumbl file does not exist")
    if args.to_bum and args.mumfile:
        if not os.path.exists(args.mumfile):
            parser.error("MUM file does not exist")            
    return args

def mum_to_bum(mumfile, outfile, verbose=False):
    length_dtype = np.uint32
    start_dtype = np.int64
    
    parser = stream_mums(mumfile, verbose=verbose, return_blocks=True)
    mum_count = 0
    lengths_out = open(outfile + '.len', 'wb')
    starts_out = open(outfile + '.starts', 'wb')
    strands_out = open(outfile + '.strands', 'wb')
    strands_array = []
    blocks_list = []
    is_partial = False
    for l, starts, strands, block in parser:
        # Write length as uint64
        if l > np.iinfo(length_dtype).max:
            raise ValueError("MUM length must be less than 2^32")
        lengths_out.write(length_dtype(l).tobytes())
        
        # Write starts as packed int64 array
        # Convert -1 to max value for the start_dtype
        starts = np.array(starts)
        if not is_partial and -1 in starts:
            is_partial = True
        starts_out.write(starts.astype(start_dtype).tobytes())
        
        # Write strands as packed bitvector
        strands_array.append(strands)
        
        if block is not None:
            blocks_list.append(block)
        mum_count += 1
    
    strands_array = np.array(strands_array, dtype=bool)
    strands_array = np.packbits(strands_array)
    strands_out.write(strands_array.tobytes())    
    
    lengths_out.close()
    starts_out.close()
    strands_out.close()
    
    flags = pack_flags({'partial': is_partial, 'coll_blocks': len(blocks_list) > 0, 'length32': True})
    
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
        if len(blocks_list) > 0:
            blocks = deserialize_coll_blocks(blocks_list)
            out.write(np.uint64(len(blocks)).tobytes())
            block_idx = np.array(blocks, dtype=np.uint32)
            out.write(block_idx.tobytes())
        
    # Clean up temporary files
    os.remove(outfile + '.len')
    os.remove(outfile + '.starts') 
    os.remove(outfile + '.strands')

def bum_to_mum(bumfile, outfile, verbose=False, chunk_size=8):
    if outfile == "-":
        outfile = sys.stdout
    else:
        outfile = open(outfile, 'w')
    
    try:
        with open(bumfile, "rb") as f:
            flags = int.from_bytes(f.read(2), byteorder='little')
            flags = unpack_flags(flags)
        
        # Use parse_bumbl_generator to read the data
        for lengths, starts, strands, blocks in parse_bumbl_generator(bumfile, verbose=verbose, return_chunk=True, return_blocks=True, chunksize=chunk_size):
            chunk_size = len(lengths)
            if not flags['coll_blocks']:
                outfile.write('\n'.join([f"{lengths[i]}\t{','.join(starts[i,:].astype(str))}\t{','.join(np.where(strands[i, :], '+', '-'))}" for i in range(chunk_size)]) + '\n')
            else:
                outfile.write('\n'.join([f"{lengths[i]}\t{','.join(starts[i,:].astype(str))}\t{','.join(np.where(strands[i, :], '+', '-'))}\t{blocks[i]}" for i in range(chunk_size)]) + '\n')
            
    except BrokenPipeError:
        # Python flushes standard streams on exit; redirect remaining output
        # to devnull to avoid another BrokenPipeError at shutdown
        sys.stdout = os.fdopen(sys.stdout.fileno(), 'wb', 0)
        sys.exit(0)

def length_upsize(args):
    if args.out == "-":
        args.out = args.bumfile.replace('.bumbl', '_32bit.bumbl')
    mum_lengths, mum_starts, mum_strands, blocks = MUMdata.parse_bums(args.bumfile, length_dtype=np.uint16, offset_dtype=np.int64)
    mum_lengths = mum_lengths.astype(np.uint32)
    mumdata = MUMdata.from_arrays(mum_lengths, mum_starts, mum_strands, blocks=blocks)
    if blocks is not None:
        mumdata.write_bums(args.out)
    else:
        mumdata.write_bums(args.out)

def main(args): 
    if args.length_upsize:
        length_upsize(args)
    elif args.to_bum:
        mum_to_bum(args.mumfile, args.out, verbose=args.verbose)
    else:
        bum_to_mum(args.bumfile, args.out, verbose=args.verbose, chunk_size=args.chunk_size)
    
if __name__ == "__main__":
    args = parse_arguments()
    main(args)
