#!/usr/bin/env python3
"""
Convert mum/bumbl to a mum-like plaintext with three BED columns (contig, start, end)
for a chosen sequence, then bgzip and tabix index.
"""

import os
import sys
import argparse
import subprocess
import tempfile
import numpy as np
from tqdm.auto import tqdm

from mumemto.utils import (
    get_sequence_lengths,
    get_contig_names,
    parse_bumbl_generator,
    unpack_flags,
)


def get_lengths_file(path):
    """Get the corresponding lengths file path for a .mums or .bumbl file."""
    base = os.path.splitext(path)[0]
    lengths_file = f"{base}.lengths"
    if not os.path.exists(lengths_file):
        raise FileNotFoundError(f"Lengths file {lengths_file} not found")
    return lengths_file


def find_chr_one(start, length, lengths):
    """Return (contig_idx, rel_start, rel_end) for one interval. BED 0-based."""
    offsets = np.cumsum(lengths)
    contig_idx = np.searchsorted(offsets, start, side='right')
    if contig_idx >= len(offsets):
        contig_idx = len(offsets) - 1
    left_start = np.hstack((0, offsets[:-1]))
    rel_start = start - left_start[contig_idx]
    rel_end = rel_start + length
    return contig_idx, int(rel_start), int(rel_end)


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description="Convert mum/bumbl to mum-like plaintext + contig/start/end columns, then bgzip and tabix."
    )
    parser.add_argument('input', nargs='?', help='Path to .mums or .bumbl file')
    parser.add_argument('--mums', '-m', dest='mumfile', help='Path to .mums file')
    parser.add_argument('--bumbl', '-b', dest='bumfile', help='Path to .bumbl file')
    parser.add_argument('--seq-idx', '-s', type=int, default=0,
                        help='Sequence index for BED coordinates (default: 0)')
    parser.add_argument('--output', '-o', dest='out', help='Output path (default: input.mum.bed.gz)')
    parser.add_argument('--lengths-file', '-l', dest='lengths_file', help='Path to .lengths file (optional)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Verbose / progress')
    parser.add_argument('--no-tabix', action='store_true', help='Only write .bed.gz, do not run tabix')
    parser.add_argument('--chunk-size', '-c', type=int, default=1024,
                        help='Chunk size when reading bumbl (default: 1024)')

    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    # Resolve input: positional or -m/-b
    if args.input and (args.mumfile or args.bumfile):
        parser.error("Use either positional input or --mums/--bumbl, not both")
    if args.input:
        if args.input.endswith('.bumbl'):
            args.bumfile = args.input
            args.mumfile = None
        elif args.input.endswith('.mums'):
            args.mumfile = args.input
            args.bumfile = None
        else:
            parser.error("Input must end with .mums or .bumbl")
    if not args.mumfile and not args.bumfile:
        parser.error("Provide an input file (positional, --mums, or --bumbl)")
    if args.mumfile and args.bumfile:
        parser.error("Provide only one of --mums or --bumbl")

    inp = args.mumfile or args.bumfile
    if not os.path.exists(inp):
        parser.error(f"Input file does not exist: {inp}")

    if args.lengths_file is None:
        args.lengths_file = get_lengths_file(inp)

    if args.out is None:
        base = os.path.splitext(inp)[0]
        if inp.endswith('.bumbl'):
            base = os.path.splitext(base)[0]  # file.bumbl -> file
        args.out = base + '.mum.bed.gz'

    return args


def mum_to_tabix_mums(mumfile, lengths_file, seq_idx, out_gz, verbose=False, run_tabix=True):
    """Read .mums line by line; append contig, start, end for seq_idx; write plain then bgzip + tabix."""
    lengths_list = get_sequence_lengths(lengths_file, multilengths=True)
    if seq_idx >= len(lengths_list):
        raise ValueError(f"seq_idx {seq_idx} >= number of sequences {len(lengths_list)}")
    lengths = np.array(lengths_list[seq_idx])
    contig_names = get_contig_names(lengths_file)[seq_idx]

    fd, plain_path = tempfile.mkstemp(suffix='.mum.bed', prefix='mum_to_tabix_')
    try:
        with os.fdopen(fd, 'w') as out:
            with open(mumfile, 'r') as f:
                lines = tqdm(f, desc='mums', disable=not verbose)
                for line in lines:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split('\t')
                    if len(parts) < 3:
                        continue
                    length = int(parts[0])
                    starts_str = parts[1].split(',')
                    if seq_idx >= len(starts_str) or not starts_str[seq_idx]:
                        continue
                    start_str = starts_str[seq_idx].strip()
                    if start_str == '' or start_str == '-1':
                        continue
                    start = int(start_str)
                    contig_i, rel_start, rel_end = find_chr_one(start, length, lengths)
                    contig = contig_names[contig_i]
                    out.write(f"{line}\t{contig}\t{rel_start}\t{rel_end}\n")
    except Exception:
        os.unlink(plain_path)
        raise

    _compress_and_index(plain_path, out_gz, run_tabix=run_tabix)
    os.unlink(plain_path)


def mum_to_tabix_bumbl(bumfile, lengths_file, seq_idx, out_gz, verbose=False, run_tabix=True, chunk_size=1024):
    """Read .bumbl in chunks; for each MUM append contig, start, end for seq_idx; write plain then bgzip + tabix."""
    lengths_list = get_sequence_lengths(lengths_file, multilengths=True)
    if seq_idx >= len(lengths_list):
        raise ValueError(f"seq_idx {seq_idx} >= number of sequences {len(lengths_list)}")
    lengths = np.array(lengths_list[seq_idx])
    contig_names = get_contig_names(lengths_file)[seq_idx]

    with open(bumfile, 'rb') as f:
        flags = unpack_flags(int.from_bytes(f.read(2), byteorder='little'))
    has_blocks = flags.get('coll_blocks', False)

    fd, plain_path = tempfile.mkstemp(suffix='.mum.bed', prefix='mum_to_tabix_')
    try:
        with os.fdopen(fd, 'w') as out:
            for lengths_chunk, starts_chunk, strands_chunk, blocks_chunk in parse_bumbl_generator(
                bumfile, seq_idx=None, verbose=verbose, return_chunk=True, return_blocks=True, chunksize=chunk_size
            ):
                n = len(lengths_chunk)
                for i in range(n):
                    start = int(starts_chunk[i, seq_idx])
                    if start == -1:
                        continue
                    length = int(lengths_chunk[i])
                    contig_i, rel_start, rel_end = find_chr_one(start, length, lengths)
                    contig = contig_names[contig_i]
                    # Mum-like line: length \t starts_csv \t strands_csv [\t block]
                    starts_csv = ','.join(starts_chunk[i, :].astype(str))
                    strands_csv = ','.join('+' if b else '-' for b in strands_chunk[i, :])
                    if has_blocks and blocks_chunk is not None:
                        block_str = str(blocks_chunk[i]) if blocks_chunk[i] is not None else '-'
                        out.write(f"{length}\t{starts_csv}\t{strands_csv}\t{block_str}\t{contig}\t{rel_start}\t{rel_end}\n")
                    else:
                        out.write(f"{length}\t{starts_csv}\t{strands_csv}\t{contig}\t{rel_start}\t{rel_end}\n")
    except Exception:
        os.unlink(plain_path)
        raise

    _compress_and_index(plain_path, out_gz, run_tabix=run_tabix)
    os.unlink(plain_path)


def _compress_and_index(plain_path, out_gz, run_tabix=True):
    """Run bgzip on plain_path -> out_gz, then tabix if run_tabix."""
    with open(plain_path, 'r') as f:
        first = f.readline()
    ncols = len(first.split('\t'))
    seq_col = ncols - 3  # 0-based
    start_col = ncols - 2
    end_col = ncols - 1

    try:
        with open(out_gz, 'wb') as dest:
            subprocess.run(
                ['bgzip', '-c', '-f', plain_path],
                check=True,
                stdout=dest,
                stderr=subprocess.DEVNULL,
            )
    except FileNotFoundError:
        # Fallback: write with gzip if bgzip not available
        import gzip
        with open(plain_path, 'rb') as src:
            with gzip.open(out_gz, 'wb') as dst:
                dst.writelines(src)
        if run_tabix:
            sys.stderr.write("bgzip not found; wrote gzip. Tabix requires bgzip; skipping index.\n")
        return
    except Exception as e:
        sys.stderr.write(f"bgzip failed: {e}\n")
        raise

    if run_tabix:
        try:
            subprocess.run(
                ['tabix', '-s', str(seq_col + 1), '-b', str(start_col + 1), '-e', str(end_col + 1), '-f', out_gz],
                check=True,
                capture_output=True,
            )
        except FileNotFoundError:
            sys.stderr.write("tabix not found; index not created.\n")
        except subprocess.CalledProcessError as e:
            sys.stderr.write(f"tabix failed: {e.stderr.decode() if e.stderr else e}\n")


def main(args):
    if args.mumfile:
        mum_to_tabix_mums(
            args.mumfile,
            args.lengths_file,
            args.seq_idx,
            args.out,
            verbose=args.verbose,
            run_tabix=not args.no_tabix,
        )
    else:
        mum_to_tabix_bumbl(
            args.bumfile,
            args.lengths_file,
            args.seq_idx,
            args.out,
            verbose=args.verbose,
            run_tabix=not args.no_tabix,
            chunk_size=args.chunk_size,
        )
    if args.verbose:
        print(f"Wrote {args.out}", file=sys.stderr)
        if not args.no_tabix and os.path.exists(args.out + '.tbi'):
            print(f"Index: {args.out}.tbi", file=sys.stderr)


if __name__ == '__main__':
    args = parse_arguments()
    main(args)
