#!/usr/bin/env python3
"""Convert mumemto .mums/.bumbl + multilengths .lengths to ASCII ONEcode .1mum.

Primary type is ``mum``. This script writes ASCII ONEcode content to a ``.1mum``
file for convenience; ONEcode's usual ASCII/binary pair would be ``.mum`` /
``.1mum``, and binary form can later be produced with ONEview/ONElib.

Schema (embedded in each output file)::

    P 3 mum
    O L 1 3 INT              # MUM match Length
    D P 1 8 INT_LIST         # Positions (abs starts per genome; -1 = absent)
    D S 1 6 STRING           # Strands as +/- string
    D B 1 3 INT              # collinear Block id (-1 if none)
    D X 1 6 STRING           # optional eXtra fields from .mums
    D G 1 3 INT              # Genome total length
    D F 1 6 STRING           # Fasta path
    D C 1 8 INT_LIST         # Contig lengths
    D N 1 11 STRING_LIST     # contig Names

Note: ONEcode ``>`` / ``<`` header lines reference other ONEcode files (forward /
backward object deps), not a list of fasta paths. Genome paths are stored as
``F`` data lines instead.

Requires a multilengths .lengths file; simple lengths format is rejected.
"""

from __future__ import annotations

import argparse
import os
import sys
from datetime import datetime

try:
    from utils import (
        MUMdata,
        get_contig_names,
        get_seq_paths,
        get_sequence_lengths,
        serialize_coll_blocks,
    )
except ImportError:
    from mumemto.utils import (
        MUMdata,
        get_contig_names,
        get_seq_paths,
        get_sequence_lengths,
        serialize_coll_blocks,
    )

SCRIPT_NAME = "mums_to_1mum"
SCRIPT_VERSION = "0.1.0"


def one_string(s: str) -> str:
    """Encode a STRING as length-prefixed ONEcode tokens."""
    return f"{len(s)} {s}"


def one_int_list(values) -> str:
    """Encode an INT_LIST as length-prefixed ONEcode tokens."""
    vals = list(values)
    return f"{len(vals)} " + " ".join(str(int(v)) for v in vals)


def one_string_list(strings) -> str:
    """Encode a STRING_LIST as length-prefixed ONEcode tokens."""
    parts = [str(len(strings))]
    for s in strings:
        parts.append(one_string(s))
    return " ".join(parts)


def parse_arguments(args=None):
    parser = argparse.ArgumentParser(
        description="Convert .mums/.bumbl + multilengths .lengths to ASCII ONEcode .1mum"
    )
    parser.add_argument(
        "mumfile",
        help="path to *.mums or *.bumbl file from mumemto",
    )
    parser.add_argument(
        "-l",
        "--lengths",
        dest="lengths",
        help="multilengths .lengths file (default: <prefix>.lengths)",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        help="output .1mum path (default: <prefix>.1mum)",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="verbose progress on stderr",
    )

    args = parser.parse_args(args)

    mumfile = args.mumfile
    if not mumfile.endswith(".mums") and not mumfile.endswith(".bumbl"):
        if os.path.exists(mumfile + ".mums"):
            mumfile = mumfile + ".mums"
        elif os.path.exists(mumfile + ".bumbl"):
            mumfile = mumfile + ".bumbl"
        else:
            print(f"MUM file {args.mumfile} not found.", file=sys.stderr)
            sys.exit(1)
    elif not os.path.exists(mumfile):
        print(f"MUM file {mumfile} not found.", file=sys.stderr)
        sys.exit(1)
    args.mumfile = mumfile

    prefix = os.path.splitext(args.mumfile)[0]
    if args.lengths is None:
        args.lengths = prefix + ".lengths"
    if not os.path.exists(args.lengths):
        print(f"Lengths file {args.lengths} not found.", file=sys.stderr)
        sys.exit(1)

    if args.output is None:
        args.output = prefix + ".1mum"

    return args


def write_1mum(outfile, mums, paths, contig_lengths, contig_names, command, verbose=False):
    n_genomes = len(paths)
    n_mums = mums.num_mums

    if mums.num_seqs != n_genomes:
        raise ValueError(
            f"MUM file has {mums.num_seqs} sequences but lengths file has {n_genomes} genomes"
        )

    genome_totals = [sum(c) for c in contig_lengths]
    has_blocks = mums.blocks is not None
    has_extras = mums.extra_fields is not None

    block_ids = None
    if has_blocks:
        block_ids = [
            -1 if b == "-" else int(b)
            for b in serialize_coll_blocks(mums.blocks, n_mums)
        ]

    # Size stats for STRING / STRING_LIST headers
    path_lens = [len(p) for p in paths]
    name_line_sums = [sum(len(n) for n in names) for names in contig_names]
    total_contigs = sum(len(c) for c in contig_lengths)
    max_contigs = max((len(c) for c in contig_lengths), default=0)
    extra_lens = [len(x) for x in mums.extra_fields] if has_extras else []

    if verbose:
        print(
            f"Writing {n_mums} MUMs across {n_genomes} genomes to {outfile}",
            file=sys.stderr,
        )

    with open(outfile, "w") as f:
        # Header: primary type + provenance + schema + counts
        f.write(f"1 {one_string('mum')} 1 0\n")
        date = datetime.now().strftime("%Y-%m-%d_%H:%M:%S")
        f.write(
            f"! {one_string(SCRIPT_NAME)} {one_string(SCRIPT_VERSION)} "
            f"{one_string(command)} {one_string(date)}\n"
        )

        f.write(f"~ P {one_string('mum')}\n")
        f.write("~ O L 1 3 INT\n")
        f.write("~ D P 1 8 INT_LIST\n")
        f.write("~ D S 1 6 STRING\n")
        f.write("~ D B 1 3 INT\n")
        f.write("~ D X 1 6 STRING\n")
        f.write("~ D G 1 3 INT\n")
        f.write("~ D F 1 6 STRING\n")
        f.write("~ D C 1 8 INT_LIST\n")
        f.write("~ D N 1 11 STRING_LIST\n")

        f.write(f"# G {n_genomes}\n")
        f.write(f"# F {n_genomes}\n")
        f.write(f"# C {n_genomes}\n")
        f.write(f"# N {n_genomes}\n")
        f.write(f"# L {n_mums}\n")
        f.write(f"# P {n_mums}\n")
        f.write(f"# S {n_mums}\n")
        if has_blocks:
            f.write(f"# B {n_mums}\n")
        if has_extras:
            f.write(f"# X {n_mums}\n")

        f.write(f"@ P {n_genomes}\n")
        f.write(f"@ S {n_genomes}\n")
        f.write(f"@ C {max_contigs}\n")
        if path_lens:
            f.write(f"@ F {max(path_lens)}\n")
        if name_line_sums:
            f.write(f"@ N {max(name_line_sums)}\n")
        if extra_lens:
            f.write(f"@ X {max(extra_lens)}\n")

        f.write(f"+ P {n_mums * n_genomes}\n")
        f.write(f"+ S {n_mums * n_genomes}\n")
        f.write(f"+ C {total_contigs}\n")
        f.write(f"+ F {sum(path_lens)}\n")
        f.write(f"+ N {sum(name_line_sums)}\n")
        if extra_lens:
            f.write(f"+ X {sum(extra_lens)}\n")

        # Data: genome catalog, then MUM objects
        for i in range(n_genomes):
            f.write(f"G {genome_totals[i]}\n")
            f.write(f"F {one_string(paths[i])}\n")
            f.write(f"C {one_int_list(contig_lengths[i])}\n")
            f.write(f"N {one_string_list(contig_names[i])}\n")

        for i in range(n_mums):
            starts = mums.starts[i]
            strands = mums.strands[i]
            strand_str = "".join(
                ("+" if bool(strands[j]) else "-") if int(starts[j]) != -1 else "-"
                for j in range(n_genomes)
            )
            f.write(f"L {int(mums.lengths[i])}\n")
            f.write(f"P {one_int_list(starts)}\n")
            f.write(f"S {one_string(strand_str)}\n")
            if has_blocks:
                f.write(f"B {block_ids[i]}\n")
            if has_extras:
                f.write(f"X {one_string(mums.extra_fields[i])}\n")


def main(args):
    try:
        contig_lengths = get_sequence_lengths(args.lengths, multilengths=True)
    except ValueError:
        print(
            "Multi-FASTA / multilengths input required; "
            f"simple lengths format not supported: {args.lengths}",
            file=sys.stderr,
        )
        sys.exit(1)

    paths = get_seq_paths(args.lengths)
    contig_names = get_contig_names(args.lengths)

    mums = MUMdata(args.mumfile, sort=False, verbose=args.verbose)

    command = " ".join([SCRIPT_NAME] + sys.argv[1:])
    try:
        write_1mum(
            args.output,
            mums,
            paths,
            contig_lengths,
            contig_names,
            command=command,
            verbose=args.verbose,
        )
    except ValueError as e:
        print(str(e), file=sys.stderr)
        sys.exit(1)

    if args.verbose:
        print(f"Wrote {args.output}", file=sys.stderr)


if __name__ == "__main__":
    main(parse_arguments())
